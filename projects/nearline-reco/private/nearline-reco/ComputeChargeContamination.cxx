#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <limits>
#include <string>
#include <fstream>
#include <numeric>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3ConditionalModule.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <icetray/I3Int.h>
#include <dataclasses/physics/CCMEventHeader.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3String.h>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

typedef std::tuple<CCMPMTKey, CCMRecoPulse> PMTKeyPulsePair;
typedef std::vector<PMTKeyPulsePair> PMTKeyPulseVector;

class ComputeChargeContamination: public I3ConditionalModule {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;

    double tau_s_;
    double tau_t_;
    double delta_t_;

    double beta_s_;
    double beta_t_;

    std::string event_input_prefix_;
    std::string unfolding_input_prefix_;
    std::string output_prefix_;

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

public:
    void Geometry(I3FramePtr frame);
    ComputeChargeContamination(const I3Context&);
    void Configure();
    void Physics(I3FramePtr frame);
};

I3_MODULE(ComputeChargeContamination);

ComputeChargeContamination::ComputeChargeContamination(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_("") {
    I3Vector<double> default_time_windows;
    default_time_windows.push_back(90.0);
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("TauSinglet", "Time constant for singlet light", 8.13);
    AddParameter("TauTriplet", "Time constant for triplet light", 743);
    AddParameter("BinWidth", "Width of the time bins", 2.0);
    AddParameter("InputEventPrefix", "Prefix for the event inputs", std::string(""));
    AddParameter("InputUnfoldingPrefix", "Prefix for the charge unfolding inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
}

void ComputeChargeContamination::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("TauSinglet", tau_s_);
    GetParameter("TauTriplet", tau_t_);
    GetParameter("BinWidth", delta_t_);
    GetParameter("InputEventPrefix", event_input_prefix_);
    GetParameter("InputUnfoldingPrefix", unfolding_input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);

    beta_s_ = exp(-delta_t_ / tau_s_);
    beta_t_ = exp(-delta_t_ / tau_t_);
}

void ComputeChargeContamination::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    pmt_keys.clear();
    if(geo_seen) {
        std::set<CCMOMGeo::OMType> allowed_pmt_types(pmt_types.begin(), pmt_types.end());
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            if(allowed_pmt_types.count(p.second.omtype) == 0)
                continue;
            pmt_keys.insert(p.first);
        }
    }
    PushFrame(frame);
}

void ComputeChargeContamination::Physics(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    I3DoubleConstPtr min_time = frame->Get<I3DoubleConstPtr>(unfolding_input_prefix_ + "ChargeUnfoldingMinTime");
    I3VectorDoubleConstPtr singlet = frame->Get<I3VectorDoubleConstPtr>(unfolding_input_prefix_ + "ChargeUnfoldingTotalSinglet");
    I3VectorDoubleConstPtr triplet = frame->Get<I3VectorDoubleConstPtr>(unfolding_input_prefix_ + "ChargeUnfoldingTotalTriplet");
    I3DoubleConstPtr start_time = frame->Get<I3DoubleConstPtr>(event_input_prefix_ + "EventStartTime");
    I3DoubleConstPtr end_time = frame->Get<I3DoubleConstPtr>(event_input_prefix_ + "EventEndTime");

    if(not min_time)
        log_fatal(("I3Double " + unfolding_input_prefix_ + "ChargeUnfoldingMinTime" + " not found in the frame.").c_str());
    if(not singlet)
        log_fatal(("I3VectorDouble " + unfolding_input_prefix_ + "ChargeUnfoldingTotalSinglet" + " not found in the frame.").c_str());
    if(not triplet)
        log_fatal(("I3VectorDouble " + unfolding_input_prefix_ + "ChargeUnfoldingTotalTriplet" + " not found in the frame.").c_str());
    if(not start_time)
        log_fatal(("I3Double " + event_input_prefix_ + "EventStartTime" + " not found in the frame.").c_str());
    if(not end_time)
        log_fatal(("I3Double " + event_input_prefix_ + "EventEndTime" + " not found in the frame.").c_str());

    double total_charge = 0.0;
    I3DoubleConstPtr event_charge = frame->Get<I3DoubleConstPtr>(event_input_prefix_ + "EventTotalCharge");
    if(event_charge) {
        total_charge = event_charge->value;
    } else {
        CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(event_input_prefix_ + "EventPulses");
        if(not pulses)
            log_fatal(("CCMRecoPulseSeriesMapMask " + event_input_prefix_ + "EventPulses" + " not found in the frame.").c_str());

        // Compute total charge
        for(CCMRecoPulseSeriesMap::const_iterator i = pulses->begin(); i != pulses->end(); i++) {
            if(pmt_keys.count(i->first) == 0)
                continue;
            for(CCMRecoPulse const & pulse: i->second) {
                if(pulse.GetTime() < start_time->value)
                    continue;
                if(pulse.GetTime() > end_time->value)
                    break;
                total_charge += pulse.GetCharge();
            }
        }
    }

    if(singlet->size() != triplet->size())
        log_fatal("Singlet and triplet vectors have different sizes.");

    if(singlet->size() == 0)
        log_fatal("Singlet and triplet vectors are empty.");

    // Compute contamination
    size_t first_bin = std::max(std::floor((start_time->value - min_time->value) / delta_t_), 0.0);
    size_t last_bin = std::max(std::floor((end_time->value - min_time->value) / delta_t_), 0.0);
    first_bin = std::min(first_bin + 1, singlet->size()) - 1;
    last_bin = std::min(last_bin + 1, singlet->size()) - 1;
    size_t n_bins = last_bin - first_bin + 1;

    double starting_singlet = 0.0;
    double starting_triplet = 0.0;
    if(first_bin == 0) {
        starting_singlet = singlet->at(0);
        starting_triplet = triplet->at(0);
        starting_singlet /= beta_s_;
        starting_triplet /= beta_t_;
    } else {
        starting_singlet = singlet->at(first_bin - 1);
        starting_triplet = triplet->at(first_bin - 1);
    }

    double total_singlet = starting_singlet * (1.0 - std::pow(beta_s_, n_bins)) / (1.0 / beta_s_ - 1.0);
    double total_triplet = starting_triplet * (1.0 - std::pow(beta_t_, n_bins)) / (1.0 / beta_t_ - 1.0);

    I3DoublePtr total_contamination = boost::make_shared<I3Double>(total_singlet + total_triplet);
    I3DoublePtr total_contamination_fraction = boost::make_shared<I3Double>((total_singlet + total_triplet) / total_charge);

    frame->Put(output_prefix_ + "ChargeContamination", total_contamination);
    frame->Put(output_prefix_ + "ChargeContaminationFraction", total_contamination_fraction);

    if(not event_charge) {
        I3DoublePtr total_charge_ptr = boost::make_shared<I3Double>(total_charge);
        frame->Put(event_input_prefix_ + "EventTotalCharge", total_charge_ptr);
    }

    PushFrame(frame);
}

