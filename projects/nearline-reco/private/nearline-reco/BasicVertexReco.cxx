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
#include <dataclasses/I3Double.h>
#include <dataclasses/I3UInt32.h>
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

class BasicVertexReco: public I3ConditionalModule {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    std::string input_prefix_;
    std::string output_prefix_;

    double start_time_;
    double end_time_;

    bool check_masked_pulses_;
    bool check_raw_pulses_;
    std::string raw_pulses_name_;
    std::string pulses_mask_name_;
    bool charge_sq_weighting_;

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

    public:
    void Geometry(I3FramePtr frame);
    BasicVertexReco(const I3Context&);
    void Configure();
    void Physics(I3FramePtr frame);
};

I3_MODULE(BasicVertexReco);

BasicVertexReco::BasicVertexReco(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_("") {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("StartTime", "Time from event start to begin using pulses for reconstruction", double(8.0));
    AddParameter("EndTime", "Time from event start to stop using pulses for reconstruction", double(20.0));
    AddParameter("InputPulsesMaskName", "Name of the input pulses mask", std::string(""));
    AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
    AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
    AddParameter("ChargeSquaredWeighting", "Weight by charge squared", bool(false));
}


void BasicVertexReco::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("StartTime", start_time_);
    GetParameter("EndTime", end_time_);
    GetParameter("InputPulsesMaskName", pulses_mask_name_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("ChargeSquareWeighting", charge_sq_weighting_);

    check_masked_pulses_ = false;
    check_raw_pulses_ = (raw_pulses_name_ != "");
    check_masked_pulses_ = (pulses_mask_name_ != "");
    if(not check_masked_pulses_) {
        pulses_mask_name_ = input_prefix_ + "EventPulses";
        if(not check_raw_pulses_) {
            check_masked_pulses_ = true;
        }
    }
}

void BasicVertexReco::Geometry(I3FramePtr frame) {
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

void BasicVertexReco::Physics(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    bool raw_pulses = false;
    CCMRecoPulseSeriesMapConstPtr pulses;
    if(check_raw_pulses_) {
        pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(raw_pulses_name_);
        raw_pulses = bool(pulses);
    }
    if(check_masked_pulses_) {
        CCMRecoPulseSeriesMapMaskConstPtr mask = frame->Get<CCMRecoPulseSeriesMapMaskConstPtr>(pulses_mask_name_);
        if(mask) {
            raw_pulses = false;
            pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_mask_name_);
        }
    }

    if(not pulses) {
        std::stringstream ss;
        ss << "Could not find ";
        if(check_raw_pulses_) {
            ss << raw_pulses_name_;
            if(check_masked_pulses_) {
                ss << " or ";
            }
        }
        if(check_masked_pulses_) {
            ss << pulses_mask_name_;
        }
        ss << " in the DAQ frame.";
        log_fatal("%s", ss.str().c_str());
        PushFrame(frame);
        return;
    }

    I3DoubleConstPtr start_time_ptr = frame->Get<I3DoubleConstPtr>(input_prefix_ + "EventStartTime");
    I3DoubleConstPtr end_time_ptr = frame->Get<I3DoubleConstPtr>(input_prefix_ + "EventEndTime");

    if(not start_time_ptr) {
        log_fatal("Count not find event start time \"%s\" in the frame.", (input_prefix_ + "EventStartTime").c_str());
    }

    if(not end_time_ptr) {
        log_fatal("Count not find event end time \"%s\" in the frame.", (input_prefix_ + "EventEndTime").c_str());
    }

    double start_time = start_time_ptr->value;
    double end_time = end_time_ptr->value;

    I3PositionPtr reco_vertex = boost::make_shared<I3Position>(0, 0, 0);
    double total_weight = 0.0;

    double min_time = start_time + start_time_;
    double max_time = std::min(end_time, start_time + end_time_);

    for(CCMRecoPulseSeriesMap::const_iterator p = pulses->begin();
            p != pulses->end(); p++) {
        if(pmt_keys.count(p->first) == 0)
            continue;
        I3Position const & pos = geo->pmt_geo.at(p->first).position;
        double total_charge = 0.0;
        for(CCMRecoPulseSeries::const_iterator i = p->second.begin(); i != p->second.end(); ++i) {
            double const & time = i->GetTime();
            if(min_time > time)
                continue;
            if(max_time < time)
                break;
            total_charge += i->GetCharge();
        }
        if(charge_sq_weighting_)
            total_charge = total_charge*total_charge;

        total_weight += total_charge;
        (*reco_vertex) += pos * total_charge;
    }

    (*reco_vertex) /= total_weight;

    frame->Put(output_prefix_ + "RecoVertex" + std::to_string(int(start_time_)) + "to" + std::to_string(int(end_time_)), reco_vertex);

    PushFrame(frame);
}

