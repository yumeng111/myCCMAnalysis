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

class LargestPMTFraction: public I3ConditionalModule {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    I3Vector<double> time_windows_;
    std::string input_prefix_;
    std::string output_prefix_;

    bool check_masked_pulses_;
    bool check_raw_pulses_;
    std::string raw_pulses_name_;
    std::string pulses_mask_name_;

    I3Vector<CCMOMGeo::OMType> pmt_types;
    std::set<CCMPMTKey> pmt_keys;

    public:
    void Geometry(I3FramePtr frame);
    LargestPMTFraction(const I3Context&);
    void Configure();
    void Physics(I3FramePtr frame);
};

I3_MODULE(LargestPMTFraction);

LargestPMTFraction::LargestPMTFraction(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_(""), pmt_types(I3Vector<CCMOMGeo::OMType>{CCMOMGeo::OMType::CCM8inCoated, CCMOMGeo::OMType::CCM8inUncoated}) {
    I3Vector<double> default_time_windows;
    default_time_windows.push_back(90.0);
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("TimeWindows", "Time window for charge estimate", default_time_windows);
    AddParameter("InputPulsesMaskName", "Name of the input pulses mask", std::string(""));
    AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
    AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
}

void LargestPMTFraction::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("TimeWindows", time_windows_);
    GetParameter("InputPulsesMaskName", pulses_mask_name_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);

    std::sort(time_windows_.begin(), time_windows_.end());

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

void LargestPMTFraction::Geometry(I3FramePtr frame) {
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

void LargestPMTFraction::Physics(I3FramePtr frame) {
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

    std::vector<double> fractions(time_windows_.size(), 0);
    std::vector<double> max_charge(time_windows_.size(), 0);
    std::vector<double> total_charge(time_windows_.size(), 0);
    for(CCMRecoPulseSeriesMap::const_iterator p = pulses->begin();
            p != pulses->end(); p++) {
        if(pmt_keys.count(p->first) == 0)
            continue;
        double charge = 0.0;
        size_t time_bin = 0;
        for(CCMRecoPulseSeries::const_iterator i = p->second.begin(); i != p->second.end(); ++i) {
            double const & time = i->GetTime();
            if(start_time > time)
                continue;
            if(end_time < time)
                break;

            charge += i->GetCharge();

            double dt = time - start_time;

            if(time_windows_[time_bin] < dt) {
                max_charge[time_bin] = std::max(max_charge[time_bin], charge);
                total_charge[time_bin] += charge;
                ++time_bin;
                if(time_bin == time_windows_.size())
                    break;
            }
        }
        if(time_bin < time_windows_.size()) {
            max_charge[time_bin] = std::max(max_charge[time_bin], charge);
            total_charge[time_bin] += charge;
        }
    }
    for(size_t time_bin = 0; time_bin < time_windows_.size(); ++time_bin) {
        if(total_charge[time_bin] > 0.0) {
            fractions[time_bin] = max_charge[time_bin] / total_charge[time_bin];
        }
        frame->Put(output_prefix_ + "LargestPMTFraction" + std::to_string(int(time_windows_[time_bin])) + "NS", boost::make_shared<I3Double>(fractions[time_bin]));
    }

    PushFrame(frame);
}

