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
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <icetray/I3Int.h>
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

class IntervalChargeSum: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;
    double timeWindow_;
    std::string input_prefix_;
    std::string output_prefix_;

    bool process_physics_frames_;
    bool process_daq_frames_;


    bool check_masked_pulses_;
    bool check_raw_pulses_;
    std::string raw_pulses_name_;
    std::string pulses_mask_name_;

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

    public:
    void Geometry(I3FramePtr frame);
    IntervalChargeSum(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Physics(I3FramePtr frame);
};

I3_MODULE(IntervalChargeSum);

IntervalChargeSum::IntervalChargeSum(const I3Context& context) : I3Module(context),
    geo_seen(false), geometry_name_("") {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
        AddParameter("ProcessPhysicsFrames", "Process physics frames", true);
        AddParameter("ProcessDAQFrames", "Process DAQ frames", true);
        AddParameter("TimeWindow", "Time window for charge estimate", 90.0);
        AddParameter("InputPulsesMaskName", "Name of the input pulses mask", std::string(""));
        AddParameter("InputRawPulsesName", "Name of the input raw pulses", std::string(""));
        AddParameter("InputEventPrefix", "Prefix for the inputs", std::string(""));
        AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
}

void IntervalChargeSum::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("ProcessPhysicsFrames", process_physics_frames_);
    GetParameter("ProcessDAQFrames", process_daq_frames_);
    GetParameter("TimeWindow", timeWindow_);
    GetParameter("InputPulsesMaskName", pulses_mask_name_);
    GetParameter("InputRawPulsesName", raw_pulses_name_);
    GetParameter("InputEventPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);

    if(not (process_physics_frames_ or process_daq_frames_)) {
        log_fatal("At least one of ProcessPhysicsFrames or ProcessDAQFrames must be true.");
    }

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

void IntervalChargeSum::Geometry(I3FramePtr frame) {
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

void IntervalChargeSum::DAQ(I3FramePtr frame) {
    if(not process_daq_frames_) {
        PushFrame(frame);
        return;
    }
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

    I3VectorDoubleConstPtr event_start_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "EventStartTimes");
    I3VectorDoubleConstPtr event_end_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "EventEndTimes");

    PMTKeyPulseVector pulse_list;
    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin();
            i != pulses->end(); i++) {
        if(pmt_keys.count(i->first) == 0)
            continue;
        for(CCMRecoPulse const & pulse: i->second) {
            pulse_list.push_back(PMTKeyPulsePair(i->first, pulse));
        }
    }

    std::sort(pulse_list.begin(), pulse_list.end(), [](auto const & t0, auto const & t1){return std::get<1>(t0).GetTime() < std::get<1>(t1).GetTime();});

    I3VectorDoublePtr event_charge = boost::make_shared<I3VectorDouble>(event_start_times->size(), 0.0);

    size_t event_idx = 0;
    for (PMTKeyPulseVector::const_iterator i = pulse_list.begin(); i != pulse_list.end() and event_idx < event_start_times->size(); ++i) {
        if(event_start_times->at(event_idx) > std::get<1>(*i).GetTime()) {
            continue;
        }
        double total_charge = 0.0;
        double start_time = event_start_times->at(event_idx);
        double end_time = event_end_times->at(event_idx);
        for (PMTKeyPulseVector::const_iterator j = i; j != pulse_list.end(); ++j) {
            double time = std::get<1>(*j).GetTime();
            if(time - start_time > timeWindow_ or time > end_time) {
                break;
            }
            total_charge += std::get<1>(*j).GetCharge();
        }
        event_charge->at(event_idx) = total_charge;
        ++event_idx;
    }

    frame->Put(output_prefix_ + "EventCharges", event_charge);

    PushFrame(frame);
}

void IntervalChargeSum::Physics(I3FramePtr frame) {
    if(not process_physics_frames_) {
        PushFrame(frame);
        return;
    }

    if(not process_physics_frames_) {
        PushFrame(frame);
        return;
    }
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    if(process_daq_frames_ and frame->Has(output_prefix_ + "EventCharges") and frame->Has(input_prefix_ + "EventIndex")) {
        int event_index = frame->Get<I3IntConstPtr>(input_prefix_ + "EventIndex")->value;
        frame->Put(output_prefix_ + "EventCharge", boost::make_shared<I3Double>(frame->Get<I3VectorDoubleConstPtr>(output_prefix_ + "EventCharges")->at(event_index)));
        PushFrame(frame);
        return;
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

    PMTKeyPulseVector pulse_list;
    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin();
            i != pulses->end(); i++) {
        if(pmt_keys.count(i->first) == 0)
            continue;
        for(CCMRecoPulse const & pulse: i->second) {
            pulse_list.push_back(PMTKeyPulsePair(i->first, pulse));
        }
    }

    std::sort(pulse_list.begin(), pulse_list.end(), [](auto const & t0, auto const & t1){return std::get<1>(t0).GetTime() < std::get<1>(t1).GetTime();});

    double start_time = frame->Get<I3DoubleConstPtr>(input_prefix_ + "EventStartTime")->value;
    double end_time = frame->Get<I3DoubleConstPtr>(input_prefix_ + "EventEndTime")->value;

    size_t event_idx = 0;
    double total_charge = 0.0;
    for (PMTKeyPulseVector::const_iterator i = pulse_list.begin(); i != pulse_list.end(); ++i) {
        if(start_time > std::get<1>(*i).GetTime()) {
            continue;
        }
        for (PMTKeyPulseVector::const_iterator j = i; j != pulse_list.end(); ++j) {
            double time = std::get<1>(*j).GetTime();
            if(time - start_time > timeWindow_ or time > end_time) {
                break;
            }
            total_charge += std::get<1>(*j).GetCharge();
        }
        break;
    }

    I3DoublePtr event_charge = boost::make_shared<I3Double>(total_charge);

    frame->Put(output_prefix_ + "EventCharge", event_charge);

    PushFrame(frame);
}

