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

class IntervalChargeSumQ: public I3ConditionalModule {
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

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

    public:
    void Geometry(I3FramePtr frame);
    IntervalChargeSumQ(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
};

I3_MODULE(IntervalChargeSumQ);

IntervalChargeSumQ::IntervalChargeSumQ(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_("") {
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

void IntervalChargeSumQ::Configure() {
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

void IntervalChargeSumQ::Geometry(I3FramePtr frame) {
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

void IntervalChargeSumQ::DAQ(I3FramePtr frame) {
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

    std::vector<I3VectorDoublePtr> event_charge_list;
    for (size_t i = 0; i < time_windows_.size(); ++i) {
        I3VectorDoublePtr event_charge = boost::make_shared<I3VectorDouble>(event_start_times->size(), 0.0);
        event_charge_list.push_back(event_charge);
    }

    size_t event_idx = 0;
    for (PMTKeyPulseVector::const_iterator i = pulse_list.begin(); i != pulse_list.end() and event_idx < event_start_times->size(); ++i) {
        if(event_start_times->at(event_idx) > std::get<1>(*i).GetTime()) {
            continue;
        }
        size_t time_bin = 0;
        double total_charge = 0.0;
        double start_time = event_start_times->at(event_idx);
        double end_time = event_end_times->at(event_idx);
        for (PMTKeyPulseVector::const_iterator j = i; j != pulse_list.end(); ++j) {
            double time = std::get<1>(*j).GetTime();
            if(time - start_time > time_windows_[time_bin] or time > end_time) {
                event_charge_list[time_bin]->at(event_idx) = total_charge;
                ++time_bin;
                if(time_bin == time_windows_.size()) {
                    break;
                }
            }
            total_charge += std::get<1>(*j).GetCharge();
        }
        if(time_bin < time_windows_.size())
            event_charge_list[time_bin]->at(event_idx) = total_charge;
        time_bin = 0;
        ++event_idx;
    }

    for (size_t i = 0; i < time_windows_.size(); ++i) {
        std::string output_key = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_[i])) + "NS";
        frame->Put(output_key, event_charge_list[i]);
    }

    PushFrame(frame);
}

class IntervalChargeSumP: public I3ConditionalModule {
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

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

    public:
    void Geometry(I3FramePtr frame);
    IntervalChargeSumP(const I3Context&);
    void Configure();
    void Physics(I3FramePtr frame);
};

I3_MODULE(IntervalChargeSumP);

IntervalChargeSumP::IntervalChargeSumP(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_("") {
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

void IntervalChargeSumP::Configure() {
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

void IntervalChargeSumP::Geometry(I3FramePtr frame) {
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

void IntervalChargeSumP::Physics(I3FramePtr frame) {
    std::vector<std::string> keys = frame->keys();

    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    std::string input_key = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_.at(0))) + "NS";

    size_t event_index = frame->Get<CCMEventHeaderConstPtr>("CCMEventHeader")->GetSubEventID();

    if(frame->Has(input_key)) {
        for(size_t i = 0; i < time_windows_.size(); ++i) {
            std::string input_key = output_prefix_ + "EventCharges" + std::to_string(int(time_windows_[i])) + "NS";
            std::string output_key = output_prefix_ + "EventCharge" + std::to_string(int(time_windows_[i])) + "NS";
            if(not frame->Has(input_key)) {
                log_fatal("Could not find \"%s\" in the frame.", input_key.c_str());
            }
            frame->Put(output_key, boost::make_shared<I3Double>(frame->Get<I3VectorDoubleConstPtr>(input_key)->at(event_index)));
        }
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

    std::vector<double> total_charge_list(time_windows_.size(), 0.0);

    size_t time_bin = 0;
    size_t event_idx = 0;
    double total_charge = 0.0;
    for(PMTKeyPulseVector::const_iterator i = pulse_list.begin(); i != pulse_list.end(); ++i) {
        if(start_time > std::get<1>(*i).GetTime()) {
            continue;
        }
        for (PMTKeyPulseVector::const_iterator j = i; j != pulse_list.end(); ++j) {
            double time = std::get<1>(*j).GetTime();
            if(time > end_time)
                break;
            if(time - start_time > time_windows_[time_bin]) {
                ++time_bin;
                total_charge_list[time_bin] = total_charge;
                if(time_bin == time_windows_.size()) {
                    break;
                }
            }
            total_charge += std::get<1>(*j).GetCharge();
        }
        if(time_bin < time_windows_.size())
            total_charge_list[time_bin] = total_charge;
        break;
    }

    for(size_t i = 0; i < time_windows_.size(); ++i) {
        I3DoublePtr charge = boost::make_shared<I3Double>(total_charge_list[i]);
        std::string key = output_prefix_ + "EventCharge" + std::to_string(int(time_windows_[i])) + "NS";
        frame->Put(key, charge);
    }

    PushFrame(frame);
}
