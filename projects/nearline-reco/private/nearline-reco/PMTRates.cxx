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

class PMTRates: public I3ConditionalModule {
    bool geo_seen;
    std::string geometry_name_;
    CCMGeometryConstPtr geo;

    I3Vector<CCMOMGeo::OMType> pmt_types;
    std::set<CCMPMTKey> pmt_keys;

    std::string pulses_;

    double min_time_;
    double max_time_;

    std::string output_prefix_;

    size_t num_frames_to_average_;
    size_t process_every_n_frames_;

    size_t frame_mod = 0;

    std::map<CCMPMTKey, std::deque<double>> running_charge_list_;

    std::string postfix_;

    public:
    void Geometry(I3FramePtr frame);
    PMTRates(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
};

I3_MODULE(PMTRates);

PMTRates::PMTRates(const I3Context& context) : I3ConditionalModule(context),
    geo_seen(false), geometry_name_(""), pmt_types(I3Vector<CCMOMGeo::OMType>{CCMOMGeo::OMType::CCM8inCoated, CCMOMGeo::OMType::CCM8inUncoated}) {
    I3Vector<double> default_time_windows;
    default_time_windows.push_back(90.0);
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
    AddParameter("Pulses", "Name of the input raw pulses", std::string("WavedeformPulses"));
    AddParameter("MinTime", "Minimum time from which to sum pulses", -std::numeric_limits<double>::infinity());
    AddParameter("MaxTime", "Maximum time from which to sum pulses", std::numeric_limits<double>::infinity());
    AddParameter("OutputPrefix", "Prefix for the outputs", std::string(""));
    AddParameter("NumFramesToAverage", "Number of frames to average", size_t(1));
    AddParameter("ProcessEveryNFrames", "Process only one frame out of this many frames", size_t(1));
}

void PMTRates::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("Pulses", pulses_);
    GetParameter("MinTime", min_time_);
    GetParameter("MaxTime", max_time_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("NumFramesToAverage", num_frames_to_average_);
    GetParameter("ProcessEveryNFrames", process_every_n_frames_);

    std::string min_s;
    if(std::isinf(min_time_)) {
        min_s = "-Inf";
    } else {
        min_s = std::to_string(int(min_time_));
    }

    std::string max_s;
    if(std::isinf(max_time_)) {
        max_s = "Inf";
    } else {
        max_s = std::to_string(int(max_time_));
    }

    postfix_ = min_s + "to" + max_s;
}

void PMTRates::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    pmt_keys.clear();
    running_charge_list_.clear();
    frame_mod = 0;
    if(geo_seen) {
        std::set<CCMOMGeo::OMType> allowed_pmt_types(pmt_types.begin(), pmt_types.end());
        for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
            if(allowed_pmt_types.count(p.second.omtype) == 0)
                continue;
            pmt_keys.insert(p.first);
            running_charge_list_.insert(std::make_pair(p.first, std::deque<double>()));
        }
    }
    PushFrame(frame);
}

void PMTRates::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen yet.");
    }

    frame_mod += 1;
    if(frame_mod % process_every_n_frames_ != 0) {
        PushFrame(frame);
        return;
    } else {
        frame_mod = 0;
    }

    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_);

    if(not pulses) {
        std::stringstream ss;
        ss << "Could not find ";
        ss << pulses_;
        ss << " in the DAQ frame.";
        log_fatal("%s", ss.str().c_str());
        PushFrame(frame);
        return;
    }

    I3MapPMTKeyDoublePtr total_charges = boost::make_shared<I3MapPMTKeyDouble>();

    for(CCMRecoPulseSeriesMap::const_iterator p = pulses->begin();
            p != pulses->end(); p++) {
        if(pmt_keys.count(p->first) == 0)
            continue;
        double charge = 0.0;
        for(CCMRecoPulseSeries::const_iterator i = p->second.begin(); i != p->second.end(); ++i) {
            double const & time = i->GetTime();
            if(min_time_ > time)
                continue;
            if(max_time_ < time)
                break;
            charge += i->GetCharge();
        }
        std::deque<double> & list = running_charge_list_[p->first];
        list.push_back(charge);
        if(list.size() > num_frames_to_average_) {
            list.pop_front();
        }
        double average_charge = 0.0;
        for(double const & x : list) {
            average_charge += x;
        }
        average_charge /= list.size();
        total_charges->insert(std::make_pair(p->first, average_charge));
    }
    frame->Put(output_prefix_ + "PMTRates" + postfix_, total_charges);

    PushFrame(frame);
}

