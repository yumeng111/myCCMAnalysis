// this module counts the number of cosmic triggers in a given dataset

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
#include <string>
#include <fstream>
#include <iostream>
#include <limits>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

typedef std::tuple<CCMPMTKey, CCMRecoPulse> PMTKeyPulsePair;
typedef std::vector<PMTKeyPulsePair> PMTKeyPulseVector;

class EventFinder: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string nim_pulses_name_;
    CCMGeometryConstPtr geo;
    double timeWindow_;
    double event_charge_threshold_;
    bool allow_overlapping_events_;
    std::string output_prefix_;
    std::string pulses_;

    I3Vector<CCMOMGeo::OMType> pmt_types = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    std::set<CCMPMTKey> pmt_keys;

    public:
    void Geometry(I3FramePtr frame);
    EventFinder(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(EventFinder);

EventFinder::EventFinder(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
        AddParameter("TimeWindow", "Size of sliding time window to examine.",
                2 * I3Units::ns);
        AddParameter("EventChargeThreshold", "Charge threshold in window to define an event", double(5.0));
        AddParameter("Pulses", "Name of pulse series to use", "OfflinePulses");
        AddParameter("PMTTypes", "PMT types to use for event finding", pmt_types);
        AddParameter("AllowOverlappingEvents", "False -> merge overlapping event windows. True -> allow overlapping event windows", bool(false));
        AddParameter("Ouptut", "Prefix for the outputs", std::string(""));
    }


void EventFinder::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("TimeWindow", timeWindow_);
    GetParameter("EventChargeThreshold", event_charge_threshold_);
    GetParameter("Pulses", pulses_);
    GetParameter("PMTTypes", pmt_types);
    GetParameter("AllowOverlappingEvents", allow_overlapping_events_);
    GetParameter("Ouptut", output_prefix_);
}


void EventFinder::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
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

void EventFinder::DAQ(I3FramePtr frame) {
    if(not frame->Has(nim_pulses_name_)) {
        log_warn(("No key named " + nim_pulses_name_ + " present in frame").c_str());
    }
    boost::shared_ptr<NIMLogicPulseSeriesMap const> nim_pulses = frame->Get<boost::shared_ptr<NIMLogicPulseSeriesMap const>>(nim_pulses_name_);
    CCMRecoPulseSeriesMapConstPtr pulses =
        frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulses_);
    if (!pulses) {
        PushFrame(frame);
        return;
    }
    if (!geo) {
        log_fatal("No Geometry frame seen yet.");
    }
    PMTKeyPulseVector pulse_list;

    for (CCMRecoPulseSeriesMap::const_iterator i = pulses->begin();
            i != pulses->end(); i++) {
        if(pmt_keys.count(i->first) == 0)
            continue;
        double nim_pulse_time = nim_pulses->at(geo->trigger_copy_map.at(i->first)).at(0).GetNIMPulseTime();
        for(CCMRecoPulse const & pulse: i->second) {
            pulse_list.push_back(PMTKeyPulsePair(i->first, pulse));
            std::get<1>(pulse_list.back()).SetTime(std::get<1>(pulse_list.back()).GetTime() - nim_pulse_time);
        }
    }

    std::sort(pulse_list.begin(), pulse_list.end(), [](auto const & t0, auto const & t1){return std::get<1>(t0).GetTime() < std::get<1>(t1).GetTime();});

    double event_start_time = NAN;
    double event_end_time = NAN;
    double max_event_charge = 0.0;
    double max_event_charge_time = NAN;
    std::vector<std::tuple<double, double, double, double>> events;
    bool have_event = false;
    double total_charge = 0.0;
    for (PMTKeyPulseVector::const_iterator i = pulse_list.begin(), j = pulse_list.begin(); j != pulse_list.end(); j++) {
        total_charge += std::get<1>(*j).GetCharge();
        while(std::get<1>(*j).GetTime() - std::get<1>(*i).GetTime() > timeWindow_) {
            total_charge -= std::get<1>(*i).GetCharge();
            i++;
        }

        if(total_charge > event_charge_threshold_) {
            if(have_event) {
                if(total_charge > max_event_charge) {
                    max_event_charge = total_charge;
                    max_event_charge_time = std::get<1>(*i).GetTime();
                }
            } else {
                max_event_charge = total_charge;
                max_event_charge_time = std::get<1>(*i).GetTime();
                event_start_time = std::get<1>(*i).GetTime();
                have_event = true;
            }
        } else {
            if(have_event) {
                event_end_time = std::get<1>(*j).GetTime();

                if((not allow_overlapping_events_) and events.size() > 0 and std::get<1>(events.back()) > event_start_time) {
                    if(max_event_charge > std::get<2>(events.back())) {
                        max_event_charge = std::get<2>(events.back());
                        max_event_charge_time = std::get<3>(events.back());
                    }
                    events.back() = std::make_tuple(std::get<0>(events.back()), event_end_time, max_event_charge, max_event_charge_time);
                } else {
                    events.emplace_back(event_start_time, event_end_time, max_event_charge, max_event_charge_time);
                }
                have_event = false;
                max_event_charge = 0.0;
            } else {
            }
        }
    }

    boost::shared_ptr<I3VectorDouble> event_start_times = boost::make_shared<I3VectorDouble>(events.size());
    boost::shared_ptr<I3VectorDouble> event_end_times = boost::make_shared<I3VectorDouble>(events.size());
    boost::shared_ptr<I3VectorDouble> max_event_charges = boost::make_shared<I3VectorDouble>(events.size());
    boost::shared_ptr<I3VectorDouble> max_event_charge_times = boost::make_shared<I3VectorDouble>(events.size());
    for(size_t i=0; i< events.size(); ++i) {
        event_start_times->at(i) = std::get<0>(events[i]);
        event_end_times->at(i) = std::get<1>(events[i]);
        max_event_charges->at(i) = std::get<2>(events[i]);
        max_event_charge_times->at(i) = std::get<3>(events[i]);
    }

    frame->Put(output_prefix_ + "EventStartTimes", event_start_times);
    frame->Put(output_prefix_ + "EventEndTimes", event_end_times);
    frame->Put(output_prefix_ + "MaxEventCharges", max_event_charges);
    frame->Put(output_prefix_ + "MaxEventChargeTimes", max_event_charge_times);

    PushFrame(frame);
}

void EventFinder::Finish() {
}

