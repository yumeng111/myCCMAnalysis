/**
 *  $Id$
 *  
 *  Copyright (C) 2011
 *  Nathan Whitehorn <nwhitehorn@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */

#include <icetray/I3ConditionalModule.h>
#include <dataclasses/I3String.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/I3MapOMKeyMask.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <phys-services/CCMSplitter.h>

#include <set>
#include <string>
#include <algorithm>

class CCMEventSplitter : public I3ConditionalModule, public CCMSplitter {
public:
    CCMEventSplitter(const I3Context& context);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
private:
    bool geo_seen;
    std::string geometry_name_;
    std::string input_prefix_, output_prefix_;
    std::string input_pulses_, output_pulses_;

    I3Vector<CCMOMGeo::OMType> trim_pmt_types_ = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    I3Vector<CCMOMGeo::OMType> keep_pmt_types_{{CCMOMGeo::OMType::CCM1in}};
    I3Vector<CCMOMGeo::OMType> exclude_pmt_types_ = {};

    std::set<CCMPMTKey> pmt_keys_trim;
    std::set<CCMPMTKey> pmt_keys_keep;
    std::set<CCMPMTKey> pmt_keys_exclude;

    double pre_padding_;
    double post_padding_;
        
};

I3_MODULE(CCMEventSplitter);

CCMEventSplitter::CCMEventSplitter(const I3Context& context) :
    I3ConditionalModule(context), CCMSplitter(configuration_),
    geo_seen(false), geometry_name_("") {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("PMTTypesToTrim", "List of PMT types to trim the pulse series of", trim_pmt_types_);
    AddParameter("PMTTypesToKeep", "List of PMT types to keep the full pulse series of", keep_pmt_types_);
    AddParameter("PMTTypesToExclude", "List of PMT types to exclude the pulse series of", exclude_pmt_types_);
    AddParameter("PreEventPaddingInNS", "Additional time before the event to avoid trimming", double(0.0));
    AddParameter("PostEventPaddingInNS", "Additional time after the event to avoid trimming", double(0.0));
    AddParameter("InputPrefix", "The prefix name of the events to split on", "");
    AddParameter("OutputPrefix", "The prefix name to output event properties with. Defaults to value of \"InputPrefix\"", "");
	AddParameter("InputPulseSeriesName", "The name of the pulse series to mask "
	    "(optional)", "");
	AddParameter("OutputPulseSeriesMaskName", "The name of the output pulse "
	    "series mask (optional)", "");
	AddParameter("SubEventStreamName",
		     "The name of the SubEvent stream.",
		     configuration_.InstanceName());
}

void CCMEventSplitter::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("PMTTypesToTrim", trim_pmt_types_);
    GetParameter("PMTTypesToKeep", keep_pmt_types_);
    GetParameter("PMTTypesToExclude", exclude_pmt_types_);
    GetParameter("PreEventPaddingInNS", pre_padding_);
    GetParameter("PostEventPaddingInNS", post_padding_);
    GetParameter("InputPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);
	GetParameter("InputPulseSeriesName", input_pulses_);
	GetParameter("OutputPulseSeriesMaskName", output_pulses_);
	GetParameter("SubEventStreamName", sub_event_stream_name_);

    if(input_prefix_.empty())
        log_fatal("You must specify an \"InputPrefix\" parameter for CCMEventSplitter");
    if(output_prefix_.empty()) {
        output_prefix_ = input_prefix_;
        log_debug("Setting output_prefix_ equal to input_prefix_ with value \"%s\"", input_prefix_.c_str());
    }

    std::vector<CCMOMGeo::OMType> pmt_types_intersection;
    std::set_intersection(trim_pmt_types_.begin(), trim_pmt_types_.end(),
                          keep_pmt_types_.begin(), keep_pmt_types_.end(),
                          std::back_inserter(pmt_types_intersection));
    if(pmt_types_intersection.size() > 0)
        log_fatal("The PMT types to trim and keep must be disjoint sets");
    pmt_types_intersection.clear();
    std::set_intersection(trim_pmt_types_.begin(), trim_pmt_types_.end(),
                          exclude_pmt_types_.begin(), exclude_pmt_types_.end(),
                          std::back_inserter(pmt_types_intersection));
    if(pmt_types_intersection.size() > 0)
        log_fatal("The PMT types to trim and exclude must be disjoint sets");
    pmt_types_intersection.clear();
    std::set_intersection(keep_pmt_types_.begin(), keep_pmt_types_.end(),
                          exclude_pmt_types_.begin(), exclude_pmt_types_.end(),
                          std::back_inserter(pmt_types_intersection));
    if(pmt_types_intersection.size() > 0)
        log_fatal("The PMT types to keep and exclude must be disjoint sets");
}

void CCMEventSplitter::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometryConstPtr geo = frame->Get<CCMGeometryConstPtr>(geometry_name_);
    geo_seen = bool(geo);
    pmt_keys_trim.clear();
    pmt_keys_keep.clear();
    pmt_keys_exclude.clear();
    if(geo_seen) {
        std::vector<
            std::tuple<I3Vector<CCMOMGeo::OMType> const &,
            std::set<CCMPMTKey> &>
        > list = {
            std::make_tuple(std::cref(trim_pmt_types_),    std::ref(pmt_keys_trim)),
            std::make_tuple(std::cref(keep_pmt_types_),    std::ref(pmt_keys_keep)),
            std::make_tuple(std::cref(exclude_pmt_types_), std::ref(pmt_keys_exclude))
        };
        for(std::tuple<I3Vector<CCMOMGeo::OMType> const &, std::set<CCMPMTKey> &> & p : list ) {
            I3Vector<CCMOMGeo::OMType> const & pmt_types = std::get<0>(p);
            std::set<CCMOMGeo::OMType> allowed_pmt_types(pmt_types.begin(), pmt_types.end());
            std::set<CCMPMTKey> & pmt_keys = std::get<1>(p);
            for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : geo->pmt_geo) {
                if(allowed_pmt_types.count(p.second.omtype) == 0)
                    continue;
                pmt_keys.insert(p.first);
            }
        }
    }
    PushFrame(frame);
}

void CCMEventSplitter::DAQ(I3FramePtr frame) {
	PushFrame(frame);

    std::string pulses_name;
    if(not input_pulses_.empty()) {
        pulses_name = input_pulses_;
    } else if(frame->Has(input_prefix_ + "EventsPulsesName")) {
        pulses_name = frame->Get<I3String>(input_prefix_ + "EventsPulsesName").value;
    } else {
        log_fatal("Frame must have \"%sEventsPulsesName\" or the parameter \"InputPulseSeriesName\" must be set", input_prefix_.c_str());
    }

    I3VectorDoubleConstPtr event_start_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "EventStartTimes");
    I3VectorDoubleConstPtr event_end_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "EventEndTimes");
    I3VectorDoubleConstPtr max_event_charges = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "MaxEventCharges");
    I3VectorDoubleConstPtr max_event_charge_times = frame->Get<I3VectorDoubleConstPtr>(input_prefix_ + "MaxEventChargeTimes");
    
    size_t n_events = event_start_times->size();

    for(size_t i=0; i<n_events; ++i) {
        I3FramePtr physics_frame = GetNextSubEvent(frame);
        double start_time = event_start_times->at(i);
        double end_time = event_end_times->at(i);
        double const & max_event_charge = max_event_charges->at(i);
        double const & max_event_charge_time = max_event_charge_times->at(i);
        physics_frame->Put(output_prefix_ + "EventStartTime", boost::make_shared<I3Double>(start_time));
        physics_frame->Put(output_prefix_ + "EventEndTime", boost::make_shared<I3Double>(end_time));
        physics_frame->Put(output_prefix_ + "MaxEventCharge", boost::make_shared<I3Double>(max_event_charge));
        physics_frame->Put(output_prefix_ + "MaxEventChargeTime", boost::make_shared<I3Double>(max_event_charge_time));

        start_time -= pre_padding_;
        end_time += post_padding_;
        boost::function<bool (const CCMPMTKey&, size_t, const CCMRecoPulse&)> predicate = [&] (const CCMPMTKey& pmt_key, size_t pulse_idx, const CCMRecoPulse & pulse) {
            if(pmt_keys_exclude.count(pmt_key) > 0)
                return false;
            else if(pmt_keys_keep.count(pmt_key) > 0)
                return true;
            else if(pmt_keys_trim.count(pmt_key) > 0)
                return pulse.GetTime() >= start_time and pulse.GetTime() <= end_time;
            else
                return false;
        };
        CCMRecoPulseSeriesMapMaskPtr mask = boost::make_shared<CCMRecoPulseSeriesMapMask>(*frame, pulses_name, predicate);
        physics_frame->Put(output_prefix_ + "EventPulses", mask);
        PushFrame(physics_frame);
    }
}

