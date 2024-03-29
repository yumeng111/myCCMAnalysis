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
#include <dataclasses/I3MapOMKeyMask.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <phys-services/CCMSplitter.h>

#include <string>

class CCMEventSplitter : public I3ConditionalModule, public CCMSplitter {
public:
    CCMEventSplitter(const I3Context& context);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
private:
    std::string input_prefix_, output_prefix_;
    std::string input_pulses_, output_pulses_;

    I3Vector<CCMOMGeo::OMType> trim_pmt_types_ = {CCMOMGeo::OMType::CCM8inUncoated, CCMOMGeo::OMType::CCM8inCoated};
    I3Vector<CCMOMGeo::OMType> keep_pmt_types_ = {CCMOMGeo::OMType::CCM1in};
    I3Vector<CCMOMGeo::OMType> exclude_pmt_types_ = {};

    double pre_padding_;
    double post_padding_;
        
};

I3_MODULE(CCMEventSplitter);

CCMEventSplitter::CCMEventSplitter(const I3Context& context) :
  I3ConditionalModule(context), CCMSplitter(configuration_) {
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
        log_debug("Setting output_prefix_ equal to input_prefix_ with value \"%s\"", input_prefix_.c_str())
    }

    //TODO check that the intersections of the pmt_type lists are empty
}

void CCMEventSplitter::Geometry(I3FramePtr frame) {
    //TODO build lists of PMT keys
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
        boost::function<bool (const CCMPMTKey&, size_t, const CCMRecoPulse&)> predicate = [this, &start_time, &end_time, i] (const CCMPMTKey& pmt_key, size_t pulse_idx, const CCMRecoPulse & pulse) {
            //TODO add the appropriate behavior for the keep, trim, and exclude lists
            if(pmt_keys.count(pmt_key) == 0)
                return false;
            return pulse.GetTime() >= start_time and pulse.GetTime() <= end_time;
        };
        CCMRecoPulseSeriesMapMaskPtr mask = boost::make_shared<CCMRecoPulseSeriesMapMask>(*frame, pulses_name, predicate);
        physics_frame->Put(output_prefix_ + "EventPulses", mask);
        PushFrame(physics_frame);
    }
}

