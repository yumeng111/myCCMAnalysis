// This modules reads in the data files and selects only the cosmic triggers to process.
// This sums the waveforms across all PMTs for each cosmic trigger and saves to the frame.
// This will also save the fixed position of the summed waveforms that can be used for timing analyses.

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
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

class CosmicTriggerSumWaveforms: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    bool IsCosmicFrame(I3FramePtr frame);
    std::deque<double> SumWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position);
public:
    CosmicTriggerSumWaveforms(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(CosmicTriggerSumWaveforms);

CosmicTriggerSumWaveforms::CosmicTriggerSumWaveforms(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void CosmicTriggerSumWaveforms::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void CosmicTriggerSumWaveforms::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    bool found_bcm_key = false;
    for(auto const & key_om : geo.pmt_geo) {
        if(key_om.second.omtype == CCMOMGeo::OMType::BeamCurrentMonitor) {
            bcm_key = key_om.first;
            found_bcm_key = true;
        }
    }
    if(not found_bcm_key) {
        log_fatal("CCMGeometry does not contain a channel corresponding to a BeamCurrentMonitor");
    }
    bcm_channel = geo.pmt_channel_map.at(bcm_key);
    cosmic_trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::CosmicTrigger, 1);
    PushFrame(frame);
}

bool CosmicTriggerSumWaveforms::IsCosmicFrame(I3FramePtr frame) {
    if(not frame->Has(nim_pulses_name_)) {
        log_warn(("No key named " + nim_pulses_name_ + " present in frame").c_str());
        return false;
    }
    boost::shared_ptr<NIMLogicPulseSeriesMap const> nim_pulses = frame->Get<boost::shared_ptr<NIMLogicPulseSeriesMap const>>(nim_pulses_name_);
    if(not nim_pulses) {
        log_warn(("No NIMLogicPulseSeriesMap named " + nim_pulses_name_ + " present in frame").c_str());
        return false;
    }
    NIMLogicPulseSeriesMap::const_iterator it = nim_pulses->find(cosmic_trigger_key);
    if(it == nim_pulses->end()) {
        log_warn(("NIMLogicPulseSeriesMap named " + nim_pulses_name_ + " does not contain the CosmicTrigger key").c_str());
        return false;
    }
    if(it->second.size() < 1) {
        // No NIM pulse on the cosmic trigger means we skip this frame
        return false;
    }

    // We found a NIM pulse on the cosmic trigger,
    // therefore this is a cosmic frame
    return true;
}


std::deque<double> CosmicTriggerSumWaveforms::SumWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position){
    
    //reset summed_waveforms_ !!!
    WaveformAccumulator summed_waveforms_;

    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second; 
        double mode = value.baseline; //baseline mode is negative!!!
	uint32_t channel = pmt_channel_map_[key];
        CCMWaveformUInt16 const & waveform = waveforms->at(channel);
        std::vector<short unsigned int> const & samples = waveform.GetWaveform();
        std::vector<double> wf_minus_baseline(samples.size());

        for(size_t i = 0; i < samples.size(); ++i){
            wf_minus_baseline[i] = samples[i] + mode;
        }
	CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
        double pos = frame->Get<NIMLogicPulseSeriesMap>("NIMPulses").at(geo.trigger_copy_map.at(key)).at(0).GetNIMPulseTime() / 2.0 ;
        summed_waveforms_.AddWaveform(wf_minus_baseline, pos);
    }

    std::deque<double> summed_wf = summed_waveforms_.GetSummedWaveform();
    fixed_position = summed_waveforms_.GetFixedPosition();
    return summed_wf;

}


void CosmicTriggerSumWaveforms::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }

    if(not IsCosmicFrame(frame)) {
        PushFrame(frame);
        return;
    }

    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    if(IsCosmicFrame(frame)) {

    // ptr to vector of all waveforms and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    // a shared pointer to store the total charge per pmt for each event
    boost::shared_ptr<I3Vector<int>> summed_wf_reference_time(new I3Vector<int>(1));

    // let's add up all the charge in the detector
    int fixed_position;
    std::deque<double> summed_wf = SumWaveforms(frame, waveforms, baseline_mode, fixed_position);
    boost::shared_ptr<I3Vector<double>> summed_wf_per_frame = boost::make_shared<I3Vector<double>>(summed_wf.begin(), summed_wf.end());
    
    for(size_t wf_it = 0; wf_it < summed_wf.size(); ++wf_it){
       summed_wf_per_frame->operator[](wf_it) = summed_wf[wf_it];
    }
    
    summed_wf_reference_time->operator[](0) = fixed_position;
    
    
    frame->Put("SummedWaveform", summed_wf_per_frame);
    frame->Put("SummedWaveformFixedPosition", summed_wf_reference_time);
    PushFrame(frame);
    }
}

void CosmicTriggerSumWaveforms::Finish() {
}

