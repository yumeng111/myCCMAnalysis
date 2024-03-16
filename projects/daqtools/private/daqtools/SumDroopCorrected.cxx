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
#include "daqtools/WaveformAccumulator.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

class SumDroopCorrected: public I3Module {
    bool geo_seen;
    WaveformAccumulator summed_waveforms_;
    std::string geometry_name_;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    std::deque<double> SumWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position);
    void ProcessWaveform(CCMWaveformUInt16 const & waveform, double const & mode, size_t const & michel_start_index, size_t const & michel_end_index, double & total_charge);
public:
    SumDroopCorrected(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(SumDroopCorrected);

SumDroopCorrected::SumDroopCorrected(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void SumDroopCorrected::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void SumDroopCorrected::Geometry(I3FramePtr frame) {
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


std::deque<double> SumDroopCorrected::SumWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position){
    
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

	for(std::pair<CCMPMTKey const, CCMTriggerKey> const & it : geo.trigger_copy_map){
           CCMPMTKey key = it.first;
           CCMTriggerKey value = it.second;
	   std::cout << "key = " << key << std::endl;
	   std::cout << "trigger value = " << value << std::endl;

	}

        double pos = frame->Get<NIMLogicPulseSeriesMap>("NIMPulses").at(geo.trigger_copy_map.at(key)).at(0).GetNIMPulseTime() / 2.0 ;
        summed_waveforms_.AddWaveform(wf_minus_baseline, pos);
    }

    std::deque<double> summed_wf = summed_waveforms_.GetSummedWaveform();
    fixed_position = summed_waveforms_.GetFixedPosition();
    std::cout << "fixed position = " << summed_waveforms_.GetFixedPosition() << std::endl;
    return summed_wf;

}


void SumDroopCorrected::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    size_t size = waveforms->size();

    // a vector storing the droop-corrected wfs for each channel
    boost::shared_ptr<I3Vector<I3Vector<double>>> droop_correced_wf(new I3Vector<I3Vector<double>>(size));

    // loop over each pmt
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second;
        double mode = value.baseline; // baseline mode is negative fyi
        uint32_t channel = pmt_channel_map_[key];

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);

        I3Vector<double> droop_corrected_wf_per_channel(waveform.GetWaveform().size());

        ProcessWaveform(waveform, mode, droop_corrected_wf_per_channel);
        droop_correced_wf->operator[](channel) = droop_corrected_wf_per_channel;
    }

    frame->Put("DroopCorrectedWf", droop_correced_wf);
    std::cout << "finished droop correcting wf" << std::endl;
    PushFrame(frame);
}

void SumDroopCorrected::ProcessWaveform(CCMWaveformUInt16 const & waveform, double const & mode, I3Vector<double> & droop_corrected_wf_per_channel){

    // get the vector of samples from the CCMWaveform object;
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();

    if (samples.size() == 0) {
	    return;
    }
    // let's subtract off the baseline and invert our waveform before passing it to the droop correction function
    std::vector<double> inv_wf (samples.size());

    for (size_t i = 0; i < samples.size(); ++i){
        inv_wf[i] = -1*(double(samples[i]) + mode);
    }

    // let's call our droop correction function and then we're done!
    CorrectDroop(inv_wf, mode, droop_corrected_wf_per_channel);
}


void SumDroopCorrected::Finish() {
}


