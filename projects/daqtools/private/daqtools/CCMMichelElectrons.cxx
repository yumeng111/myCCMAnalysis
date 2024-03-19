// this module identifies regions where we want to look for michel electrons
// reconstructs the vertex of the event by calling VertexReconstruction.h
// if the event is near the center of the detector, we sum up charge in each pmt and save to frame
// we ultimately want to save cherenkov light per PMT to frame, but don't have pulse finding...
// so just saving total charge of event for now as a place holder

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
//#include "daqtools/VertexReconstruction.h"

class CCMMichelElectrons: public I3Module {
    WaveformAccumulator summed_waveforms_;
    std::string geometry_name_;
    bool geo_seen;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    bool IsCosmicFrame(I3FramePtr frame);
    std::deque<double> SumWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position);
    void FindMicheleRegion(std::deque<double> const & summed_wf, size_t & michel_start_index, size_t & michel_end_index);
    void ProcessWaveform(CCMWaveformUInt16 const & waveform, double const & mode, size_t const & michel_start_index, size_t const & michel_end_index, double & total_charge);
public:
    CCMMichelElectrons(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(CCMMichelElectrons);

CCMMichelElectrons::CCMMichelElectrons(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void CCMMichelElectrons::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void CCMMichelElectrons::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
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

bool CCMMichelElectrons::IsCosmicFrame(I3FramePtr frame) {
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

std::deque<double> CCMMichelElectrons::SumWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position){
    
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


void CCMMichelElectrons::FindMicheleRegion(std::deque<double> const & summed_wf, size_t & michel_start_index, size_t & michel_end_index){
    // so this frame loops over the summed wf and sees if there's a peak of michele singlet light ~ 2microseconds after the trigger
    double diff_between_bins_threshold = 0.5;
    double baseline_threshold = 3;

    // let's iterate over the summed_wf to find the maximum value aka our muon triplet peak 
    double muon_peak_value = summed_wf[0];
    size_t muon_peak_value_position = 0;

    for(size_t wf_it = 0; wf_it < summed_wf.size(); ++wf_it){
       if(summed_wf[wf_it] < muon_peak_value){    //CHECK IF WF IS INVERTED OR NOT
          muon_peak_value = summed_wf[wf_it];
	  muon_peak_value_position = wf_it;
       }
    
    }

    // now let's loop over the waveform and see what the max value is after ~2microseconds or more after the max_value_position
    size_t two_usec_after_muon_peak = std::min(summed_wf.size(), muon_peak_value_position+1000); // 1000 bins = 2000 ns = 2 microseconds
    double michel_peak_value = summed_wf[two_usec_after_muon_peak];
    size_t michel_peak_value_position = two_usec_after_muon_peak;
    size_t possible_michel_end_index = 0;

    for(size_t i = two_usec_after_muon_peak; i < summed_wf.size(); ++i){
	if(summed_wf[i] < michel_peak_value){ //NEED TO CHECK IF WF IS INVERTED
	  michel_peak_value = summed_wf[i];
	  michel_peak_value_position = i;
	}
    }

    if(michel_peak_value > muon_peak_value*0.15){
        // so there's a michel electron in this spectrum! yay!
        // now we need to find where it starts and ends by looping over wf
        
        // let's start by looping over the wf before michel_peak_value_position and define the
        // start of the peak as where the wf returns to the baseline
        for(size_t michel_it = michel_peak_value_position; michel_it > michel_peak_value_position - 150; --michel_it){
            if(summed_wf[michel_it] < baseline_threshold and summed_wf[michel_it] > (-1*baseline_threshold)){
              // so we're within a few counts of the baseline
              michel_start_index = michel_it;
              break;
            }
          }
          // now let's loop over the wf after michel_peak_value_position and define the
          // end of the peak as where the wf drops off less steeply
          for(size_t michel_it = michel_peak_value_position; michel_it < std::min(michel_peak_value_position + 150, summed_wf.size()-1); ++michel_it){
              double diff_between_bins = summed_wf[michel_it] - summed_wf[michel_it+1];
              if(diff_between_bins < diff_between_bins_threshold and diff_between_bins > (-1*diff_between_bins_threshold)){
                // now we should be within the triplet light region
                // this might need some fine tuning
                michel_end_index = michel_it;
                break;
              }
            }

        }

}


void CCMMichelElectrons::DAQ(I3FramePtr frame) {
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
    std::cout << "wf size = " << waveforms->size() << std::endl;

    // a shared pointer to store the total charge per pmt for each event
    boost::shared_ptr<I3Map<CCMPMTKey, CCMWaveformUInt16>> total_charge_per_pmt = boost::make_shared<I3Map<CCMPMTKey, CCMWaveformUInt16>>();
    //boost::shared_ptr<I3Vector<double>> summed_wf_per_frame(new I3Vector<double>(8000)); //hard coding in 8000 samples...need to fix at some point
    boost::shared_ptr<I3Vector<size_t>> michel_start_index_per_frame(new I3Vector<size_t>(1));
    boost::shared_ptr<I3Vector<size_t>> michel_end_index_per_frame(new I3Vector<size_t>(1));
    boost::shared_ptr<I3Vector<int>> summed_wf_reference_time(new I3Vector<int>(1));

    // let's add up all the charge in the detector
    int fixed_position;
    std::deque<double> summed_wf = SumWaveforms(frame, waveforms, baseline_mode, fixed_position);
    boost::shared_ptr<I3Vector<double>> summed_wf_per_frame = boost::make_shared<I3Vector<double>>(summed_wf.begin(), summed_wf.end());
    
    for(size_t wf_it = 0; wf_it < summed_wf.size(); ++wf_it){
       summed_wf_per_frame->operator[](wf_it) = summed_wf[wf_it];
    }
    
    summed_wf_reference_time->operator[](0) = fixed_position;
    
    // now let's see if we have a michele candidate 2 microseconds from the trigger
    size_t michel_start_index = 0;
    size_t michel_end_index = 0;
    FindMicheleRegion(summed_wf, michel_start_index, michel_end_index);
    michel_start_index_per_frame->operator[](0) = michel_start_index;
    michel_end_index_per_frame->operator[](0) = michel_end_index;

    if(michel_start_index != 0 and michel_end_index!= 0 ){
      // that means we've identified a michel > 2 microseconds after the muon peak!
      // let's reconstruct a vertex!

      double event_radius = 0;
      double event_azimuth = 0;
      double event_z = 0;

//      std::tie(event_radius, event_azimuth, event_z) = VertexReconstruction(waveforms, baseline_mode, summed_wf, michel_start_index, michel_end_index);
//                                                      // pretty sure this is the wrong way to call this function....
//                                                      // something to look into
//
//      if(event_radius < 10.0 and event_z < 5 and event_z > -5){
//        // found a michel in the middle-ish of the detector! now let's save charge per PMT
//
//        for(std::pair<CCMPMTKey const, double> const & it : baseline_mode){
//           CCMPMTKey key = it.first;
//           double mode = it.second; //baseline mode is negative!!!
//           uint32_t channel = pmt_channel_map_[key];
//           CCMWaveformUInt16 const & waveform = waveforms->at(channel);
//
//           double total_charge = 0;
//           ProcessWaveform(waveform, mode, michel_start_index, michel_end_index, total_charge);
//           total_charge_per_pmt->emplace(std::make_pair(key, total_charge));
//        }
//      }
    }
    frame->Put("SummedWaveform", summed_wf_per_frame);
    frame->Put("MichelStartIndex", michel_start_index_per_frame);
    frame->Put("MichelEndIndex", michel_end_index_per_frame);
    frame->Put("SummedWaveformFixedPosition", summed_wf_reference_time);
    std::cout << "finished finding Michele electrons!" << std::endl;
    PushFrame(frame);
    }
}

void CCMMichelElectrons::ProcessWaveform(CCMWaveformUInt16 const & waveform, double const & mode, size_t const & michel_start_index, size_t const & michel_end_index, double & total_charge){
    // get the vector of samples from the CCMWaveform object;
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();

    if (samples.size() == 0) {
    	return;
    }
    // let's subtract the baseline off wf and add up all of the wf values in our michel range
    for(size_t michel_it = michel_start_index; michel_it <= michel_end_index; ++michel_it){
       total_charge += (samples[michel_it] + mode);
    }

}

void CCMMichelElectrons::Finish() {
}


