// module to check that waveform accum is working right!
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
#include <random>

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

class TestWaveformAccum: public I3Module {
    WaveformAccumulator summed_waveforms_;
    std::string geometry_name_;
    bool geo_seen;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    std::deque<double> SumVectors(I3FramePtr frame, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode);
public:
    TestWaveformAccum(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(TestWaveformAccum);

TestWaveformAccum::TestWaveformAccum(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void TestWaveformAccum::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void TestWaveformAccum::Geometry(I3FramePtr frame) {
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


std::deque<double> TestWaveformAccum::SumVectors(I3FramePtr frame, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode){
    
    //reset summed_waveforms_ !!!
    WaveformAccumulator summed_waveforms_;

    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second; 
        double mode = value.baseline; //baseline mode is negative!!!
	uint32_t channel = pmt_channel_map_[key];
	
	// now let's make our vector of a random size that is full of all 0s except a 1 in a single position
	std::random_device rd;  // a seed source for the random number engine
    	std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
    	std::uniform_int_distribution<> vec_size(7500, 8500);
    	std::uniform_int_distribution<> fixed_pos(0, 7500);

	std::vector<double> random_vec(vec_size(gen), 0.0);
	double pos = fixed_pos(gen);
	random_vec[pos] = 1;
	
	summed_waveforms_.AddWaveform(random_vec, pos);
    }

    std::deque<double> summed_wf = summed_waveforms_.GetSummedWaveform();
    return summed_wf;

}

void TestWaveformAccum::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    // let's add up some randome vectors
    std::deque<double> summed_wf = SumVectors(frame, baseline_mode);
    boost::shared_ptr<I3Vector<double>> summed_wf_per_frame = boost::make_shared<I3Vector<double>>(summed_wf.begin(), summed_wf.end());
    
    for(size_t wf_it = 0; wf_it < summed_wf.size(); ++wf_it){
       summed_wf_per_frame->operator[](wf_it) = summed_wf[wf_it];
    }

    frame->Put("WaveformAccumCheck", summed_wf_per_frame);
    std::cout << "finished checking waveform accum!" << std::endl;
    PushFrame(frame);
   
}


void TestWaveformAccum::Finish() {
}

