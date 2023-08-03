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
#include <algorithm>
#include <iterator>

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

class SumWaveforms: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    std::set<CCMTriggerKey> allowed_trigger_keys;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    I3Map<CCMPMTKey, CCMOMGeo> pmt_geo_;
    void Geometry(I3FramePtr frame);
    bool IsCorrectTrigger(I3FramePtr frame);
    std::tuple<std::deque<double>, std::deque<unsigned int>> ComputeSummedWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position);
public:
    SumWaveforms(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(SumWaveforms);

SumWaveforms::SumWaveforms(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
    AddParameter("AllowedTriggerKeys", "Trigger keys to sum", I3Vector<CCMTriggerKey>());
}


void SumWaveforms::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    I3Vector<CCMTriggerKey> keys;
    GetParameter("AllowedTriggerKeys", keys);
    std::copy(keys.begin(), keys.end(), std::inserter(allowed_trigger_keys, allowed_trigger_keys.end()));
}


void SumWaveforms::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    pmt_geo_ = geo.pmt_geo;
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
    PushFrame(frame);
}

bool SumWaveforms::IsCorrectTrigger(I3FramePtr frame) {
    if(allowed_trigger_keys.size() == 0)
        return true;
    if(not frame->Has(nim_pulses_name_)) {
        log_warn(("No key named " + nim_pulses_name_ + " present in frame").c_str());
        return false;
    }
    boost::shared_ptr<NIMLogicPulseSeriesMap const> nim_pulses = frame->Get<boost::shared_ptr<NIMLogicPulseSeriesMap const>>(nim_pulses_name_);
    if(not nim_pulses) {
        log_warn(("No NIMLogicPulseSeriesMap named " + nim_pulses_name_ + " present in frame").c_str());
        return false;
    }
    bool found_allowed_trigger = false;
    for(auto trigger_key : allowed_trigger_keys) {
        NIMLogicPulseSeriesMap::const_iterator it = nim_pulses->find(trigger_key);
        if(it != nim_pulses->end()) {
            if(it->second.size() >= 0) {
                found_allowed_trigger = true;
                break;
            }
        }
    }
    // We found a NIM pulse on the cosmic trigger,
    // therefore this is a cosmic frame
    return found_allowed_trigger;
}


std::tuple<std::deque<double>, std::deque<unsigned int>> SumWaveforms::ComputeSummedWaveforms(I3FramePtr frame, boost::shared_ptr<const CCMWaveformUInt16Series> const & waveforms, I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode, int & fixed_position){

    //reset summed_waveforms_ !!!
    WaveformAccumulator summed_waveforms_;

    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        CCMOMGeo::OMType const & omtype = pmt_geo_[key].omtype;
        if(omtype != CCMOMGeo::OMType::CCM8inUncoated and omtype != CCMOMGeo::OMType::CCM8inCoated)
            continue;
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
    std::deque<unsigned int> summed_wf_counts = summed_waveforms_.GetCounts();
    fixed_position = summed_waveforms_.GetFixedPosition();
    return {summed_wf, summed_wf_counts};
}


void SumWaveforms::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }

    if(not IsCorrectTrigger(frame)) {
        PushFrame(frame);
        return;
    }

    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    // a shared pointer to store the total charge per pmt for each event
    boost::shared_ptr<I3Vector<int>> summed_wf_reference_time(new I3Vector<int>(1));

    // let's add up all the charge in the detector
    int fixed_position;
    std::tuple<std::deque<double>, std::deque<unsigned int>> summed_wf = ComputeSummedWaveforms(frame, waveforms, baseline_mode, fixed_position);
    boost::shared_ptr<I3Vector<double>> summed_wf_per_frame = boost::make_shared<I3Vector<double>>(std::get<0>(summed_wf).begin(), std::get<0>(summed_wf).end());
    boost::shared_ptr<I3Vector<unsigned int>> summed_wf_counts_per_frame = boost::make_shared<I3Vector<unsigned int>>(std::get<1>(summed_wf).begin(), std::get<1>(summed_wf).end());

    summed_wf_reference_time->operator[](0) = fixed_position;

    frame->Put("SummedWaveform", summed_wf_per_frame);
    frame->Put("SummedWaveformCounts", summed_wf_counts_per_frame);
    frame->Put("SummedWaveformFixedPosition", summed_wf_reference_time);
    PushFrame(frame);
}

void SumWaveforms::Finish() {
}

