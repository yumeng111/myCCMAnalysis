#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <string>
#include <random>
#include <iostream>
#include <algorithm>

#include <icetray/I3Frame.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

#include "daqtools/OnlineRobustStats.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"

class SumPulses : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;

    // Internal state
    bool geo_seen;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    std::map<CCMPMTKey, WaveformAccumulator> summed_waveforms_;

    bool seen_frame;
    I3FramePtr last_frame;

public:
    SumPulses(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(SumPulses);

void SumPulses::ProcessFrame(I3FramePtr frame) {
    I3Map<CCMPMTKey, std::vector<double>> const & pulse_samples = frame->Get<I3Map<CCMPMTKey, std::vector<double>> const>("PulseSamples");
    for(std::pair<CCMPMTKey const, std::vector<double>> const & p : pulse_samples) {
        CCMPMTKey pmt_key = p.first;
        std::vector<double> const & wf = p.second;
        int pos = std::distance(wf.begin(), std::max_element(wf.begin(), wf.end()));
        summed_waveforms_[pmt_key].AddWaveform(wf, pos);
    }
}

SumPulses::SumPulses(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_(""), geo_seen(false), seen_frame(false), last_frame(nullptr) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
}

void SumPulses::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
}

void SumPulses::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);


    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;

    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey pmt_key = p.first;
        summed_waveforms_.insert({pmt_key, WaveformAccumulator()});
    }

    geo_seen = true;
    PushFrame(frame);
}

void SumPulses::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    ProcessFrame(frame);
    if(seen_frame) {
        PushFrame(last_frame);
    }
    last_frame = frame;
    seen_frame = true;
}

void SumPulses::Finish() {
    if(seen_frame) {
        boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> summed_waveforms = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
        boost::shared_ptr<I3Map<CCMPMTKey, int>> waveform_peak_positions = boost::make_shared<I3Map<CCMPMTKey, int>>();
        boost::shared_ptr<I3Map<CCMPMTKey, std::vector<unsigned int>>> waveform_counts = boost::make_shared<I3Map<CCMPMTKey, std::vector<unsigned int>>>();
        for(std::pair<CCMPMTKey const, WaveformAccumulator> const & p : summed_waveforms_) {
            CCMPMTKey pmt_key = p.first;
            std::deque<double> deque_wf = p.second.GetSummedWaveform();
            std::deque<unsigned int> deque_counts = p.second.GetCounts();
            std::vector<double> wf(deque_wf.begin(), deque_wf.end());
            std::vector<unsigned int> counts(deque_counts.begin(), deque_counts.end());
            summed_waveforms->insert({pmt_key, wf});
            waveform_peak_positions->insert({pmt_key, p.second.GetFixedPosition()});
            waveform_counts->insert({pmt_key, counts});
        }
        last_frame->Put("SummedPulses", summed_waveforms);
        last_frame->Put("SummedPulsesPeakPositions", waveform_peak_positions);
        last_frame->Put("SummedPulsesCounts", waveform_counts);
        PushFrame(last_frame);
    }
    seen_frame = false;
    last_frame = nullptr;
    Flush();
}
