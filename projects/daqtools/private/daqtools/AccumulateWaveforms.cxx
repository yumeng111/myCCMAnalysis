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
#include <icetray/I3PODHolder.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

#include "daqtools/OnlineRobustStats.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"

class AccumulateWaveforms : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string input_prefix_;
    std::string output_prefix_;
    double min_counts_;
    double max_counts_;
    bool consume_frames_;

    // Internal state
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    WaveformAccumulator summed_waveform_;

public:
    AccumulateWaveforms(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(AccumulateWaveforms);

void AccumulateWaveforms::ProcessFrame(I3FramePtr frame) {
    std::vector<double> const & samples = frame->Get<I3Vector<double> const>(input_prefix_ + "SummedWaveform");
    std::vector<unsigned int> const & counts = frame->Get<I3Vector<unsigned int> const>(input_prefix_ + "SummedWaveformCounts");
    I3PODHolder<int> const & fixed_position = frame->Get<I3PODHolder<int> const>(input_prefix_ + "SummedWaveformFixedPosition");
    summed_waveform_.AddWaveform(samples, fixed_position.value, counts);
}

AccumulateWaveforms::AccumulateWaveforms(const I3Context& context) : I3Module(context),
    geometry_name_("") {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("InputPrefix", "Prefix for the output of SumWaveforms", std::string(""));
    AddParameter("OutputPrefix", "Prefix for the module output", std::string(""));
    AddParameter("MinCountsThreshold", "Minimum ADC counts to accept for pulses", double(20));
    AddParameter("MaxCountsThreshold", "Maximum ADC counts to accept for pulses", double(50));
    AddParameter("ConsumeFrames", "Consume frames used as input?", bool(true));
}

void AccumulateWaveforms::Configure() {
    GetParameter("MinCountsThreshold", min_counts_);
    GetParameter("MaxCountsThreshold", max_counts_);
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("InputPrefix", input_prefix_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("ConsumeFrames", consume_frames_);
}

void AccumulateWaveforms::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    PushFrame(frame);
}

void AccumulateWaveforms::DAQ(I3FramePtr frame) {
    ProcessFrame(frame);
    if(not consume_frames_)
        PushFrame(frame);
}

void AccumulateWaveforms::Finish() {

    boost::shared_ptr<I3PODHolder<int32_t>> p_fixed_position = boost::make_shared<I3PODHolder<int32_t>>();

    std::deque<double> deque_samples = summed_waveform_.GetSummedWaveform();
    std::deque<unsigned int> deque_counts = summed_waveform_.GetCounts();
    int fixed_position = summed_waveform_.GetFixedPosition();

    boost::shared_ptr<I3Vector<double>> p_samples = boost::make_shared<I3Vector<double>>(deque_samples.begin(), deque_samples.end());
    boost::shared_ptr<I3Vector<unsigned int>> p_counts = boost::make_shared<I3Vector<unsigned int>>(deque_counts.begin(), deque_counts.end());

    p_fixed_position->value = fixed_position;

    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::DAQ);

    frame->Put((output_prefix_ + "AccumulatedWaveform").c_str(), p_samples);
    frame->Put((output_prefix_ + "AccumulatedWaveformCounts").c_str(), p_counts);
    frame->Put((output_prefix_ + "AccumulatedWaveformFixedPosition").c_str(), p_fixed_position);

    PushFrame(frame);
    Flush();
}
