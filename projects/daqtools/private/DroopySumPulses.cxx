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

class DroopySumPulses : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;
    std::string pulses_name_;
    std::string droop_pulses_name_;
    bool already_summed_;
    std::string counts_name_;
    std::string peak_positions_name_;

    // Internal state
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    std::map<CCMPMTKey, WaveformAccumulator> summed_waveforms_;
    std::map<CCMPMTKey, WaveformAccumulator> droop_summed_waveforms_;

public:
    DroopySumPulses(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(DroopySumPulses);

void DroopySumPulses::ProcessFrame(I3FramePtr frame) {
    I3Map<CCMPMTKey, std::vector<double>> const & pulse_samples = frame->Get<I3Map<CCMPMTKey, std::vector<double>> const>(pulses_name_);
    I3Map<CCMPMTKey, std::vector<double>> const & pulse_samples_droop_corrected = frame->Get<I3Map<CCMPMTKey, std::vector<double>> const>(droop_pulses_name_);
    if(already_summed_) {
        I3Map<CCMPMTKey, int> const & peak_positions = frame->Get<I3Map<CCMPMTKey, int> const>(peak_positions_name_);
        I3Map<CCMPMTKey, std::vector<unsigned int>> const & counts = frame->Get<I3Map<CCMPMTKey, std::vector<unsigned int>> const>(counts_name_);
        for(std::pair<CCMPMTKey const, std::vector<double>> const & p : pulse_samples) {
            CCMPMTKey pmt_key = p.first;
            std::vector<double> const & wf = p.second;
            int pos = peak_positions.at(pmt_key);
            std::vector<unsigned int> count = counts.at(pmt_key);
            summed_waveforms_[pmt_key].AddWaveform(wf, pos, count);
        }
        for(std::pair<CCMPMTKey const, std::vector<double>> const & p : pulse_samples_droop_corrected) {
            CCMPMTKey pmt_key = p.first;
            std::vector<double> const & wf_droop_corrected = p.second;
            int pos = peak_positions.at(pmt_key);
            std::vector<unsigned int> count = counts.at(pmt_key);
            droop_summed_waveforms_[pmt_key].AddWaveform(wf_droop_corrected, pos, count);
        }
    } else {
        for(std::pair<CCMPMTKey const, std::vector<double>> const & p : pulse_samples) {
            CCMPMTKey pmt_key = p.first;
            std::vector<double> const & wf = p.second;
            int pos = std::distance(wf.begin(), std::max_element(wf.begin(), wf.end()));
            summed_waveforms_[pmt_key].AddWaveform(wf, pos);
        }
        for(std::pair<CCMPMTKey const, std::vector<double>> const & p : pulse_samples_droop_corrected) {
            CCMPMTKey pmt_key = p.first;
            std::vector<double> const & wf_droop_corrected = p.second;
            int pos = std::distance(wf_droop_corrected.begin(), std::max_element(wf_droop_corrected.begin(), wf_droop_corrected.end()));
            droop_summed_waveforms_[pmt_key].AddWaveform(wf_droop_corrected, pos);
        }
    }
}

DroopySumPulses::DroopySumPulses(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_("") {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("PulseSamplesName", "Key for pulse samples", std::string(""));
    AddParameter("PulseSamplesNameDroopCorrected", "Key for pulse samples", std::string(""));
    AddParameter("AlreadySummed", "Is the input pulses that have already been added together? If so, use the counts and positons stored in the frame", bool(false));
    AddParameter("CountsName", "Key for summed pulses counts", std::string(""));
    AddParameter("PeakPositionsName", "Key for summed pulses reference positions", std::string(""));
}

void DroopySumPulses::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("PulseSamplesName", pulses_name_);
    GetParameter("PulseSamplesNameDroopCorrected", droop_pulses_name_);
    GetParameter("AlreadySummed", already_summed_);
    if(pulses_name_ == "") {
        if(already_summed_) {
            pulses_name_ = "SummedPulses";
        } else {
            pulses_name_ = "PulseSamples";
        }
    }
    if(droop_pulses_name_ == "") {
        if(already_summed_) {
            droop_pulses_name_ = "SummedPulses";
        } else {
            droop_pulses_name_ = "PulseSamplesDroopCorrected";
        }
    }
    GetParameter("CountsName", counts_name_);
    if(counts_name_ == "" and already_summed_) {
        log_warn("CountsName not specified, so assuming default.");
        counts_name_ = pulses_name_ + "Counts";
    }
    GetParameter("PeakPositionsName", peak_positions_name_);
    if(peak_positions_name_ == "" and already_summed_) {
        log_warn("PeakPositionsName not specified, so assuming default.");
        peak_positions_name_ = pulses_name_ + "PeakPositions";
    }
}

void DroopySumPulses::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    PushFrame(frame);
}

void DroopySumPulses::DAQ(I3FramePtr frame) {
    ProcessFrame(frame);
}

void DroopySumPulses::Finish() {
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> summed_waveforms = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> summed_waveforms_droop_corrected = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
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
    for(std::pair<CCMPMTKey const, WaveformAccumulator> const & p : droop_summed_waveforms_) {
        CCMPMTKey pmt_key = p.first;
        std::deque<double> deque_wf = p.second.GetSummedWaveform();
        std::deque<unsigned int> deque_counts = p.second.GetCounts();
        std::vector<double> wf(deque_wf.begin(), deque_wf.end());
        std::vector<unsigned int> counts(deque_counts.begin(), deque_counts.end());
        summed_waveforms_droop_corrected->insert({pmt_key, wf});
    }
    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::DAQ);
    frame->Put("SummedPulses", summed_waveforms);
    frame->Put("SummedPulsesDroopCorrected", summed_waveforms_droop_corrected);
    frame->Put("SummedPulsesPeakPositions", waveform_peak_positions);
    frame->Put("SummedPulsesCounts", waveform_counts);
    PushFrame(frame);
    Flush();
}
