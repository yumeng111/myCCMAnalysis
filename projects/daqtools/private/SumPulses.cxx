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
    std::string pulses_name_;
    bool already_summed_;
    std::string counts_name_;
    std::string peak_positions_name_;
    double min_counts_;
    double max_counts_;
    std::string summed_pulses_name_;
    std::string summed_pulses_peak_pos_name_;
    std::string summed_pulses_counts_name_;

    // Internal state
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    std::map<CCMPMTKey, WaveformAccumulator> summed_waveforms_;

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
    if(already_summed_) {
        I3Map<CCMPMTKey, std::vector<double>> const & pulse_samples = frame->Get<I3Map<CCMPMTKey, std::vector<double>> const>(pulses_name_);
        I3Map<CCMPMTKey, std::vector<double>> const & summed_pulse_samples = frame->Get<I3Map<CCMPMTKey, std::vector<double>> const>(pulses_name_);
        I3Map<CCMPMTKey, int> const & peak_positions = frame->Get<I3Map<CCMPMTKey, int> const>(peak_positions_name_);
        I3Map<CCMPMTKey, std::vector<unsigned int>> const & counts = frame->Get<I3Map<CCMPMTKey, std::vector<unsigned int>> const>(counts_name_);
        for(std::pair<CCMPMTKey const, std::vector<double>> const & p : summed_pulse_samples) {
            CCMPMTKey pmt_key = p.first;
            std::vector<double> const & wf = p.second;
            int pos = peak_positions.at(pmt_key);
            std::vector<unsigned int> count = counts.at(pmt_key);
            summed_waveforms_[pmt_key].AddWaveform(wf, pos, count);
        }
    } else {
        I3Map<CCMPMTKey, std::vector<std::vector<double>>> const & pulse_samples = frame->Get<I3Map<CCMPMTKey, std::vector<std::vector<double>>> const>(pulses_name_);
        for(std::pair<CCMPMTKey const, std::vector<std::vector<double>>> const & p : pulse_samples) {
            CCMPMTKey pmt_key = p.first;
            std::vector<std::vector<double>> const & vector_of_wfs = p.second;
            // now let's loop over our vector of pulse samples
            for (size_t wfs_it = 0; wfs_it < vector_of_wfs.size(); ++wfs_it){
                // let's find the max adc counts of each pulse
                std::vector<double> current_wf = vector_of_wfs[wfs_it];
                double max_adc_counts = 0;
                for (size_t wf_it = 0; wf_it < current_wf.size(); ++wf_it){
                    max_adc_counts = std::max(max_adc_counts, current_wf[wf_it]);
                }
                // now let's check if we're within the desired range
                if (max_adc_counts < max_counts_ and max_adc_counts > min_counts_){
                    // so we're within our desired range! yay
                    int pos = std::distance(current_wf.begin(), std::max_element(current_wf.begin(), current_wf.end()));
                    summed_waveforms_[pmt_key].AddWaveform(current_wf, pos);
                }
            }
        }
    }
}

SumPulses::SumPulses(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_("") {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("PulseSamplesName", "Key for pulse samples", std::string(""));
    AddParameter("AlreadySummed", "Is the input pulses that have already been added together? If so, use the counts and positons stored in the frame", bool(false));
    AddParameter("CountsName", "Key for summed pulses counts", std::string(""));
    AddParameter("PeakPositionsName", "Key for summed pulses reference positions", std::string(""));
    AddParameter("MinCountsThreshold", "Minimum ADC counts to accept for pulses", double(20));
    AddParameter("MaxCountsThreshold", "Maximum ADC counts to accept for pulses", double(50));
    AddParameter("SummedPulsesName", "Name for summed pulses", std::string("SummedPulses"));
    AddParameter("SummedPulsesPeakPosName", "Name for summed pulses peak pos", std::string("SummedPulsesPeakPositions"));
    AddParameter("SummedPulsesCountsName", "Name for summed pulses counts", std::string("SummedPulsesCounts"));
}

void SumPulses::Configure() {
    GetParameter("SummedPulsesName", summed_pulses_name_);
    GetParameter("SummedPulsesPeakPosName", summed_pulses_peak_pos_name_);
    GetParameter("SummedPulsesCountsName", summed_pulses_counts_name_);
    GetParameter("MinCountsThreshold", min_counts_);
    GetParameter("MaxCountsThreshold", max_counts_);
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("PulseSamplesName", pulses_name_);
    GetParameter("AlreadySummed", already_summed_);
    if(pulses_name_ == "") {
        if(already_summed_) {
            pulses_name_ = "SummedPulses";
        } else {
            pulses_name_ = "PulseSamples";
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

void SumPulses::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    PushFrame(frame);
}

void SumPulses::DAQ(I3FramePtr frame) {
    ProcessFrame(frame);
}

void SumPulses::Finish() {
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
    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::DAQ);
    frame->Put(summed_pulses_name_, summed_waveforms);
    frame->Put(summed_pulses_peak_pos_name_, waveform_peak_positions);
    frame->Put(summed_pulses_counts_name_, waveform_counts);
    PushFrame(frame);
    Flush();
}
