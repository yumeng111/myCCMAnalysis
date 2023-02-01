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

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Orientation.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

class CCMTriggerMerger : public I3Module {
    std::string geometry_name_;
    std::string digital_readout_name_;
    std::string triggers_name_;
    std::string daq_config_name_;
    std::string waveforms_output_;

    boost::shared_ptr<const CCMAnalysis::Binary::CCMDAQConfig> last_daq_config_;
    boost::shared_ptr<const CCMGeometry> last_geometry_;
    std::map<size_t, size_t> board_idx_to_machine_idx;
    std::map<size_t, std::pair<size_t, size_t>> board_idx_to_channel_idxs;
    std::deque<I3FramePtr> daq_frame_cache_;

    size_t triggers_seen;
    size_t merged_triggers_output;
    size_t total_triggers_output;

    bool TriggersOverlap(boost::shared_ptr<const std::vector<CCMTrigger>> trigger0, boost_shared_ptr<const std::vector<CCMTrigger> trigger1, size_t idx, uint32_t extra_samples = 0);
    bool TriggersOverlap(boost::shared_ptr<const std::vector<CCMTrigger>> trigger0, boost_shared_ptr<const std::vector<CCMTrigger> trigger1);
    void CacheGeometryFrame(I3FramePtr frame);
    void CacheDAQFrame(I3FramePtr frame);
    void PushMergedTriggerFrames(bool last_frame = false);
    I3FramePtr MergeTriggerFrames(std::vector<I3FramePtr> & frames);
public:
    CCMTriggerMerger(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMTriggerMerger);

namespace detail {
    std::string tolower(std::string const & s) {
        std::string r(s);
        std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c){return std::tolower(c);});
        return r;
    }

    bool ChannelEmpty(boost::shared_ptr<const std::vector<CCMTrigger>> trigger, size_t channel_idx) {
        return (*trigger)[0].channel_sizes[channel_idx] == 0 or (*trigger)[0].channel_masks[channel_idx] == 0;
    }

    bool ChannelsEmpty(boost::shared_ptr<const std::vector<CCMTrigger>> trigger, size_t board_idx) {
        std::pair<size_t, size_t> channel_idxs = board_idx_to_channel_idxs[board_idx];
        bool empty = true;
        for(size_t channel_idx=channesl_idxs.first; channel_idx<channel_idxs.second; ++channel_idx) {
            if(not ChannelEmpty(trigger, channel_idx)) {
                empty = false;
                break;
            }
        }
        return empty;
    }

}

CCMTriggerMerger::CCMTriggerMerger(const I3Context& context) : I3Module(context) {
    triggers_seen = 0;
    merged_triggers_output = 0;
    total_triggers_output = 0;
    // CCMWaveforms output parameter
    // Criteria for merging
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDigitalReadoutName", "Key for vector of vector of samples that represents the digitizer readout", std::string("CCMDigitalReadout"));
    AddParameter("CCMTriggersName", "Key for vector of CCMTriggers", std::string("CCMTriggers"));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMDAQConfig>::value()));
    AddParameter("CCMWaveformsOutput", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
}

void CCMTriggerMerger::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMDigitalReadout", digital_readout_name_);
    GetParameter("CCMTriggersName", triggers_name_);
    GetParameter("CCMDAQConfigName", daq_config_name_);
    GetParameter("CCMWaveformsOutput", waveforms_output_);
}

bool CCMTriggerMerger::TriggersOverlap(boost::shared_ptr<const std::vector<CCMTrigger>> trigger0, boost_shared_ptr<const std::vector<CCMTrigger> trigger1, size_t idx, uint32_t extra_samples) {
    constexpr long double ns_per_cycle = 8;
    constexpr long double ns_per_sample = 2;
    constexpr long double samples_per_cycle = ns_per_cyle / ns_per_sample;

    int32_t cycle_start_diff = CCMAnalysis::Binary::subtract_times((*trigger1)[0].board_times[idx], (*trigger0)[0].board_times[idx]);
    long double samples_start_diff = cycle_diff * samples_per_cycle;
    long double expected_samples = board_idx_to_machine_idx[idx] + extra_samples;
    return samples_start_diff < expected_samples;
}

bool CCMTriggerMerger::TriggersOverlap(boost::shared_ptr<const std::vector<CCMTrigger>> trigger0, boost_shared_ptr<const std::vector<CCMTrigger> trigger1) {
    bool merge_needed = false;
    bool merge_possible = true;
    for(size_t board_idx=0; board_idx<board_idx_to_machine_idx.size(); ++board_idx) {
        bool trigger0_empty = ChannelsEmpty(trigger0, board_idx);
        bool trigger1_empty = ChannelsEmpty(trigger1, board_idx);
        bool triggers_overlap = false;
        if((not trigger0_empty) and (not trigger1_empty))
            triggers_overlap = TriggersOverlap(trigger0, trigger1, board_idx);
        merge_needed |= triggers_overlap;
        // Exclude the case we have two full non-overlapping triggers
        merge_possible &= not (triggers_overlap or trigger0_empty or trigger1_empty);
    }

    // Merging two full length triggers that do not overlap should not be possible...
    if(merge_needed and not merge_possible) {
        // Check that triggers that would exclude a merge are at least close to overlapping
        constexpr uint32_t num_extra_samples_allowed = 12;
        bool merge_allowed = true;
        for(size_t board_idx=0; board_idx<board_idx_to_machine_idx.size(); ++board_idx) {
            bool trigger0_empty = ChannelsEmpty(trigger0, board_idx);
            bool trigger1_empty = ChannelsEmpty(trigger1, board_idx);
            bool triggers_overlap = false;
            if((not trigger0_empty) and (not trigger1_empty))
                triggers_overlap = TriggersOverlap(trigger0, trigger1, board_idx, num_extra_samples_allowed);
            // Exclude the case we have two full non-overlapping triggers
            merge_allowed &= not (triggers_overlap or trigger0_empty or trigger1_empty);
        }
        if(not merge_allowed)
            log_fatal("Need to merge two triggers, but merge is not possible because triggers from some boards overlap and others do not.");
    }
    return merge_needed;
}

void CCMTriggerMerger::CacheGeometryFrame(I3FramePtr frame) {
    if(not frame->Has(daq_config_name_))
        log_fatal("Frame does not contain daq condig" + daq_config_name_);
    if(not frame->Has(geometry_name_))
        log_fatal("Frame does not contain geometry: " + geometry_name_);
    last_daq_config_ = frame->Get<boost::shared_ptr<const CCMAnalysis::Binary::CCMDAQConfig>>(daq_config_name_);
    last_geometry_ = frame->Get<boost::shared_ptr<const CCMGeometry>>(geometry_name_);

    std::vector<size_t> n_boards;
    size_t total_boards = 0;
    size_t board_idx = 0;
    size_t total_channels = 0;
    size_t channel_idx = 0;
    board_idx_to_machine_idx.clear();
    board_idx_to_channel_idxs.clear();
    for(size_t machine_idx=0; machine_idx<last_daq_config_->machine_configurations.size(); ++machine_idx) {
        size_t n_boards = last_daq_config_->machine_configurations[i].num_digitizer_boards;
        for(total_boards += n_boards; board_idx<total_boards; ++board_idx) {
            size_t num_channels = last_daq_config_->digitizer_boards[board_idx].channels.size();
            board_idx_to_machine_idx[board_idx] = machine_idx;
            board_idx_to_channel_idxs[board_idx] = {channel_idx, channel_idx + num_channels};
            channel_idx += num_channels;
        }
    }
    last_daq_config_->machine_configurations;
}

void CCMTriggerMerger::CacheDAQFrame(I3FramePtr frame) {
    if(not frame->Has(digital_readout_name_))
        log_fatal("Frame does not contain digital readout: " + digital_readout_name_);
    if(not frame->Has(triggers_name_))
        log_fatal("Frame does not contain triggers: " + triggers_name_);
    daq_frame_cache_.push_back(frame);
}

void CCMTriggerMerger::PushMergedTriggerFrames(bool last_frame) {
    if(daq_frame_cache_.size() == 0 or (daq_frame_cache_.size() < 2 and not last_frame))
        return;

    std::vector<std::vector<I3FramePtr>> grouped_frames({{daq_frame_cache_[0]}});
    I3FramePtr last_frame = daq_frame_cache_[0];

    for(size_t frame_idx=1; frame_idx<daq_frame_cache_.size(); ++frame_idx) {
        I3FramePtr current_frame = daq_frame_cache_[frame_idx];
        bool triggers_overlap = TriggersOverlap(
                last_frame->Get<boost::shared_ptr<const I3Vector<CCMTrigger>>>(triggers_name_),
                current_frame->Get<boost::shared_ptr<const I3Vector<CCMTrigger>>>(triggers_name_));
        if(triggers_overlap) {
            grouped_frames.back().push_back(current_frame)
        } else {
            grouped_frames.emplace_back();
            grouped_frames.back().push_back(current_frame);
        }
        last_frame = current_frame;
    }
    for(size_t group_idx=0; group_idx + (1 - last_frame) < grouped_frames.size(); ++group_idx) {
        I3FramePtr frame = MergeTriggerFrames(grouped_frames[group_idx]);
        if(grouped_frames[group_idx].size() > 1)
            merged_triggers_output += 1;
        total_triggers_output += 1;
        PushFrame(frame);
    }
    size_t num_to_pop = daq_frame_cache_.size() - (last_frame ? 0 : grouped_frames.back().size());
    for(size_t i=0; i<num_to_pop; ++i)
        daq_frame_cache_.pop_front();
}

I3FramePtr CCMTriggerMerger::MergeTriggerFrames(std::vector<I3FramePtr> & frames) {
    if (frames.size() == 0) {
        return new I3Frame(I3Frame::DAQ);
    } else if(frames.size() == 1) {
        return frames[0];
    }

    I3FramePtr first_frame = frames[0];
    boost::shared_ptr<const std::vector<std::vector<uint16_t>>> first_readout = first_frame->Get<boost::shared_ptr<const std::vector<std::vector<uint16_t>>>>(digital_readout_name_);
    std::vector<CCMWaveform> output_waveforms;
    output_waveforms.reserve(first_readout->size());
    std::vector<CCTrigger> output_triggers;
    output_triggers.reserve(frames.size());

    std::vector<uint32_t> ref_times;
    std::vector<double> start_times;
    for(size_t f_idx=0; f_idx<frames.size(); ++i) {
        size_t readout_size = input_waveform->size();
        for(size_t i=0; i<readout_size; ++i) {
            if(not ChannelsEmpty()
        }
    }

    for(size_t w_idx=0; w_idx<first_readout->size(); ++w_idx) {
        output_waveforms.emplace_back();
        CCMWaveform & w = output_waveforms.back();
        size_t total_size = 0;
        for(size_t i=0; i<frames.size(); ++i) {
            total_size += frames[i]->Get<boost::shared_ptr<const std::vector<std::vector<uint16_t>>>>(digitial_readout_name_)->size();
        }
        std::vector<uint16_t> waveform(total_size);
        size_t last_pos = 0;
        bool recorded_start_time = false;
        uint32_t start_time = 0;
        for(size_t i=0; i<frames.size(); ++i) {
            boost::shared_ptr<const std::vector<std::vector<uint16_t>>> input_waveform = frames[i]->Get<boost::shared_ptr<const std::vector<std::vector<uint16_t>>>>(digitial_readout_name_);
            boost::shared_ptr<const std::vector<CCMAnalysis::Binary::CCMTrigger>> input_trigger = frames[i]->Get<boost::shared_ptr<const std::vector<CCMAnalysis::Binary::CCMTrigger>>>(triggers_name_);
            if(input_waveform.size() > 0) {
                if(not recorded_start_time) {
                    start_time
                }
                std::copy(input_waveform.begin(), input_waveform.end(), waveform.begin() + last_pos);
                last_pos += input_waveform.size();
            }
        }
        w.SetWaveform(waveform);
        w.SetBinWidth(2);
        w.SetStartTime(0);
        w.SetWaveformInformation();
        w.SetSource(CCMSource::V1730);
    }

    I3FramePtr ouput_frame = new I3Frame(I3Frame::DAQ);
    output_frame->Put(waveforms_output_, output_waveforms, I3Frame::DAQ);
    output_frame->Put(triggers_name_, output_triggers, I3Frame::DAQ);
}

void CCMTriggerMerger::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() == I3Frame::Geometry) {
        CacheGeometryFrame();
        PushMergedTriggerFrames(true);
        PushFrame(frame);
        return;
    }

    if(frame->GetStop() != I3Frame::DAQ) {
        PushFrame(frame);
        return;
    }

    triggers_seen += 1;
    CacheDAQFrame(frame);
    PushMergedTriggerFrames();
}

void CCMTriggerMerger::Finish() {
    PushMergedTriggerFrames(true);
    Flush();
    log_notice_stream(
            "Triggers seen: " << triggers_seen << "\n" <<
            "Merged triggers ouput: " << merged_triggers_output << "\n" << 
            "Total triggers output: " << total_triggers_output
    );
    // log_notice_stream(
    //     "Merged " << offsets.size() << " DAQ streams into " << counter << " separate triggers. Encountered "
    //     << counter-incomplete_counter << " complete triggers and " << incomplete_counter << " incomplete triggers");
}
