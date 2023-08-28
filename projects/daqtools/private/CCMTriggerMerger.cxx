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
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"


// Detects adjacent DAQ frames that overlap in time
// Merges overlapping DAQ frames into a single frame by concatenating waveforms and triggers
// Converts vector of raw digitizer readouts into vector of CCMWaveform objects

class CCMTriggerMerger : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string digital_readout_name_;
    std::string triggers_name_;
    std::string daq_config_name_;
    std::string waveforms_output_;
    std::string board_time_offsets_name_;
    std::string trigger_times_name_;
    std::string first_trigger_time_output_;
    uint32_t num_extra_samples_allowed_;

    // Important constants
    constexpr static long double ns_per_cycle = 8;
    constexpr static long double ns_per_sample = 2;
    constexpr static long double samples_per_cycle = ns_per_cycle / ns_per_sample;

    // Information to cache for each geometry frame
    boost::shared_ptr<const CCMAnalysis::Binary::CCMDAQConfig> last_daq_config_;
    boost::shared_ptr<const CCMGeometry> last_geometry_;
    boost::shared_ptr<const I3Vector<I3Vector<int64_t>>> last_board_time_offsets_;
    size_t num_boards;
    std::vector<size_t> num_boards_by_machine;
    std::vector<size_t> num_channels_by_board;
    std::vector<size_t> num_samples_by_board;
    std::map<size_t, size_t> board_idx_to_machine_idx;
    std::map<size_t, std::pair<size_t, size_t>> board_idx_to_channel_idxs;
    std::map<size_t, size_t> channel_idx_to_board_idx;
    std::vector<int64_t> offsets_;

    // DAQ frames yet to be merged
    std::deque<I3FramePtr> daq_frame_cache_;

    // Counters
    size_t triggers_seen;
    size_t merged_triggers_output;
    size_t total_triggers_output;

    // Convenience functions
    bool TriggersOverlap(boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times0, boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times1, size_t idx, uint32_t extra_samples = 0);
    bool TriggersOverlap(boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times0, boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times1);
    void CacheGeometryFrame(I3FramePtr frame);
    void CacheDAQFrame(I3FramePtr frame);
    void PushMergedTriggerFrames(bool last_frame = false);
    bool ChannelsEmpty(boost::shared_ptr<const std::vector<CCMAnalysis::Binary::CCMTrigger>> trigger, size_t board_idx);
    bool ChannelsEmpty(boost::shared_ptr<const std::vector<std::pair<bool, int64_t>>> trigger, size_t board_idx);
    I3FramePtr MergeTriggerFrames(std::vector<I3FramePtr> & frames);
    std::vector<std::pair<bool, int64_t>> GetStartTimes(std::vector<I3FramePtr> const & frames);
    std::tuple<std::vector<std::pair<bool, int64_t>>, int64_t> ComputeRelativeStartTimes(std::vector<std::pair<bool, int64_t>> const & input_times);
    I3FramePtr TransformSingleTriggerFrame(I3FramePtr frame);
    I3FramePtr TransformMultipleTriggerFrames(std::vector<I3FramePtr> const & frames);
public:
    CCMTriggerMerger(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMTriggerMerger);

namespace detail {
    bool ChannelEmpty(boost::shared_ptr<const std::vector<CCMAnalysis::Binary::CCMTrigger>> trigger, size_t channel_idx) {
        return (*trigger)[0].channel_sizes[channel_idx] == 0 or (*trigger)[0].channel_masks[channel_idx] == 0;
    }
}

bool CCMTriggerMerger::ChannelsEmpty(boost::shared_ptr<const std::vector<CCMAnalysis::Binary::CCMTrigger>> trigger, size_t board_idx) {
    std::pair<size_t, size_t> channel_idxs = board_idx_to_channel_idxs[board_idx];
    bool empty = true;
    for(size_t channel_idx=channel_idxs.first; channel_idx<channel_idxs.second; ++channel_idx) {
        if(not detail::ChannelEmpty(trigger, channel_idx)) {
            empty = false;
            break;
        }
    }
    return empty;
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
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMAnalysis::Binary::CCMDAQConfig>::value()));
    AddParameter("BoardTimeOffsetsName", "Key for CCMBoardOffsets", std::string("BoardTimeOffsets"));
    AddParameter("CCMWaveformsOutput", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("TriggerTimesName", "Key for trigger times", std::string("TriggerTimes"));
    AddParameter("FirstTriggerTimeOutput", "Key for time of first trigger in frame relative to run start", std::string("FirstTriggerTime"));
    AddParameter("NumExtraSamplesForOverlap", "Number of extra samples allowed between subsequent board triggers when checking if a frame merge is allowed.", uint32_t(12));
}

void CCMTriggerMerger::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMDigitalReadoutName", digital_readout_name_);
    GetParameter("CCMTriggersName", triggers_name_);
    GetParameter("CCMDAQConfigName", daq_config_name_);
    GetParameter("BoardTimeOffsetsName", board_time_offsets_name_);
    GetParameter("CCMWaveformsOutput", waveforms_output_);
    GetParameter("TriggerTimesName", trigger_times_name_);
    GetParameter("FirstTriggerTimeOutput", first_trigger_time_output_);
    GetParameter("NumExtraSamplesForOverlap", num_extra_samples_allowed_);
}

bool CCMTriggerMerger::TriggersOverlap(boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times0, boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times1, size_t idx, uint32_t extra_samples) {
    int64_t cycle_start_diff = (*trigger_times1)[idx].second - (*trigger_times0)[idx].second;
    long double samples_start_diff = cycle_start_diff * samples_per_cycle;
    long double expected_samples = num_samples_by_board[idx] + extra_samples;
    return samples_start_diff < expected_samples;
}

bool CCMTriggerMerger::TriggersOverlap(boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times0, boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> trigger_times1) {
    bool merge_needed = false;
    bool merge_possible = true;
    for(size_t board_idx=0; board_idx<board_idx_to_machine_idx.size(); ++board_idx) {
        bool trigger0_empty = (*trigger_times0)[board_idx].first;
        bool trigger1_empty = (*trigger_times1)[board_idx].first;
        bool triggers_overlap = false;
        if((not trigger0_empty) and (not trigger1_empty))
            triggers_overlap = TriggersOverlap(trigger_times0, trigger_times1, board_idx);
        merge_needed |= triggers_overlap;
        // Exclude the case we have two full non-overlapping triggers
        merge_possible &= triggers_overlap or trigger0_empty or trigger1_empty;
    }

    // Merging two full length triggers that do not overlap should not be possible...
    if(merge_needed and not merge_possible) {
        // Check that triggers that would exclude a merge are at least close to overlapping
        bool merge_allowed = true;
        for(size_t board_idx=0; board_idx<board_idx_to_machine_idx.size(); ++board_idx) {
            bool trigger0_empty = (*trigger_times0)[board_idx].first;
            bool trigger1_empty = (*trigger_times1)[board_idx].first;
            bool triggers_overlap = false;
            if((not trigger0_empty) and (not trigger1_empty))
                triggers_overlap = TriggersOverlap(trigger_times0, trigger_times1, board_idx, num_extra_samples_allowed_);
            // Exclude the case we have two full non-overlapping triggers
            merge_allowed &= triggers_overlap or trigger0_empty or trigger1_empty;
        }
        if(not merge_allowed) {
            std::cerr << "Issue merging triggers. Triggers from some boards overlap, but others do not" << std::endl;
            std::cerr << "Trigger times:" << std::endl;
            for(size_t board_idx=0; board_idx<board_idx_to_machine_idx.size(); ++board_idx) {
                bool trigger0_empty = (*trigger_times0)[board_idx].first;
                bool trigger1_empty = (*trigger_times1)[board_idx].first;
                int64_t trigger0_time = (*trigger_times0)[board_idx].second;
                int64_t trigger1_time = (*trigger_times1)[board_idx].second;
                std::cerr << "Board " << board_idx << ": (" <<
                    (trigger0_empty ? "No trigger" : "Have trigger");
                if(not trigger0_empty)
                    std::cerr << " " << trigger0_time;
                std::cerr << ") (";
                std::cerr << (trigger1_empty ? "No trigger" : "Have trigger");
                if(not trigger1_empty)
                    std::cerr << " " << trigger1_time;
                std::cerr << ")";
                bool triggers_overlap = false;
                if((not trigger0_empty) and (not trigger1_empty)) {
                    triggers_overlap = TriggersOverlap(trigger_times0, trigger_times1, board_idx, num_extra_samples_allowed_);
                }
                std::cerr << " Overlap? " << (triggers_overlap ? "Yes" : "No");
                std::cerr << std::endl;
            }
            log_fatal("Need to merge two triggers, but merge is not possible because triggers from some boards overlap and others do not.");
        }
    }
    return merge_needed;
}

void CCMTriggerMerger::CacheGeometryFrame(I3FramePtr frame) {
    if(not frame->Has(daq_config_name_))
        log_fatal(("Frame does not contain daq config" + daq_config_name_).c_str());
    if(not frame->Has(geometry_name_))
        log_fatal(("Frame does not contain geometry: " + geometry_name_).c_str());
    last_daq_config_ = frame->Get<boost::shared_ptr<const CCMAnalysis::Binary::CCMDAQConfig>>(daq_config_name_);
    last_geometry_ = frame->Get<boost::shared_ptr<const CCMGeometry>>(geometry_name_);
    last_board_time_offsets_ = frame->Get<boost::shared_ptr<const I3Vector<I3Vector<int64_t>>>>(board_time_offsets_name_);

    offsets_.clear();
    for(size_t board_idx=0; board_idx<last_board_time_offsets_->size(); ++board_idx) {
        std::copy(
                (*last_board_time_offsets_)[board_idx].begin(),
                (*last_board_time_offsets_)[board_idx].end(),
                std::back_inserter(offsets_));
    }

    std::vector<size_t> n_boards;
    size_t total_boards = 0;
    size_t board_idx = 0;
    size_t total_channels = 0;
    size_t channel_idx = 0;
    board_idx_to_machine_idx.clear();
    board_idx_to_channel_idxs.clear();
    num_boards_by_machine.clear();
    num_channels_by_board.clear();
    num_samples_by_board.clear();
    channel_idx_to_board_idx.clear();
    for(size_t machine_idx=0; machine_idx<last_daq_config_->machine_configurations.size(); ++machine_idx) {
        size_t n_boards = last_daq_config_->machine_configurations[machine_idx].num_digitizer_boards;
        num_boards_by_machine.push_back(n_boards);
        for(total_boards += n_boards; board_idx<total_boards; ++board_idx) {
            size_t num_channels = last_daq_config_->digitizer_boards[board_idx].channels.size();
            num_channels_by_board.push_back(num_channels);
            num_samples_by_board.push_back(last_daq_config_->machine_configurations[machine_idx].num_samples);
            board_idx_to_machine_idx[board_idx] = machine_idx;
            board_idx_to_channel_idxs[board_idx] = {channel_idx, channel_idx + num_channels};
            channel_idx_to_board_idx[channel_idx] = board_idx;
            channel_idx += num_channels;
        }
    }
    last_daq_config_->machine_configurations;
}

void CCMTriggerMerger::CacheDAQFrame(I3FramePtr frame) {
    if(not frame->Has(digital_readout_name_))
        log_fatal(("Frame does not contain digital readout: " + digital_readout_name_).c_str());
    if(not frame->Has(triggers_name_))
        log_fatal(("Frame does not contain triggers: " + triggers_name_).c_str());
    daq_frame_cache_.push_back(frame);
}

void CCMTriggerMerger::PushMergedTriggerFrames(bool last_frame) {
    // The last_frame parameter indicates if there are no more frames in the stream to be added
    // Do nothing if we have no frames cached
    // Do nothing if last_frame == false and we have less than two frames cached
    if(daq_frame_cache_.size() == 0 or (daq_frame_cache_.size() < 2 and not last_frame))
        return;

    // Create a vector of frame groups
    // The last entry indicates the group currently being assembled
    std::vector<std::vector<I3FramePtr>> grouped_frames({{daq_frame_cache_[0]}});
    I3FramePtr prev_frame = daq_frame_cache_[0];

    // Iterate over remaining frames and add to the last group or create a new group accordingly
    for(size_t frame_idx=1; frame_idx<daq_frame_cache_.size(); ++frame_idx) {
        I3FramePtr current_frame = daq_frame_cache_[frame_idx];
        bool triggers_overlap = TriggersOverlap(
                prev_frame->Get<boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>>>(trigger_times_name_),
                current_frame->Get<boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>>>(trigger_times_name_));

        if(triggers_overlap) {
            // Overlapping triggers belong in the current group
            grouped_frames.back().push_back(current_frame);
        } else {
            // Non-overlapping triggers start a new group
            grouped_frames.emplace_back();
            grouped_frames.back().push_back(current_frame);
        }
        prev_frame = current_frame;
    }

    // Merge the groups and push the newly created frames
    // Only process the last group if there are no more input frames expected
    for(size_t group_idx=0; group_idx + (1 - last_frame) < grouped_frames.size(); ++group_idx) {
        I3FramePtr frame = MergeTriggerFrames(grouped_frames[group_idx]);
        if(grouped_frames[group_idx].size() > 1)
            merged_triggers_output += 1;
        total_triggers_output += 1;
        PushFrame(frame);
    }

    // Remove all the processed frames from the queue
    size_t num_to_pop = daq_frame_cache_.size() - (last_frame ? 0 : grouped_frames.back().size());
    for(size_t i=0; i<num_to_pop; ++i)
        daq_frame_cache_.pop_front();
}

std::vector<std::pair<bool, int64_t>> CCMTriggerMerger::GetStartTimes(std::vector<I3FramePtr> const & frames) {
    std::vector<std::pair<bool, int64_t>> start_times;
    for(size_t f_idx=0; f_idx<frames.size(); ++f_idx) {
        boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>> input_waveform = frames[f_idx]->Get<boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>>>(digital_readout_name_);
        boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> input_times = frames[f_idx]->Get<boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>>>(trigger_times_name_);
        size_t readout_size = input_waveform->size();
        // Initialize the start_times vector if it has not been already
        if(start_times.size() == 0)
            start_times = std::vector<std::pair<bool, int64_t>>(readout_size, std::pair<bool, int64_t>(false, 0));

        for(size_t i=0; i<readout_size; ++i) {
            // Enure start_times[i] contains the first valid start time in the sequence of frames
            if(not start_times[i].first)
               start_times[i] = (*input_times)[i];
        }
    }
    return start_times;
}

std::tuple<std::vector<std::pair<bool, int64_t>>, int64_t> CCMTriggerMerger::ComputeRelativeStartTimes(std::vector<std::pair<bool, int64_t>> const & input_times) {
    // Subtract the minimum valid start_time from all valid start_times
    // Leave invalid entries in start_time alone

    std::vector<std::pair<bool, int64_t>> start_times = input_times;

    bool found_min_start_time = false;
    int64_t min_start_time = 0;
    for(size_t i=0; i<start_times.size(); ++i) {
        if(not start_times[i].first)
            continue;
        if(not found_min_start_time) {
            min_start_time = start_times[i].second;
            found_min_start_time = true;
            continue;
        }
        if(start_times[i].second < min_start_time)
            min_start_time = start_times[i].second;
    }

    for(size_t i=0; i<start_times.size(); ++i)
        if(start_times[i].first)
            start_times[i] = {start_times[i].first, start_times[i].second - min_start_time};

    return {start_times, min_start_time};
}

I3FramePtr CCMTriggerMerger::TransformSingleTriggerFrame(I3FramePtr frame) {
    boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>> first_readout = frame->Get<boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>>>(digital_readout_name_);
    boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>> input_times = frame->Get<boost::shared_ptr<const I3Vector<std::pair<bool, int64_t>>>>(trigger_times_name_);
    boost::shared_ptr<I3Vector<CCMWaveformUInt16>> output_waveforms(new I3Vector<CCMWaveformUInt16>(first_readout->size()));

    std::tuple<std::vector<std::pair<bool, int64_t>>, int64_t> relative_start_times = ComputeRelativeStartTimes(*input_times);
    std::vector<std::pair<bool, int64_t>> start_times = std::get<0>(relative_start_times);
    int64_t absolute_start_time = std::get<1>(relative_start_times);

    // Copy raw waveform data into waveform object and set metadata accordingly
    for(size_t i=0; i<first_readout->size(); ++i) {
        CCMWaveformUInt16 & w = (*output_waveforms)[i];
        if((*first_readout)[i].size() > 0) {
            size_t board_idx = channel_idx_to_board_idx[i];
            double t = start_times[board_idx].second * ns_per_cycle;
            if(start_times[board_idx].first)
                w.SetStartTime(t);
            else
                w.SetStartTime(0);
            w.SetWaveform((*first_readout)[i]);
        }
        w.SetBinWidth(2);
        w.SetWaveformInformation({});
        w.SetSource(CCMSource::V1730);
    }
    frame->Delete(digital_readout_name_);
    frame->Put(first_trigger_time_output_, boost::make_shared<I3PODHolder<int64_t>>(absolute_start_time), I3Frame::DAQ);
    frame->Put(waveforms_output_, output_waveforms, I3Frame::DAQ);
    return frame;
}

I3FramePtr CCMTriggerMerger::TransformMultipleTriggerFrames(std::vector<I3FramePtr> const & frames) {
    I3FramePtr first_frame = frames[0];
    boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>> first_readout = first_frame->Get<boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>>>(digital_readout_name_);
    boost::shared_ptr<I3Vector<CCMWaveformUInt16>> output_waveforms(new I3Vector<CCMWaveformUInt16>());
    output_waveforms->reserve(first_readout->size());
    boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>> output_triggers(new I3Vector<CCMAnalysis::Binary::CCMTrigger>());
    output_triggers->reserve(frames.size());

    std::vector<std::pair<bool, int64_t>> start_times = GetStartTimes(frames);
    std::tuple<std::vector<std::pair<bool, int64_t>>, int64_t> relative_start_times = ComputeRelativeStartTimes(start_times);
    start_times = std::get<0>(relative_start_times);
    int64_t absolute_start_time = std::get<1>(relative_start_times) * ns_per_cycle;

    // Concatenate the output triggers
    for(size_t i=0; i<frames.size(); ++i) {
        output_triggers->emplace_back(frames[i]->Get<boost::shared_ptr<const I3Vector<CCMAnalysis::Binary::CCMTrigger>>>(triggers_name_)->operator[](0));
    }

    for(size_t w_idx=0; w_idx<first_readout->size(); ++w_idx) {
        // Get the waveform to save things to
        output_waveforms->emplace_back();
        CCMWaveformUInt16 & w = output_waveforms->back();

        // Compute the size of the waveform after merging
        size_t total_size = 0;
        for(size_t i=0; i<frames.size(); ++i) {
            total_size += frames[i]->Get<boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>>>(digital_readout_name_)->operator[](w_idx).size();
        }

        // Copy the input information into a single waveform
        std::vector<uint16_t> waveform(total_size);
        size_t last_pos = 0;
        for(size_t i=0; i<frames.size(); ++i) {
            I3Vector<uint16_t> const & input_waveform = frames[i]->Get<boost::shared_ptr<const I3Vector<I3Vector<uint16_t>>>>(digital_readout_name_)->operator[](w_idx);
            if(input_waveform.size() > 0) {
                if(start_times[w_idx].first)
                    w.SetStartTime(start_times[w_idx].second * ns_per_cycle);
                else
                    w.SetStartTime(0);
                std::copy(input_waveform.begin(), input_waveform.end(), waveform.begin() + last_pos);
                last_pos += input_waveform.size();
            }
        }

        // Assign all the information to the waveform object
        w.SetWaveform(waveform);
        w.SetBinWidth(2);
        w.SetWaveformInformation({});
        w.SetSource(CCMSource::V1730);
    }

    I3FramePtr output_frame(new I3Frame(I3Frame::DAQ));
    output_frame->Put(first_trigger_time_output_, boost::make_shared<I3PODHolder<int64_t>>(absolute_start_time), I3Frame::DAQ);
    output_frame->Put(waveforms_output_, output_waveforms, I3Frame::DAQ);
    output_frame->Put(triggers_name_, output_triggers, I3Frame::DAQ);

    return output_frame;
}

I3FramePtr CCMTriggerMerger::MergeTriggerFrames(std::vector<I3FramePtr> & frames) {
    if (frames.size() == 0) {
        return boost::shared_ptr<I3Frame>(new I3Frame(I3Frame::DAQ));
    } else if(frames.size() == 1) {
        return TransformSingleTriggerFrame(frames[0]);
    } else {
        return TransformMultipleTriggerFrames(frames);
    }
}

void CCMTriggerMerger::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() == I3Frame::Geometry) {
        CacheGeometryFrame(frame);
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
