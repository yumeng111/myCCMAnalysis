#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <fstream>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>

#include "dataio/I3FileStager.h"
#include "dataio/I3FrameSequence.h"

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

int32_t subtract_times(uint32_t t1, uint32_t t0) {
    // Compute the difference of two unsigned time counters with one overflow bit
    // Cover the trivial case
    if(t1 == t0) {
        return 0;
    }

    // We know the 32nd bit is an overflow bit so we must mask it out
    bool overflow_0 = t0 & (0x1 << 31);
    bool overflow_1 = t1 & (0x1 << 31);
    bool different_overflow = overflow_0 ^ overflow_1;
    int32_t time_0 = uint32_t(t0 & 0x7FFFFFFF);
    int32_t time_1 = uint32_t(t1 & 0x7FFFFFFF);

    int32_t diff;
    if(different_overflow) {
        // We know which time is greater
        // But one time has crossed the overflow boundary
        if(overflow_0) {
            // Convert t1 to a negative number in the frame of t0's zero then subtract
            diff = (time_1 - 0x7FFFFFFE) - time_0;
        } else {
            // Convert t0 to a negative number in the frame of t1's zero then subtract
            diff = time_1 - int32_t(int32_t(time_0) - int32_t(0x7FFFFFFE));
        }
    } else {
        // We don't know which time is greater
        // Use the difference with the smallest magnitude

        // Cover the trivial case
        int32_t diff_0 = time_1 - time_0;

        // Cover the case where the difference crosses the overflow boundary
        int32_t diff_1;
        if(time_0 < time_1) {
            // Convert t1 to a negative number in the frame of t0's zero then subtract
            diff_1 = (time_1 - 0x7FFFFFFE) - time_0;
        } else {
            // Convert t0 to a negative number in the frame of t1's zero then subtract
            diff_1 = time_1 - (time_0 - 0x7FFFFFFE);
        }
        if(std::abs(diff_0) < std::abs(diff_1)) {
            diff = diff_0;
        } else {
            diff = diff_1;
        }
        // diff = std::min(diff_0, diff_1);
    }
    return diff;
}

std::vector<uint8_t> empty_mask(I3FramePtr frame) {
    CCMAnalysis::Binary::CCMDAQConfig const & config = frame->Get<CCMAnalysis::Binary::CCMDAQConfig>("CCMDAQConfig");
    std::vector<uint16_t> const & channel_sizes = frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers[0].channel_sizes;
    size_t n_boards = config.digitizer_boards.size();
    std::vector<uint8_t> mask(n_boards, 0);
    size_t last_idx = 0;
    for(size_t i=0; i<n_boards; ++i) {
        size_t n_channels = config.digitizer_boards[i].channels.size();
        size_t next_idx = last_idx + n_channels;
        for(size_t j=last_idx; j<next_idx; ++j) {
            mask[i] |= channel_sizes[j] > 0;
        }
        last_idx = next_idx;
    }
    return mask;
}

std::vector<uint32_t> read_times(I3FramePtr frame) {
    return frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers[0].board_times;
}

class TimeReader {
    dataio::I3FrameSequence frame_seq;
    I3FramePtr current_frame;
    size_t n_boards;
    std::vector<std::deque<int64_t>> time_cache;
    std::vector<uint32_t> last_raw_time;

    bool FrameMeetsRequirements(I3FramePtr frame) {
        return frame->Has("CCMTriggerReadout") and frame->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout")->triggers.size() > 0;
    }

    bool PopFrame() {
        if(not frame_seq.more())
            return false;
        current_frame = frame_seq.pop_daq();
        while(not FrameMeetsRequirements(current_frame)) {
            if(not frame_seq.more())
                return false;
            current_frame = frame_seq.pop_daq();
        }
        return true;
    }

    bool PopTimes() {
        bool result = PopFrame();
        if(not result)
            return false;
        std::vector<uint32_t> time_read = read_times(current_frame);
        std::vector<uint8_t> mask = empty_mask(current_frame);
        for(size_t i=0; i<n_boards; ++i) {
            if(not mask[i])
                continue;
            uint32_t raw_time = time_read[i];
            int64_t abs_time = time_cache[i].back() + subtract_times(raw_time, last_raw_time[i]);
            time_cache[i].push_back(abs_time);
            last_raw_time[i] = raw_time;
        }
        return true;
    }

public:
    TimeReader(std::vector<std::string> const & file_names, size_t n_skip) :
    frame_seq(file_names) {
        PopFrame();
        CCMAnalysis::Binary::CCMDAQConfig const & config = current_frame->Get<CCMAnalysis::Binary::CCMDAQConfig>("CCMDAQConfig");
        this->n_boards = config.digitizer_boards.size();
        time_cache.resize(n_boards);
        last_raw_time.resize(n_boards);
        std::vector<uint32_t> time_read = read_times(current_frame);
        std::vector<uint8_t> mask = empty_mask(current_frame);
        std::fill(last_raw_time.begin(), last_raw_time.end(), 0);
        for(size_t i=0; i< n_boards; ++i) {
            if(not mask[i])
                continue;
            time_cache[i].push_back(time_read[i]);
            last_raw_time[i] = time_read[i];
        }
        while(true) {
            bool all_skipped = true;
            for(size_t i=0; i<n_boards; ++i) {
                all_skipped &= time_cache[i].size() > n_skip;
            }
            if(all_skipped)
                break;
            PopTimes();
        }
        for(size_t i=0; i<n_boards; ++i) {
            size_t size = time_cache[i].size();
            size = std::min(size, n_skip);
            for(size_t j=0; j<size; ++j)
                time_cache[i].pop_front();
        }
    }
    TimeReader(std::vector<std::string> const & file_names) : TimeReader(file_names, 0) {}

    std::vector<int64_t> GetTimes(size_t N, size_t board_idx) {
        while(time_cache[board_idx].size() < N) {
            PopTimes();
        }
        std::vector<int64_t> result(N);
        std::copy(std::begin(time_cache[board_idx]), std::begin(time_cache[board_idx]) + N, std::begin(result));
        return result;
    }

    size_t NBoards() const {
        return this->n_boards;
    }
};

struct IDX {
    bool present = false;
    int64_t value = 0;
    IDX() {}
    IDX(IDX const & idx) : present(idx.present), value(idx.value) {}
    IDX(int64_t v) : present(true), value(v) {}
    IDX & operator=(int64_t v) {
        present = true;
        value = v;
        return *this;
    }
    IDX & operator=(IDX const & v) {
        present = v.present;
        value = v.value;
        return *this;
    }
    operator int64_t() {
        if(not present)
            throw std::runtime_error("No value present");
        return value;
    }
    operator bool() {
        return present;
    }
    bool operator<(int64_t v) {
        if(not present)
            throw std::runtime_error("No value present");
        return value < v;
    }
    bool operator>(int64_t v) {
        if(not present)
            throw std::runtime_error("No value present");
        return value > v;
    }
    IDX operator+(int64_t v) {
        if(not present)
            throw std::runtime_error("No value present");
        return IDX(value + v);
    }
};

std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> find_pairs(std::vector<int64_t> const & times0, std::vector<int64_t> const & times1, int64_t offset1, int64_t max_delta) {
    std::vector<std::tuple<IDX, IDX>> good_pairs;
    std::vector<std::tuple<IDX, IDX>> orphans;
    std::vector<std::tuple<IDX, IDX>> pairs;
    int64_t i = 0;
    int64_t j = 0;
    IDX last_i;
    IDX last_j;
    while(i < int64_t(times0.size()) or j < int64_t(times1.size())) {
        if(i == int64_t(times0.size()) or j == int64_t(times1.size())) {
            i = times0.size();
            j = times1.size();
            if(not last_i)
                last_i = -1;
            if(not last_j)
                last_j = -1;
            std::vector<std::tuple<IDX, IDX, int64_t>> triplets;
            for(int64_t _i = last_i + int64_t(1); _i<i; ++_i)
                triplets.emplace_back(_i, IDX(), times0[_i]);
            for(int64_t _j = last_j + int64_t(1); _j<j; ++_j)
                triplets.emplace_back(_j, IDX(), times1[_j] + offset1);
            for(size_t _i=0; _i<triplets.size(); ++_i) {
                orphans.emplace_back(std::get<0>(triplets[_i]), std::get<1>(triplets[_i]));
                pairs.emplace_back(std::get<0>(triplets[_i]), std::get<1>(triplets[_i]));
            }
            last_j = j;
            last_i = i;
            continue;
        }

        int64_t time_diff = times0[i] - (times1[j] + offset1);
        if(std::abs(time_diff) <= max_delta) {
            if((last_i and last_i < i-1) or (last_j and last_j < j-1)) {
                std::vector<std::tuple<IDX, IDX, int64_t>> triplets;
                for(int64_t _i = last_i + int64_t(1); _i<i; ++_i)
                    triplets.emplace_back(_i, IDX(), times0[_i]);
                for(int64_t _j = last_j + int64_t(1); _j<j; ++_j)
                    triplets.emplace_back(_j, IDX(), times1[_j] + offset1);
                std::sort(triplets.begin(), triplets.end(), [](std::tuple<IDX, IDX, int64_t> const & a, std::tuple<IDX, IDX, int64_t> const & b){return std::get<2>(a) < std::get<2>(b);});
                for(size_t _i=0; _i<triplets.size(); ++_i) {
                    orphans.emplace_back(std::get<0>(triplets[_i]), std::get<1>(triplets[_i]));
                    pairs.emplace_back(std::get<0>(triplets[_i]), std::get<1>(triplets[_i]));
                }
            }
            pairs.emplace_back(i, j);
            good_pairs.emplace_back(pairs.back());
            last_i = i;
            last_j = j;
            i += 1;
            j += 1;
        } else {
            if(times0[i] < times1[i] + offset1) {
                i += 1;
            } else {
                j += 1;
            }
        }
    }
    return {pairs, good_pairs, orphans};
}

int64_t get_time_delta(std::vector<int64_t> times0, std::vector<int64_t> times1, size_t delta_trigger) {
    if(delta_trigger < 0)
        return times0[0] - times1[-delta_trigger];
    else
        return times0[delta_trigger] - times1[0];
}

int64_t compute_trigger_offset(TimeReader & reader0, size_t board_idx0, TimeReader & reader1, size_t board_idx1, std::vector<int64_t> jitter_tests={-2, 2}, int64_t max_delta=2, size_t min_triggers=100, size_t max_triggers=2000, size_t increment=25, double threshold=0.9) {
    uint64_t n_triggers = min_triggers;
    int64_t prev_min = 0;
    int64_t prev_max = 0;
    int64_t new_min = -(n_triggers/2);
    int64_t new_max = (n_triggers/2) + 1;
    std::vector<int64_t> times0 = reader0.GetTimes(n_triggers, board_idx0);
    std::vector<int64_t> times1 = reader1.GetTimes(n_triggers, board_idx1);
    size_t x = prev_max;
    std::vector<int64_t> deltas;
    std::generate_n(std::back_inserter(deltas), new_max - prev_max, [&](){return x++;});
    x = new_min;
    std::generate_n(std::back_inserter(deltas), prev_min - new_min, [&](){return x++;});
    std::vector<std::tuple<int64_t, int64_t, size_t, size_t, size_t>> results;
    std::tuple<int64_t, int64_t, size_t, size_t, size_t> best_pair_result;

    while(true) {
        for(int64_t const & delta_trigger : deltas) {
            int64_t delta = get_time_delta(times0, times1, delta_trigger);
            std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> pair_result = find_pairs(times0, times1, delta, max_delta);
            results.emplace_back(delta_trigger, delta, std::get<0>(pair_result).size(), std::get<1>(pair_result).size(), std::get<2>(pair_result).size());
            for(int64_t t : jitter_tests) {
                std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> pair_result = find_pairs(times0, times1, delta + t, max_delta);
                results.emplace_back(delta_trigger, delta + t, std::get<0>(pair_result).size(), std::get<1>(pair_result).size(), std::get<2>(pair_result).size());
            }
        }

        best_pair_result = *std::max_element(results.begin(), results.end(), [](std::tuple<int64_t, int64_t, size_t, size_t, size_t> const & a, std::tuple<int64_t, int64_t, size_t, size_t, size_t> const & b){return std::get<3>(a) < std::get<3>(b);});

        int64_t delta_trigger = std::get<0>(best_pair_result);
        int64_t delta = std::get<1>(best_pair_result);
        size_t pairs = std::get<2>(best_pair_result);
        size_t good_pairs = std::get<3>(best_pair_result);
        size_t orphans = std::get<4>(best_pair_result);

        if(double(good_pairs) / double(n_triggers - std::abs(delta_trigger)) >= threshold)
            break;

        n_triggers += increment;
        if(n_triggers > max_triggers)
            throw std::runtime_error("Could not align triggers");

        new_max = int64_t(0.5*n_triggers) + 1;
        new_min = -int64_t(0.5*n_triggers);
        deltas.resize(0);
        x = prev_max;
        std::generate_n(std::back_inserter(deltas), new_max - prev_max, [&](){return x++;});
        x = new_min;
        std::generate_n(std::back_inserter(deltas), prev_min - new_min, [&](){return x++;});

        times0 = reader0.GetTimes(n_triggers, board_idx0);
        times1 = reader1.GetTimes(n_triggers, board_idx1);

        prev_min = new_min;
        prev_max = new_max;
    }

    return std::get<1>(best_pair_result);
}

std::vector<std::vector<int64_t>> compute_offsets(std::vector<std::vector<std::string>> file_lists, double max_time_diff) {
    int64_t max_delta = std::ceil(max_time_diff / 8.0);
    size_t n_daqs = file_lists.size();
    std::vector<std::vector<int64_t>> offsets(n_daqs);
    std::vector<size_t> n_boards;
    std::vector<TimeReader> readers;
    readers.reserve(n_daqs);
    for(size_t i=0; i<n_daqs; ++i) {
        readers.emplace_back(file_lists[i]);
        n_boards.emplace_back(readers.back().NBoards());
        offsets[i] = std::vector<int64_t>(n_boards.back(), 0);
    }
    for(size_t i=1; i<n_daqs; ++i)
        offsets[i][0] = compute_trigger_offset(readers[0], 0, readers[i], 0);
    for(size_t i=0; i<n_daqs; ++i)
        for(size_t j=1; j<n_boards[i]; ++j)
            offsets[i][j] = compute_trigger_offset(readers[0], 0, readers[i], j);
    return offsets;
}

class MergedSource : public I3Module {
    size_t max_delta;
    size_t n_daqs;
    std::vector<size_t> n_boards;
    std::vector<dataio::I3FrameSequencePtr> frame_sequences;
    std::vector<std::deque<I3FramePtr>> frame_cache;
    std::vector<std::vector<int64_t>> frame_idxs;
    std::vector<std::vector<uint8_t>> empty;
    std::vector<std::vector<std::deque<uint8_t>>> mask_cache;
    std::vector<std::vector<std::deque<long double>>> time_cache;
    std::vector<std::vector<uint32_t>> last_raw_time;
    std::vector<std::vector<int64_t>> last_time;
    std::vector<CCMAnalysis::Binary::CCMDAQConfigConstPtr> configs;
    int fill_computer_time = -1;
    bool push_config = false;

    std::vector<std::vector<int64_t>> offsets;

    size_t counter = 0;

    bool PopFrame(size_t daq_idx);
    I3FramePtr GetConfigFrame();
    bool NextTrigger(size_t daq_idx, size_t board_idx);
    bool NextTriggers();
    void ClearUnusedFrames();
    std::tuple<boost::shared_ptr<I3Vector<I3Vector<uint16_t>>>, boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>>> GetTriggerReadout();
public:
    MergedSource(const I3Context&);
    void Configure();
    void Process();

    SET_LOGGER("MergedSource");
};

I3_MODULE(MergedSource);

bool MergedSource::PopFrame(size_t daq_idx) {
    I3FramePtr frame;
    while(true) {
        if(not frame_sequences[daq_idx]->more())
            return false;
        frame = frame_sequences[daq_idx]->pop_frame();
        if(frame->GetStop() == I3Frame::Geometry) {
            std::cerr << "Setting config for " << daq_idx << std::endl;
            configs[daq_idx] = frame->Get<CCMAnalysis::Binary::CCMDAQConfigConstPtr>("CCMDAQConfig");
            push_config = true;
            continue;
        } else if(frame->GetStop() != I3Frame::DAQ) {
            continue;
        }
        if((not frame->Has("CCMTriggerReadout")) or (frame->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout")->triggers.size() == 0))
            continue;
        break;
    }
    frame_cache[daq_idx].push_back(frame);
    CCMAnalysis::Binary::CCMTriggerReadoutConstPtr readout = frame->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout");
    std::vector<uint32_t> const & time_read = readout->triggers[0].board_times;
    if(fill_computer_time < 0)
        fill_computer_time = readout->triggers[0].board_computer_times.size() > 0;
    std::vector<uint8_t> mask = empty_mask(frame);
    for(size_t board_idx=0; board_idx<n_boards[daq_idx]; ++board_idx) {
        mask_cache[daq_idx][board_idx].push_back(mask[board_idx]);
        if(not mask[board_idx]) {
            time_cache[daq_idx][board_idx].push_back(std::numeric_limits<long double>::infinity());
            continue;
        }
        uint32_t raw_time = time_read[board_idx];
        int64_t rel_time = subtract_times(raw_time, last_raw_time[daq_idx][board_idx]);
        long double abs_time = rel_time;
        if(time_cache[daq_idx][board_idx].size() > 0)
            abs_time += last_time[daq_idx][board_idx];
        time_cache[daq_idx][board_idx].push_back(abs_time);
        last_raw_time[daq_idx][board_idx] = raw_time;
        last_time[daq_idx][board_idx] = abs_time;
    }
    return true;
}

I3FramePtr MergedSource::GetConfigFrame() {
    CCMAnalysis::Binary::CCMDAQConfigPtr config = boost::make_shared<CCMAnalysis::Binary::CCMDAQConfig>();
    for(CCMAnalysis::Binary::CCMDAQConfigConstPtr c : configs) {
        std::copy(c->machine_configurations.begin(), c->machine_configurations.end(), std::back_inserter(config->machine_configurations));
        std::copy(c->digitizer_boards.begin(), c->digitizer_boards.end(), std::back_inserter(config->digitizer_boards));
    }
    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::Geometry);
    frame->Put("CCMDAQConfig", config, I3Frame::Geometry);
    return frame;
}

bool MergedSource::NextTrigger(size_t daq_idx, size_t board_idx) {
    size_t current_frame_idx = frame_idxs[daq_idx][board_idx];
    size_t frame_idx = current_frame_idx + 1;
    while(true) {
        while(frame_cache[daq_idx].size() <= frame_idx) {
            bool res = PopFrame(daq_idx);
            if(not res)
                return false;
        }
        if(mask_cache[daq_idx][board_idx][frame_idx])
            break;
        frame_idx += 1;
    }
    frame_idxs[daq_idx][board_idx] = frame_idx;
    return true;
}

bool MergedSource::NextTriggers() {
    for(size_t daq_idx=0; daq_idx < n_daqs; ++daq_idx) {
        for(size_t board_idx=0; board_idx < n_boards[daq_idx]; ++board_idx) {
            bool res = NextTrigger(daq_idx, board_idx);
            if(not res)
                return false;
        }
    }
    ClearUnusedFrames();
    return true;
}

void MergedSource::ClearUnusedFrames() {
    for(size_t daq_idx=0; daq_idx < n_daqs; ++daq_idx) {
        size_t min_idx = *std::min_element(std::begin(frame_idxs[daq_idx]), std::end(frame_idxs[daq_idx]));
        for(size_t i=0; i<min_idx; ++i) {
            frame_cache[daq_idx].pop_front();
            for(size_t board_idx=0; board_idx < n_boards[daq_idx]; ++board_idx) {
                mask_cache[daq_idx][board_idx].pop_front();
                time_cache[daq_idx][board_idx].pop_front();
                frame_idxs[daq_idx][board_idx] -= 1;
            }
        }
    }
}

inline void merge_triggers(
        boost::shared_ptr<I3Vector<I3Vector<uint16_t>>> output_samples,
        boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>> output_triggers,
        CCMAnalysis::Binary::CCMTriggerReadoutConstPtr source_trigger,
        size_t board_idx,
        size_t first_idx,
        size_t last_idx,
        bool fill_computer_time = false) {
    CCMAnalysis::Binary::CCMTrigger & o = (*output_triggers)[0];
    CCMAnalysis::Binary::CCMTrigger const & s = source_trigger->triggers[0];
    assert(s.channel_sizes.size() > first_idx and s.channel_sizes.size() >= last_idx);
    assert(s.channel_masks.size() > first_idx and s.channel_masks.size() >= last_idx);
    assert(s.channel_temperatures.size() > first_idx and s.channel_sizes.size() >= last_idx);
    assert(s.board_event_numbers.size() > board_idx);
    assert(s.board_times.size() > board_idx);
    std::copy(s.channel_sizes.begin() + first_idx, s.channel_sizes.begin() + last_idx, std::back_inserter(o.channel_sizes));
    std::copy(s.channel_masks.begin() + first_idx, s.channel_masks.begin() + last_idx, std::back_inserter(o.channel_masks));
    std::copy(s.channel_temperatures.begin() + first_idx, s.channel_temperatures.begin() + last_idx, std::back_inserter(o.channel_temperatures));
    std::copy(s.board_event_numbers.begin() + board_idx, s.board_event_numbers.begin() + board_idx + 1, std::back_inserter(o.board_event_numbers));
    std::copy(s.board_times.begin() + board_idx, s.board_times.begin() + board_idx + 1, std::back_inserter(o.board_times));
    if(fill_computer_time) {
        assert(s.board_computer_times.size() > board_idx);
        std::copy(s.board_computer_times.begin() + board_idx, s.board_computer_times.begin() + board_idx + 1, std::back_inserter(o.board_computer_times));
    }
    if(source_trigger->samples.size() == s.channel_masks.size()) {
        assert(source_trigger->samples.size() > first_idx and source_trigger->samples.size() >= last_idx);
        size_t n_channels = last_idx - first_idx;
        output_samples->reserve(output_samples->size() + n_channels);
        for(size_t i=first_idx; i<last_idx; ++i) {
            output_samples->emplace_back(source_trigger->samples[i]);
        }
        // std::copy(source_trigger->samples.begin() + first_idx, source_trigger->samples.begin() + last_idx, std::back_inserter(*static_cast<boost::shared_ptr<std::vector<I3Vector<uint16_t>>>>(output_samples)));
    } else {
        size_t new_first_idx = 0;
        for(size_t i=0; i<first_idx; ++i) {
            if(s.channel_masks[i])
                ++new_first_idx;
        }
        size_t new_last_idx = new_first_idx + last_idx - first_idx;
        assert(source_trigger->samples.size() > new_first_idx and source_trigger->samples.size() >= new_last_idx);
        size_t n_channels = new_last_idx - new_first_idx;
        output_samples->reserve(output_samples->size() + n_channels);
        for(size_t i=new_first_idx; i<new_last_idx; ++i) {
            output_samples->emplace_back(source_trigger->samples[i]);
        }
        // std::copy(source_trigger->samples.begin() + new_first_idx, source_trigger->samples.begin() + new_last_idx, std::back_inserter(*static_cast<boost::shared_ptr<std::vector<I3Vector<uint16_t>>>>(output_samples)));
    }
}

inline void merge_empty_trigger(
        boost::shared_ptr<I3Vector<I3Vector<uint16_t>>> output_samples,
        boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>> output_triggers,
        size_t n_channels,
        bool fill_computer_time = false) {
    CCMAnalysis::Binary::CCMTrigger & o = (*output_triggers)[0];
    std::fill_n(std::back_inserter(o.channel_sizes), n_channels, 0);
    std::fill_n(std::back_inserter(o.channel_masks), n_channels, 0);
    std::fill_n(std::back_inserter(o.channel_temperatures), n_channels, 0);
    output_samples->reserve(output_samples->size() + n_channels);
    for(size_t i=0; i<n_channels; ++i) {
        output_samples->emplace_back();
    }
    std::fill_n(std::back_inserter(o.board_event_numbers), 1, 0);
    std::fill_n(std::back_inserter(o.board_times), 1, 0);
    if(fill_computer_time) {
        struct timespec t; t.tv_sec = 0; t.tv_nsec = 0;
        std::fill_n(std::back_inserter(o.board_computer_times), 1, t);
    }
}

std::tuple<boost::shared_ptr<I3Vector<I3Vector<uint16_t>>>, boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>>>
    MergedSource::GetTriggerReadout() {
    CCMAnalysis::Binary::CCMTriggerReadoutPtr readout = boost::make_shared<CCMAnalysis::Binary::CCMTriggerReadout>();
    boost::shared_ptr<I3Vector<I3Vector<uint16_t>>> output_samples = boost::make_shared<I3Vector<I3Vector<uint16_t>>>();
    boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>> output_triggers = boost::make_shared<I3Vector<CCMAnalysis::Binary::CCMTrigger>>(1);
    std::vector<std::vector<long double>> times;
    times.resize(n_daqs);
    long double min_time = std::numeric_limits<long double>::infinity();
    bool all_bad = true;
    for(size_t daq_idx=0; daq_idx < n_daqs; ++daq_idx) {
        std::vector<long double> & t = times[daq_idx];
        t.resize(n_boards[daq_idx]);
        for(size_t board_idx=0; board_idx < n_boards[daq_idx]; ++board_idx) {
            size_t frame_idx = frame_idxs[daq_idx][board_idx];
            if(empty[daq_idx][board_idx]) {
                t[board_idx] = std::numeric_limits<long double>::infinity();
            } else {
                t[board_idx] = time_cache[daq_idx][board_idx][frame_idx] + offsets[daq_idx][board_idx];
            }
            min_time = std::min(min_time, t[board_idx]);
            all_bad &= empty[daq_idx][board_idx];
        }
    }

    if(all_bad or std::isinf(min_time))
        return {nullptr, nullptr};
    std::vector<std::vector<int>> state(n_daqs);
    bool print_empty = false;
    for(size_t daq_idx=0; daq_idx < n_daqs; ++daq_idx) {
        size_t last_idx = 0;
        for(size_t board_idx=0; board_idx < n_boards[daq_idx]; ++board_idx) {
            size_t n_channels = configs[daq_idx]->digitizer_boards[board_idx].channels.size();
            // channel_sizes.push_back(n_channels);
            size_t next_idx = last_idx + n_channels;
            size_t frame_idx = frame_idxs[daq_idx][board_idx];
            if(empty[daq_idx][board_idx]) {
                print_empty = true;
                merge_empty_trigger(output_samples, output_triggers, n_channels, fill_computer_time);
                state[daq_idx].push_back(0);
            } else if(std::isinf(times[daq_idx][board_idx])) {
                print_empty = true;
                state[daq_idx].push_back(1);
                throw std::runtime_error("Should not see inf here!");
            } else if(times[daq_idx][board_idx] - min_time > max_delta) {
                print_empty = true;
                merge_empty_trigger(output_samples, output_triggers, n_channels, fill_computer_time);
                state[daq_idx].push_back(2);
            } else {
                CCMAnalysis::Binary::CCMTriggerReadoutConstPtr tr = frame_cache[daq_idx][frame_idx]->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout");
                merge_triggers(output_samples, output_triggers, tr, board_idx, last_idx, next_idx, fill_computer_time);
                bool res = NextTrigger(daq_idx, board_idx);
                if(not res) {
                    times[daq_idx][board_idx] = std::numeric_limits<long double>::infinity();
                    empty[daq_idx][board_idx] = true;
                    print_empty = true;
                    state[daq_idx].push_back(4);
                } else {
                    state[daq_idx].push_back(3);
                }
            }
            last_idx = next_idx;
        }
    }
    if(print_empty) {
        std::cerr << "Min time: " << min_time << std::endl;
        std::cerr << "Empty : [" << std::endl;
        for(size_t daq_idx=0; daq_idx < n_daqs; ++daq_idx) {
            std::cerr << "\t[";
            for(size_t board_idx=0; board_idx < n_boards[daq_idx]; ++board_idx) {
                std::cerr << int(state[daq_idx][board_idx]) << " ";
            }
            std::cerr << "]," << std::endl;
        }
        std::cerr << "]" << std::endl;
    }
    ClearUnusedFrames();
    return {output_samples, output_triggers};
}

MergedSource::MergedSource(const I3Context& context) : I3Module(context) {
    AddParameter("FileLists",
            "File lists to merge",
            std::vector<std::vector<std::string>>());

    AddParameter("MaxTimeDiff",
            "Maximum time difference between associated triggers (ns)",
            16);

}

void MergedSource::Configure() {
    std::vector<std::vector<std::string>> file_lists;
    GetParameter("FileLists", file_lists);

    double max_time_diff;
    GetParameter("MaxTimeDiff", max_time_diff);
    max_delta = std::ceil(max_time_diff / 8.0);

    offsets = compute_offsets(file_lists, max_time_diff);

    n_daqs = file_lists.size();
    frame_sequences.reserve(n_daqs);

    n_boards.resize(n_daqs);
    frame_cache.resize(n_daqs);
    frame_idxs.resize(n_daqs);
    empty.resize(n_daqs);
    mask_cache.resize(n_daqs);
    time_cache.resize(n_daqs);
    last_raw_time.resize(n_daqs);
    last_time.resize(n_daqs);
    configs.resize(n_daqs);

    for(size_t daq_idx=0; daq_idx<n_daqs; ++daq_idx) {
        frame_sequences.emplace_back(boost::make_shared<dataio::I3FrameSequence>(file_lists[daq_idx]));
        size_t n = offsets[daq_idx].size();
        n_boards[daq_idx] = n;
        frame_idxs[daq_idx] = std::vector<int64_t>(n, -1);
        empty[daq_idx] = std::vector<uint8_t>(n, 0);
        last_raw_time[daq_idx] = std::vector<uint32_t>(n, 0);
        last_time[daq_idx] = std::vector<int64_t>(n, 0);
        for(size_t j=0; j<n; ++j) {
            mask_cache[daq_idx].emplace_back();
            time_cache[daq_idx].emplace_back();
        }
    }

    NextTriggers();
}

void MergedSource::Process() {
    if(push_config) {
        PushFrame(GetConfigFrame());
        push_config = false;
    }

    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::DAQ);
    std::tuple<boost::shared_ptr<I3Vector<I3Vector<uint16_t>>>, boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>>> readout = GetTriggerReadout();
    if(std::get<0>(readout) == nullptr or std::get<1>(readout) == nullptr) {
        RequestSuspension();
        return;
    }
    frame->Put("CCMDigitalReadout", std::get<0>(readout), I3Frame::DAQ);
    frame->Put("CCMTriggers", std::get<1>(readout), I3Frame::DAQ);
    PushFrame(frame);
    std::cerr << "Pushed frame " << counter << std::endl;
    ++counter;
}
