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
        return frame->Has("CCMTriggerReadout") and frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers.size() > 0;
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
            if(not all_skipped)
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
    while(i < times0.size() or j < times1.size()) {
        if(i == times0.size() or j == times1.size()) {
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
            if((not last_i and last_i < i-1) or (not last_j and last_j < j-1)) {
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
    int64_t n_triggers = min_triggers;
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
            std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> pair_result = find_pairs(times1, times0, delta, max_delta);
            results.emplace_back(delta_trigger, delta, std::get<0>(pair_result).size(), std::get<1>(pair_result).size(), std::get<2>(pair_result).size());
            for(int64_t t : jitter_tests) {
                std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> pair_result = find_pairs(times1, times0, delta + t, max_delta);
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

std::vector<int64_t> compute_offsets(std::vector<std::vector<std::string>> file_lists, double max_time_diff) {
    int64_t max_delta = std::ceil(max_time_diff / 8.0);
}

