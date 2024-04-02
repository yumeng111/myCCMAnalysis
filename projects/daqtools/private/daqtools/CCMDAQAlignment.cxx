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
#include <algorithm>

#include <ctime>
#include <cstdio>
#include <string>
#include <time.h>
#include <cstdlib>
#include <cstring>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>

#include "dataio/I3FileStager.h"
#include "dataio/I3FrameSequence.h"

#include "dataclasses/physics/CCMEventHeader.h"
#include "dataclasses/I3Time.h"

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

namespace {
std::string GetTimeStringToSecond(time_t now) {
    char buf[sizeof "2011-10-08T07:07:09Z"];
    std::strftime(buf, sizeof buf, "%FT%TZ", gmtime(&now));
    return std::string(buf);
}
}

inline int ParseInt(const char* value) {
    return std::strtol(value, nullptr, 10);
}

struct timespec ParseISO8601(std::string const & input) {
    struct timespec precision_time;
    precision_time.tv_sec = 0;
    precision_time.tv_nsec = 0;
    constexpr const size_t expected_length = sizeof("1234-12-12T12:12:12Z") - 1;
    static_assert(expected_length == 20, "Unexpected ISO 8601 date/time length");

    if (input.size() < expected_length) {
        return precision_time;
    }

    size_t size = input.size();
    char * c_str = new char[size + 1];
    std::memcpy(c_str, input.c_str(), size + 1);

    std::tm time = { 0 };
    time.tm_year = ParseInt(&c_str[0]) - 1900;
    time.tm_mon = ParseInt(&c_str[5]) - 1;
    time.tm_mday = ParseInt(&c_str[8]);
    time.tm_hour = ParseInt(&c_str[11]);
    time.tm_min = ParseInt(&c_str[14]);
    time.tm_sec = ParseInt(&c_str[17]);
    time.tm_isdst = 0;
    precision_time.tv_sec = timegm(&time);
    if(input.size() > expected_length) {
        char * start = &c_str[20];
        char * end_pos;
        long sub_second_count = std::strtol(start, &end_pos, 10);
        size_t sub_second_len = std::distance(start, end_pos);
        if(sub_second_count > 9) {
            for(unsigned int i=0; i<sub_second_len-9; ++i)
                sub_second_count /= 10;
        } else if(sub_second_count < 9) {
            for(unsigned int i=0; i<9-sub_second_len; ++i)
                sub_second_count *= 10;
        }
        precision_time.tv_nsec = sub_second_count;
    }
    delete[] c_str;
    return precision_time;
}

bool ParseFileName(std::string fname, int & run_number, int & file_number, struct timespec & spec) {
    std::string parse = fname;

    std::string delimiter = "/";
    size_t char_pos = parse.rfind(delimiter);
    size_t size = char_pos;
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" folder could not be parsed.", fname.c_str());
    } else {
        size = char_pos;
        parse = parse.substr(size + delimiter.size(), std::string::npos);
    }

    delimiter = "_";
    char_pos = parse.find(delimiter);
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" prefix could not be parsed.", fname.c_str());
        return false;
    }
    size = char_pos;
    std::string prefix = parse.substr(0, size);
    parse = parse.substr(size + delimiter.size(), std::string::npos);

    char_pos = parse.find(delimiter);
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" run_number could not be parsed.", fname.c_str());
        return false;
    }
    size = char_pos;
    std::string run_string = parse.substr(0, size);
    parse = parse.substr(size + delimiter.size(), std::string::npos);

    delimiter = "run";
    char_pos = run_string.find(delimiter);
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" run_number could not be parsed.", fname.c_str());
        return false;
    }
    size = char_pos;
    run_string = run_string.substr(size + delimiter.size(), std::string::npos);
    run_number = std::atoi(run_string.c_str());

    delimiter = "_";
    char_pos = parse.find(delimiter);
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" file_number could not be parsed.", fname.c_str());
        return false;
    }
    size = char_pos;
    std::string file_string = parse.substr(0, size);
    parse = parse.substr(size + delimiter.size(), std::string::npos);

    delimiter = "file";
    char_pos = file_string.find(delimiter);
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" file_number could not be parsed.", fname.c_str());
        return false;
    }
    size = char_pos;
    file_string = file_string.substr(size + delimiter.size(), std::string::npos);
    file_number = std::atoi(file_string.c_str());

    delimiter = ".";
    char_pos = parse.find(delimiter);
    if(char_pos == std::string::npos) {
        log_warn("Filename \"%s\" time could not be parsed.", fname.c_str());
        return false;
    }
    size = char_pos;
    std::string iso_string = parse.substr(0, size);
    parse = parse.substr(size + delimiter.size(), std::string::npos);
    spec = ParseISO8601(iso_string);

    return true;
}

std::vector<uint8_t> empty_mask(I3FramePtr frame) {
    // Create a mask that ignores empty triggers

    CCMAnalysis::Binary::CCMDAQConfig const & config = frame->Get<CCMAnalysis::Binary::CCMDAQConfig>("CCMDAQConfig");
    std::vector<uint16_t> const & channel_sizes = frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers[0].channel_sizes;
    std::vector<uint16_t> const & channel_masks = frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers[0].channel_masks;
    size_t n_boards = config.digitizer_boards.size();
    std::vector<uint8_t> mask(n_boards, 0);
    size_t last_idx = 0;
    for(size_t i=0; i<n_boards; ++i) {
        size_t n_channels = config.digitizer_boards[i].channels.size();
        size_t next_idx = last_idx + n_channels;
        for(size_t j=last_idx; j<next_idx; ++j) {
            // An empty trigger is marked by having a channel size of zero or a channel mask of zero
            bool not_empty = channel_sizes[j] > 0 and channel_masks[j] > 0;
            mask[i] |= not_empty;
        }
        last_idx = next_idx;
    }
    return mask;
}

std::vector<uint32_t> read_times(I3FramePtr frame) {
    // Get the raw board times from the frame
    return frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers[0].board_times;
}

std::vector<struct timespec> read_computer_times(I3FramePtr frame) {
    // Get the computer times from the frame
    return frame->Get<CCMAnalysis::Binary::CCMTriggerReadout>("CCMTriggerReadout").triggers[0].board_computer_times;
}

class TimeReader {
    // Frame sequence to read from
    dataio::I3FrameSequence frame_seq;
    // Currently loaded frame
    I3FramePtr current_frame;
    // Number of boards
    size_t n_boards;
    // Cache of read times
    std::vector<std::deque<int64_t>> time_cache;
    std::vector<std::deque<struct timespec>> computer_time_cache;
    // Last raw time read
    std::vector<uint32_t> last_raw_time;

    bool FrameMeetsRequirements(I3FramePtr frame) {
        // Check if the frame has the necessary information for merging
        return frame->Has("CCMTriggerReadout") and frame->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout")->triggers.size() > 0;
    }

    bool PopFrame() {
        // Get another frame and cache it
        // Skip frames if they do not meet requirements
        // Return true for a found frame
        // Return false if no more frames are available
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
        // Pop a frame and put the times into the time cache
        bool result = PopFrame();
        if(not result)
            return false;
        std::vector<uint32_t> time_read = read_times(current_frame);
        std::vector<struct timespec> computer_time_read = read_computer_times(current_frame);
        std::vector<uint8_t> mask = empty_mask(current_frame);
        for(size_t i=0; i<n_boards; ++i) {
            if(not mask[i]) // Empty triggers have no time associated so we can skip placing anything in the cache
                continue;
            uint32_t raw_time = time_read[i];
            int64_t abs_time = time_cache[i].back() + CCMAnalysis::Binary::subtract_times(raw_time, last_raw_time[i]);
            time_cache[i].push_back(abs_time);
            computer_time_cache[i].push_back(computer_time_read[i]);
            last_raw_time[i] = raw_time;
        }
        return true;
    }

public:
    TimeReader(std::vector<std::string> const & file_names, size_t n_skip) :
    frame_seq(file_names) {
        // Get the first frame
        PopFrame();
        // Assume we have the configuration from the first frame
        CCMAnalysis::Binary::CCMDAQConfig const & config = current_frame->Get<CCMAnalysis::Binary::CCMDAQConfig>("CCMDAQConfig");

        // Set up the containers and fill the cache
        this->n_boards = config.digitizer_boards.size();
        time_cache.resize(n_boards);
        computer_time_cache.resize(n_boards);
        last_raw_time.resize(n_boards);
        std::vector<uint32_t> time_read = read_times(current_frame);
        std::vector<struct timespec> computer_time_read = read_computer_times(current_frame);
        std::vector<uint8_t> mask = empty_mask(current_frame);
        std::fill(last_raw_time.begin(), last_raw_time.end(), 0);
        for(size_t i=0; i< n_boards; ++i) {
            if(not mask[i])
                continue;
            time_cache[i].push_back(time_read[i]);
            computer_time_cache[i].push_back(computer_time_read[i]);
            last_raw_time[i] = time_read[i];
        }

        // Load more times until we have at least the number we intend to skip in each cache
        while(true) {
            bool all_skipped = true;
            for(size_t i=0; i<n_boards; ++i) {
                all_skipped &= time_cache[i].size() > n_skip;
            }
            if(all_skipped)
                break;
            bool res = PopTimes();
            if(not res)
                log_fatal("Not enough frames to align data streams");
        }

        // Skip times on each board
        for(size_t i=0; i<n_boards; ++i) {
            size_t size = time_cache[i].size();
            size = std::min(size, n_skip);
            for(size_t j=0; j<size; ++j)
                time_cache[i].pop_front();
        }
    }

    // By default we should not skip any timed
    // Skipping is for testing purposes only
    TimeReader(std::vector<std::string> const & file_names) : TimeReader(file_names, 0) {}

    std::vector<int64_t> GetTimes(size_t N, size_t board_idx) {
        // Get the specified number of times for the chosen board
        // Pop frames and store times until we have enough
        while(time_cache[board_idx].size() < N) {
            bool res = PopTimes();
            if(not res)
                log_fatal("Not enough frames to align data streams");
        }
        // Copy and return the resulting times
        std::vector<int64_t> result(N);
        std::copy(std::begin(time_cache[board_idx]), std::begin(time_cache[board_idx]) + N, std::begin(result));
        return result;
    }

    std::vector<std::deque<int64_t>> & GetAllTimes() {
        return time_cache;
    }

    std::vector<std::deque<struct timespec>> & GetAllComputerTimes() {
        return computer_time_cache;
    }

    size_t NBoards() const {
        // Return the number of boards
        return this->n_boards;
    }
};

struct IDX {
    // Represents an index and replicates python-style missing data with a "None" value
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
    // Attempts to pair up times from two streams that match within a specified tolerance
    // Paired times are represented as a pair of two value indices
    std::vector<std::tuple<IDX, IDX>> good_pairs;
    // Unpaired times are represented by a pair with one "missing data" index and one valid index
    std::vector<std::tuple<IDX, IDX>> orphans;
    // All paired and unpaired times are collected in the following structure
    std::vector<std::tuple<IDX, IDX>> pairs;

    // Indicated the current position being processed in each time stream
    int64_t i = 0;
    int64_t j = 0;

    // Indicates the last position saved as a pair or orphan in each time stream
    // An invalid index indicates no times in the corresponding stream have been saved
    IDX last_i;
    IDX last_j;
    while(i < int64_t(times0.size()) or j < int64_t(times1.size())) {
        // If either time stream have no times to process, process the remaining times in the other stream as orphans
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

        // Compute the time difference between the two currently selected times (including the offset)
        int64_t time_diff = times0[i] - (times1[j] + offset1);

        // Check if the time difference is within the matching threshold
        if(std::abs(time_diff) <= max_delta) {
            // If the last saved index is invalid or the last saved index is less than the previously processed index
            //   then there are orphaned times to save
            if((last_i and last_i < i-1) or (last_j and last_j < j-1)) {
                // Grab the orphaned indices and their times
                std::vector<std::tuple<IDX, IDX, int64_t>> triplets;
                for(int64_t _i = last_i + int64_t(1); _i<i; ++_i)
                    triplets.emplace_back(_i, IDX(), times0[_i]);
                for(int64_t _j = last_j + int64_t(1); _j<j; ++_j)
                    triplets.emplace_back(_j, IDX(), times1[_j] + offset1);
                // Sort the orphaned indices by their times
                std::sort(triplets.begin(), triplets.end(), [](std::tuple<IDX, IDX, int64_t> const & a, std::tuple<IDX, IDX, int64_t> const & b){return std::get<2>(a) < std::get<2>(b);});
                // Add the orphaned indices to the relevant structures
                for(size_t _i=0; _i<triplets.size(); ++_i) {
                    orphans.emplace_back(std::get<0>(triplets[_i]), std::get<1>(triplets[_i]));
                    pairs.emplace_back(std::get<0>(triplets[_i]), std::get<1>(triplets[_i]));
                }
            }
            // Store the matching pair and increment the appropriate counters
            pairs.emplace_back(i, j);
            good_pairs.emplace_back(pairs.back());
            last_i = i;
            last_j = j;
            i += 1;
            j += 1;
        } else {
            // Increment the appropriate counter to try and find a matching time
            if(times0[i] < times1[i] + offset1) {
                i += 1;
            } else {
                j += 1;
            }
        }
    }

    // Return the matched and unmatched indices
    return {pairs, good_pairs, orphans};
}

int64_t get_time_delta(std::vector<int64_t> times0, std::vector<int64_t> times1, size_t delta_trigger) {
    // Get the time difference based on a signed index offset between the two streams
    if(delta_trigger < 0)
        return times0[0] - times1[-delta_trigger];
    else
        return times0[delta_trigger] - times1[0];
}

struct timespec subtract_trigger_timediff_from_timespec(struct timespec spec, int64_t delta_trigger) {
    constexpr long double ns_per_cycle = 8;
    if(delta_trigger < 0)
        return spec;
    struct timespec result = spec;
    for(size_t i=0; i<ns_per_cycle; ++i) {
        int64_t remaining = delta_trigger;

        if(remaining > 1e9) {
            int64_t diff = delta_trigger / 1e9;
            result.tv_nsec -= diff;
            remaining -= diff * 1e9;
        }

        if(remaining < result.tv_nsec) {
            result.tv_nsec -= remaining;
            remaining = 0;
        } else {
            int64_t nsec = 1e9 + result.tv_nsec - remaining;
            result.tv_sec -= 1;
            result.tv_nsec = nsec;
            remaining = 0;
        }
    }

    if(ns_per_cycle - size_t(ns_per_cycle) > 0) {
        int64_t remaining = delta_trigger * (ns_per_cycle - size_t(ns_per_cycle));

        if(remaining > 1e9) {
            int64_t diff = delta_trigger / 1e9;
            result.tv_nsec -= diff;
            remaining -= diff * 1e9;
        }

        if(remaining < result.tv_nsec) {
            result.tv_nsec -= remaining;
            remaining = 0;
        } else {
            int64_t nsec = 1e9 + result.tv_nsec - remaining;
            result.tv_sec -= 1;
            result.tv_nsec = nsec;
            remaining = 0;
        }
    }

    return result;
}

struct timespec subtract_ns_from_timespec(struct timespec spec, int64_t ns) {
    if(ns < 0)
        return spec;

    struct timespec result = spec;
    int64_t remaining = ns;

    if(remaining > 1e9) {
        int64_t diff = ns / 1e9;
        result.tv_nsec -= diff;
        remaining -= diff * 1e9;
    }

    if(remaining < result.tv_nsec) {
        result.tv_nsec -= remaining;
        remaining = 0;
    } else {
        int64_t nsec = 1e9 + result.tv_nsec - remaining;
        result.tv_sec -= 1;
        result.tv_nsec = nsec;
        remaining = 0;
    }

    return result;
}

struct timespec add_ns_to_timespec(struct timespec spec, int64_t ns) {
    if(ns < 0)
        return spec;

    struct timespec result = spec;
    int64_t remaining = ns;

    if(remaining > 1e9) {
        int64_t diff = ns / 1e9;
        result.tv_sec += diff;
        remaining -= diff * 1e9;
    }

    remaining += result.tv_nsec;
    if(remaining > 1e9) {
        int64_t diff = remaining / 1e9;
        result.tv_sec += diff;
        remaining -= diff * 1e9;
    }
    result.tv_nsec = remaining;

    return result;
}

struct timespec min_timespecs(std::vector<struct timespec> specs) {
    size_t N = specs.size();

    int64_t min_tv_sec = specs.at(0).tv_sec;
    int64_t min_tv_nsec = specs.at(0).tv_nsec;

    for(size_t i=1; i<N; ++i) {
        if((specs.at(i).tv_sec < min_tv_sec) or ((specs.at(i).tv_sec == min_tv_sec) and (specs.at(i).tv_nsec < min_tv_nsec))) {
            min_tv_sec = specs.at(i).tv_sec;
            min_tv_nsec = specs.at(i).tv_nsec;
        }
    }

    struct timespec result;
    result.tv_sec = min_tv_sec;
    result.tv_nsec = min_tv_nsec;

    return result;
}

struct timespec average_timespecs(std::vector<struct timespec> specs) {
    size_t N = specs.size();

    int64_t min_tv_sec = specs.at(0).tv_sec;

    for(size_t i=1; i<N; ++i) {
        min_tv_sec = std::min(min_tv_sec, int64_t(specs.at(i).tv_sec));
    }

    int64_t total_tv_sec = 0;
    int64_t total_tv_nsec = 0;
    for(size_t i=0; i<N; ++i) {
        total_tv_sec += specs.at(i).tv_sec - min_tv_sec;
        total_tv_nsec += specs.at(i).tv_nsec;
        if(total_tv_nsec > 1e9) {
            int64_t diff = total_tv_nsec / 1e9;
            total_tv_sec += diff;
            total_tv_nsec -= diff * 1e9;
        }
    }

    long double total_tv_sec_int = int64_t(total_tv_sec / N);
    long double total_tv_sec_fraction = ((long double)(total_tv_sec) / (long double)(N)) - total_tv_sec_int;
    total_tv_sec = total_tv_sec_int;
    total_tv_nsec /= N;
    total_tv_nsec += total_tv_sec_fraction * 1e9;
    if(total_tv_nsec > 1e9) {
        int64_t diff = total_tv_nsec / 1e9;
        total_tv_sec += diff;
        total_tv_nsec -= diff * 1e9;
    }

    total_tv_sec += min_tv_sec;

    struct timespec result;
    result.tv_sec = total_tv_sec;
    result.tv_nsec = total_tv_nsec;
    return result;
}

int64_t average_timediffs_in_ns(std::vector<int64_t> times) {
    constexpr long double ns_per_cycle = 8;

    size_t N = times.size();

    int64_t min_time = times.at(0);

    for(size_t i=1; i<N; ++i) {
        min_time = std::min(min_time, times.at(i));
    }

    int64_t total_time = 0;
    for(size_t i=0; i<N; ++i) {
        total_time += times.at(i) - min_time;
    }

    total_time /= N;
    total_time += min_time;

    return total_time * ns_per_cycle;
}

int64_t compute_trigger_offset(TimeReader & reader0, size_t board_idx0, TimeReader & reader1, size_t board_idx1, std::vector<int64_t> jitter_tests={-2, 2}, int64_t max_delta=2, size_t min_triggers=500, size_t max_triggers=2000, size_t increment=25, double threshold=0.9) {
    // Compute the time offset by testing various offsets and selecting the offset that results the the most valid time pairs

    // Number of triggers to test
    uint64_t n_triggers = min_triggers;

    // Boundaries of previously tested offsets
    int64_t prev_min = 0;
    int64_t prev_max = 0;

    // Boundaries of offsets to test
    int64_t new_min = -(n_triggers/2);
    int64_t new_max = (n_triggers/2) + 1;

    // Times used for testing the offsets
    std::vector<int64_t> times0 = reader0.GetTimes(n_triggers, board_idx0);
    std::vector<int64_t> times1 = reader1.GetTimes(n_triggers, board_idx1);

    // Generate offsets to test
    size_t x = prev_max;
    std::vector<int64_t> deltas;
    std::generate_n(std::back_inserter(deltas), new_max - prev_max, [&](){return x++;});
    x = new_min;
    std::generate_n(std::back_inserter(deltas), prev_min - new_min, [&](){return x++;});
    std::vector<std::tuple<int64_t, int64_t, size_t, size_t, size_t>> results;
    std::tuple<int64_t, int64_t, size_t, size_t, size_t> best_pair_result;

    while(true) {
        // Generate results for each offset
        for(int64_t const & delta_trigger : deltas) {
            int64_t delta = get_time_delta(times0, times1, delta_trigger);
            std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> pair_result = find_pairs(times0, times1, delta, max_delta);
            results.emplace_back(delta_trigger, delta, std::get<0>(pair_result).size(), std::get<1>(pair_result).size(), std::get<2>(pair_result).size());
            for(int64_t t : jitter_tests) {
                std::tuple<std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>, std::vector<std::tuple<IDX, IDX>>> j_pair_result = find_pairs(times0, times1, delta + t, max_delta);
                if(std::get<2>(j_pair_result).size() > std::get<2>(pair_result).size())
                    results.emplace_back(delta_trigger, delta + t, std::get<0>(j_pair_result).size(), std::get<1>(j_pair_result).size(), std::get<2>(j_pair_result).size());
            }
        }

        // Find the result with the most pairs
        best_pair_result = *std::max_element(results.begin(), results.end(), [](std::tuple<int64_t, int64_t, size_t, size_t, size_t> const & a, std::tuple<int64_t, int64_t, size_t, size_t, size_t> const & b){return std::get<3>(a) < std::get<3>(b);});

        int64_t delta_trigger = std::get<0>(best_pair_result);
        int64_t delta = std::get<1>(best_pair_result);
        size_t pairs = std::get<2>(best_pair_result);
        size_t good_pairs = std::get<3>(best_pair_result);
        size_t orphans = std::get<4>(best_pair_result);

        // Check if the fraction of good pairs exceeds the threshold
        if(double(good_pairs) / double(n_triggers - std::abs(delta_trigger)) >= threshold)
            break;

        // If the threshold is not met, increment the boundaries of the offsets, generate new offsets, and read the newly needed times

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

std::tuple<std::vector<std::vector<int64_t>>, struct timespec> compute_offsets(std::vector<std::vector<std::string>> file_lists, int64_t max_delta=2, std::vector<int64_t> jitter_tests={-2, 2}, size_t min_triggers=500, size_t max_triggers=2000, size_t increment=25, double threshold=0.9) {
    // Compute offsets for all boards
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
        offsets[i][0] = compute_trigger_offset(readers[0], 0, readers[i], 0, jitter_tests, max_delta, min_triggers, max_triggers, increment, threshold);
    for(size_t i=0; i<n_daqs; ++i)
        for(size_t j=1; j<n_boards[i]; ++j)
            offsets[i][j] = compute_trigger_offset(readers[0], 0, readers[i], j, jitter_tests, max_delta, min_triggers, max_triggers, increment, threshold);

    std::vector<std::vector<std::deque<int64_t>>> cached_times;
    std::vector<std::vector<std::deque<struct timespec>>> cached_computer_times;
    for(size_t i=0; i<n_daqs; ++i) {
        cached_times.push_back(readers[i].GetAllTimes());
        cached_computer_times.push_back(readers[i].GetAllComputerTimes());
    }

    bool have_computer_times;
    if(cached_computer_times.size() > 0) {
        have_computer_times = true;
        for(size_t i=0; i<n_daqs; ++i) {
            bool daq_has_times = true;
            for(size_t j=0; j<cached_times.at(i).size(); ++j) {
                if(cached_times.at(i).at(j).size() == 0) {
                    daq_has_times = false;
                }
            }
            if(not daq_has_times) {
                have_computer_times = false;
                break;
            }
        }
    } else {
        have_computer_times = false;
    }

    if(have_computer_times) {
        size_t min_size = cached_times.at(0).at(0).size();
        for(size_t i=0; i<n_daqs; ++i) {
            for(size_t j=0; j<cached_times.at(i).size(); ++j) {
                min_size = std::min(min_size, cached_times.at(i).at(j).size());
                min_size = std::min(min_size, cached_computer_times.at(i).at(j).size());
            }
        }

        for(size_t i=0; i<n_daqs; ++i) {
            for(size_t j=0; j<cached_times.at(i).size(); ++j) {
                cached_times.at(i).at(j).resize(min_size);
                cached_computer_times.at(i).at(j).resize(min_size);
            }
        }

        std::vector<struct timespec> timespec_samples;
        for(size_t i=0; i<n_daqs; ++i) {
            for(size_t j=0; j<cached_times.at(i).size(); ++j) {
                struct timespec sample = cached_computer_times.at(i).at(j).back();
                if(sample.tv_sec > 0 or sample.tv_nsec > 0) {
                    timespec_samples.push_back(sample);
                }
            }
        }

        struct timespec reference_time = min_timespecs(timespec_samples);
        for(size_t i=min_size-1; i>0; --i) {
            size_t upper = i;
            size_t lower = i-1;
            std::vector<int64_t> time_diffs;
            timespec_samples.clear();
            for(size_t j=0; j<n_daqs; ++j) {
                for(size_t k=0; k<cached_times.at(j).size(); ++k) {
                    int64_t time_diff = cached_times.at(j).at(k).at(upper) - cached_times.at(j).at(k).at(lower);
                    if(time_diff > 0)
                        time_diffs.push_back(time_diff);

                    struct timespec sample = cached_computer_times.at(j).at(k).at(lower);
                    if(sample.tv_sec > 0 or sample.tv_nsec > 0) {
                        timespec_samples.push_back(sample);
                    }
                }
            }
            int64_t ns = average_timediffs_in_ns(time_diffs);
            reference_time = subtract_ns_from_timespec(reference_time, ns);

            timespec_samples.push_back(reference_time);
            reference_time = min_timespecs(timespec_samples);
        }

        return {offsets, reference_time};
    } else {
        struct timespec reference_time;
        reference_time.tv_sec = 0;
        reference_time.tv_nsec = 0;

        return {offsets, reference_time};
    }
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
    struct timespec start_time;
    int run_number;

    bool parse_successful = false;
    int parsed_run_number = 0;
    int parsed_file_number = 0;
    struct timespec parsed_time;

    size_t counter = 0;
    size_t incomplete_counter = 0;

    // Important constants
    constexpr static long double ns_per_cycle = 8;

    bool PopFrame(size_t daq_idx);
    I3FramePtr GetConfigFrame();
    bool NextTrigger(size_t daq_idx, size_t board_idx);
    bool NextTriggers();
    void ClearUnusedFrames();
    std::tuple<boost::shared_ptr<I3Vector<I3Vector<uint16_t>>>, boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>>, boost::shared_ptr<I3Vector<std::pair<bool, int64_t>>>, boost::shared_ptr<CCMEventHeader>> GetTriggerReadout();
public:
    MergedSource(const I3Context&);
    void Configure();
    void Process();
    void Finish();

    SET_LOGGER("MergedSource");
};

I3_MODULE(MergedSource);

bool MergedSource::PopFrame(size_t daq_idx) {
    I3FramePtr frame;
    // Grab frames until we get a DAQ frame
    while(true) {
        // Fail if we do not have any more frames
        if(not frame_sequences[daq_idx]->more())
            return false;
        frame = frame_sequences[daq_idx]->pop_frame();
        if(frame->GetStop() == I3Frame::Geometry) {
            // Update stored configurations if we see a Geometry frame
            CCMAnalysis::Binary::CCMDAQConfigConstPtr new_config = frame->Get<CCMAnalysis::Binary::CCMDAQConfigConstPtr>("CCMDAQConfig");
            if(configs[daq_idx] == nullptr or *(configs[daq_idx]) != *new_config) {
                configs[daq_idx] = new_config;
                push_config = true;
                // Set the push_config flag to false if any configs are missing
                for(CCMAnalysis::Binary::CCMDAQConfigConstPtr c : configs)
                    if(c == nullptr)
                        push_config = false;
            }
            continue;
        } else if(frame->GetStop() != I3Frame::DAQ) {
            // SKip non-DAQ and non-Geometry frames
            continue;
        }
        // Skip this frame if we do not have the necessary information in the frame
        if((not frame->Has("CCMTriggerReadout")) or (frame->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout")->triggers.size() == 0))
            continue;
        break;
    }

    // Store the frame
    frame_cache[daq_idx].push_back(frame);
    CCMAnalysis::Binary::CCMTriggerReadoutConstPtr readout = frame->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout");
    std::vector<uint32_t> const & time_read = readout->triggers[0].board_times;

    // If the fill_computer_time flag is negative unset (i.e. negative), then set it
    if(fill_computer_time < 0)
        // Check if the computer times are stored in the input
        fill_computer_time = readout->triggers[0].board_computer_times.size() > 0;

    // Compute which boards have the appropriate data
    std::vector<uint8_t> mask = empty_mask(frame);
    for(size_t board_idx=0; board_idx<n_boards[daq_idx]; ++board_idx) {
        mask_cache[daq_idx][board_idx].push_back(mask[board_idx]);
        // Store inf in time cache if data is missing
        if(not mask[board_idx]) {
            time_cache[daq_idx][board_idx].push_back(std::numeric_limits<long double>::infinity());
            continue;
        }
        // Add times to cache
        uint32_t raw_time = time_read[board_idx];
        int64_t rel_time = CCMAnalysis::Binary::subtract_times(raw_time, last_raw_time[daq_idx][board_idx]);
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
    // Create a Geometry frame that contains the CCMDAQConfig composed of cached configurations from each DAQ, and the board time offsets
    CCMAnalysis::Binary::CCMDAQConfigPtr config = boost::make_shared<CCMAnalysis::Binary::CCMDAQConfig>();
    for(CCMAnalysis::Binary::CCMDAQConfigConstPtr c : configs) {
        std::copy(c->machine_configurations.begin(), c->machine_configurations.end(), std::back_inserter(config->machine_configurations));
        std::copy(c->digitizer_boards.begin(), c->digitizer_boards.end(), std::back_inserter(config->digitizer_boards));
    }
    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::Geometry);
    frame->Put("CCMDAQConfig", config, I3Frame::Geometry);
    boost::shared_ptr<I3Vector<I3Vector<int64_t>>> frame_offsets = boost::make_shared<I3Vector<I3Vector<int64_t>>>();
    frame_offsets->reserve(offsets.size());
    for(size_t i=0; i<offsets.size(); ++i)
        frame_offsets->emplace_back(offsets[i]);
    frame->Put("BoardTimeOffsets", frame_offsets, I3Frame::Geometry);
    return frame;
}

bool MergedSource::NextTrigger(size_t daq_idx, size_t board_idx) {
    // Current frame index
    size_t current_frame_idx = frame_idxs[daq_idx][board_idx];

    // Desired frame index
    size_t frame_idx = current_frame_idx + 1;

    // Add frames to the cache until the desired frame index is present
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

    // Set the current frame index to be the desired frame index
    frame_idxs[daq_idx][board_idx] = frame_idx;
    return true;
}

bool MergedSource::NextTriggers() {
    // Call NextTrigger for each board
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
    // Append trigger information from one DAQ into a combined ouput
    CCMAnalysis::Binary::CCMTrigger & o = (*output_triggers)[0];
    CCMAnalysis::Binary::CCMTrigger const & s = source_trigger->triggers[0];
    assert(s.channel_sizes.size() > first_idx and s.channel_sizes.size() >= last_idx);
    assert(s.channel_masks.size() > first_idx and s.channel_masks.size() >= last_idx);
    assert(s.channel_temperatures.size() > first_idx and s.channel_sizes.size() >= last_idx);
    assert(s.board_event_numbers.size() > board_idx);
    assert(s.board_times.size() > board_idx);

    // Copy channel_sizes, channel_masks, channel_temperatures, board_event_times, board_event_numbers, and board_times
    std::copy(s.channel_sizes.begin() + first_idx, s.channel_sizes.begin() + last_idx, std::back_inserter(o.channel_sizes));
    std::copy(s.channel_masks.begin() + first_idx, s.channel_masks.begin() + last_idx, std::back_inserter(o.channel_masks));
    std::copy(s.channel_temperatures.begin() + first_idx, s.channel_temperatures.begin() + last_idx, std::back_inserter(o.channel_temperatures));
    std::copy(s.board_event_numbers.begin() + board_idx, s.board_event_numbers.begin() + board_idx + 1, std::back_inserter(o.board_event_numbers));
    std::copy(s.board_times.begin() + board_idx, s.board_times.begin() + board_idx + 1, std::back_inserter(o.board_times));

    // Optionall copy board_computer_times
    if(fill_computer_time) {
        assert(s.board_computer_times.size() > board_idx);
        std::copy(s.board_computer_times.begin() + board_idx, s.board_computer_times.begin() + board_idx + 1, std::back_inserter(o.board_computer_times));
    }

    // Check if all expected samples are present
    if(source_trigger->samples.size() == s.channel_masks.size()) {
        // Copy all the samples in the specified range
        assert(source_trigger->samples.size() > first_idx and source_trigger->samples.size() >= last_idx);
        size_t n_channels = last_idx - first_idx;
        output_samples->reserve(output_samples->size() + n_channels);
        for(size_t i=first_idx; i<last_idx; ++i) {
            output_samples->emplace_back(source_trigger->samples[i]);
        }
    } else {
        // Iterate over channel_mask
        // Only increment counter if channel_mask is set to true
        // Counter represents the number of samples actually present in the samples vector
        size_t new_first_idx = 0;
        for(size_t i=0; i<first_idx; ++i) {
            if(s.channel_masks[i])
                ++new_first_idx;
        }
        // Assume all samples we want are present, so the last_idx is just the first_idx plus the difference
        size_t new_last_idx = new_first_idx + last_idx - first_idx;
        assert(source_trigger->samples.size() > new_first_idx and source_trigger->samples.size() >= new_last_idx);
        // Copy all the samples in the specified range
        size_t n_channels = new_last_idx - new_first_idx;
        output_samples->reserve(output_samples->size() + n_channels);
        for(size_t i=new_first_idx; i<new_last_idx; ++i) {
            output_samples->emplace_back(source_trigger->samples[i]);
        }
    }
}

inline void merge_empty_trigger(
        boost::shared_ptr<I3Vector<I3Vector<uint16_t>>> output_samples,
        boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>> output_triggers,
        size_t n_channels,
        bool fill_computer_time = false) {
    // Append an empty trigger to the output
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

std::tuple<boost::shared_ptr<I3Vector<I3Vector<uint16_t>>>, boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>>, boost::shared_ptr<I3Vector<std::pair<bool, int64_t>>>, boost::shared_ptr<CCMEventHeader>>
    MergedSource::GetTriggerReadout() {
    CCMAnalysis::Binary::CCMTriggerReadoutPtr readout = boost::make_shared<CCMAnalysis::Binary::CCMTriggerReadout>();
    boost::shared_ptr<I3Vector<I3Vector<uint16_t>>> output_samples = boost::make_shared<I3Vector<I3Vector<uint16_t>>>();
    boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>> output_triggers = boost::make_shared<I3Vector<CCMAnalysis::Binary::CCMTrigger>>(1);
    std::vector<std::vector<long double>> times;
    times.resize(n_daqs);
    // Fill the times from current frame indices
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

    // Check if we have run out of triggers
    if(all_bad or std::isinf(min_time))
        return {nullptr, nullptr, nullptr, nullptr};

    std::vector<std::vector<int>> state(n_daqs);
    bool is_incomplete = false;
    boost::shared_ptr<I3Vector<std::pair<bool, int64_t>>> output_times(new I3Vector<std::pair<bool, int64_t>>());
    std::vector<int64_t> sample_times;
    for(size_t daq_idx=0; daq_idx < n_daqs; ++daq_idx) {
        size_t last_idx = 0;
        for(size_t board_idx=0; board_idx < n_boards[daq_idx]; ++board_idx) {
            size_t n_channels = configs[daq_idx]->digitizer_boards[board_idx].channels.size();
            // channel_sizes.push_back(n_channels);
            size_t next_idx = last_idx + n_channels;
            size_t frame_idx = frame_idxs[daq_idx][board_idx];
            if(empty[daq_idx][board_idx]) {
                is_incomplete = true;
                merge_empty_trigger(output_samples, output_triggers, n_channels, fill_computer_time);
                output_times->emplace_back(false, 0);
                state[daq_idx].push_back(0);
            } else if(std::isinf(times[daq_idx][board_idx])) {
                is_incomplete = true;
                state[daq_idx].push_back(1);
                log_fatal("Should not see inf here!");
                merge_empty_trigger(output_samples, output_triggers, n_channels, fill_computer_time);
                output_times->emplace_back(false, 0);
            } else if(times[daq_idx][board_idx] - min_time > max_delta) {
                is_incomplete = true;
                merge_empty_trigger(output_samples, output_triggers, n_channels, fill_computer_time);
                output_times->emplace_back(false, 0);
                state[daq_idx].push_back(2);
            } else {
                // Fill the trigger output from the appropriate input
                CCMAnalysis::Binary::CCMTriggerReadoutConstPtr tr = frame_cache[daq_idx][frame_idx]->Get<CCMAnalysis::Binary::CCMTriggerReadoutConstPtr>("CCMTriggerReadout");
                merge_triggers(output_samples, output_triggers, tr, board_idx, last_idx, next_idx, fill_computer_time);
                output_times->emplace_back(true, times[daq_idx][board_idx]);
                sample_times.push_back(times[daq_idx][board_idx]);
                // Grab the next trigger
                bool res = NextTrigger(daq_idx, board_idx);
                if(not res) {
                    times[daq_idx][board_idx] = std::numeric_limits<long double>::infinity();
                    empty[daq_idx][board_idx] = true;
                    state[daq_idx].push_back(4);
                } else {
                    state[daq_idx].push_back(3);
                }
            }
            last_idx = next_idx;
        }
    }
    size_t max_wf_size = *std::max_element(output_triggers->at(0).channel_sizes.begin(), output_triggers->at(0).channel_sizes.end());

    int64_t avg_time_in_ns = average_timediffs_in_ns(sample_times);
    struct timespec current_time = add_ns_to_timespec(start_time, avg_time_in_ns);

    boost::shared_ptr<CCMEventHeader> header = boost::make_shared<CCMEventHeader>();
    CCMTime event_start_time(current_time.tv_sec, current_time.tv_nsec);
    CCMTime event_end_time = event_start_time + max_wf_size * 2;

    header->SetRunID(run_number);
    header->SetSubRunID(0);
    header->SetEventID(counter);
    header->SetSubEventID(0);
    header->SetStartTime(event_start_time);
    header->SetEndTime(event_end_time);

    if(is_incomplete)
        ++incomplete_counter;
    // Clear out frames that can no longer be referenced
    ClearUnusedFrames();
    return {output_samples, output_triggers, output_times, header};
}

MergedSource::MergedSource(const I3Context& context) : I3Module(context) {
    AddParameter("FileLists",
            "File lists to merge",
            std::vector<std::vector<std::string>>());

    AddParameter("MaxTimeDiff",
            "Maximum time difference between associated triggers (ns)",
            16);
    AddParameter("RunNumber",
            "Run numbner",
            int(-1));
}

void MergedSource::Configure() {
    std::vector<std::vector<std::string>> file_lists;
    GetParameter("FileLists", file_lists);

    double max_time_diff;
    GetParameter("MaxTimeDiff", max_time_diff);
    max_delta = std::ceil(max_time_diff / ns_per_cycle);

    GetParameter("RunNumber", run_number);

    parse_successful = false;
    parsed_run_number = 0;
    parsed_file_number = 0;
    parsed_time.tv_sec = 0;
    parsed_time.tv_nsec = 0;

    for(size_t i=0; i<file_lists.size(); ++i) {
        std::string fname = file_lists.at(i).front();
        int this_parsed_run_number;
        int this_parsed_file_number;
        struct timespec this_parsed_time;
        bool okay = ParseFileName(fname, this_parsed_run_number, this_parsed_file_number, this_parsed_time);
        if(not okay)
            continue;
        if(parse_successful) {
            parsed_file_number = std::min(parsed_file_number, this_parsed_file_number);
            if(parsed_time.tv_sec > this_parsed_time.tv_sec or ((parsed_time.tv_sec == this_parsed_time.tv_sec) and (parsed_time.tv_nsec > this_parsed_time.tv_nsec))) {
                parsed_time = this_parsed_time;
            }
        } else {
            parsed_run_number = this_parsed_run_number;
            parsed_file_number = this_parsed_file_number;
            parsed_time = this_parsed_time;
        }
        parse_successful = true;
    }

    // Compute the offsets to use for trigger alignment
    std::tuple<std::vector<std::vector<int64_t>>, struct timespec> offset_results = compute_offsets(file_lists, max_delta);
    offsets = std::get<0>(offset_results);
    start_time = std::get<1>(offset_results);

    if(start_time.tv_sec == 0 and start_time.tv_nsec == 0) {
        if(parse_successful) {
            log_warn("Start time cannot be determined from i3 file. Defaulting to file timestamp.");
            start_time = parsed_time;
        } else {
            log_warn("Start time cannot be determined from i3 file or file timestamp. Defaulting to zero");
        }
    }

    if(run_number == -1) {
        if(parse_successful) {
            run_number = parsed_run_number;
        } else {
            run_number = 0;
            log_warn("Run number is not set by Module Parameter or by parsed filename. Defaulting to zero.");
        }
    }

    int64_t max_offset = 0;
    int64_t min_offset = 0;
    for(size_t daq_idx=0; daq_idx < offsets.size(); ++daq_idx) {
        for(size_t board_idx=0; board_idx < offsets[daq_idx].size(); ++board_idx) {
            int64_t x = offsets[daq_idx][board_idx];
            if(x > max_offset)
                max_offset = x;
            if(x < min_offset)
                min_offset = x;
        }
    }
    std::stringstream ss;
    ss << max_offset;
    size_t width = ss.str().size();
    ss.str("");
    ss << min_offset;
    width = std::max(width, ss.str().size());
    ss.str("");
    ss << "Board offsets computed:";
    for(size_t daq_idx=0; daq_idx < offsets.size(); ++daq_idx) {
        ss << "\n[";
        bool first = true;
        for(size_t board_idx=0; board_idx < offsets[daq_idx].size(); ++board_idx) {
            if(first)
                first = false;
            else
                ss << ", ";
            ss << std::setw(width) << offsets[daq_idx][board_idx];
        }
        ss << "]";
    }
    ss << "\n";
    log_notice("%s", ss.str().c_str());

    // Set up all the data structures
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

    // Grab the first set of triggers
    NextTriggers();
}

void MergedSource::Process() {
    // Push the Geometry frame and config if needed
    if(push_config) {
        PushFrame(GetConfigFrame());
        push_config = false;
    }

    // Create a new frame
    I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::DAQ);

    // Get the readout output
    std::tuple<boost::shared_ptr<I3Vector<I3Vector<uint16_t>>>, boost::shared_ptr<I3Vector<CCMAnalysis::Binary::CCMTrigger>>, boost::shared_ptr<I3Vector<std::pair<bool, int64_t>>>, boost::shared_ptr<CCMEventHeader>> readout = GetTriggerReadout();

    // Exit if we cannot create any more output
    if(std::get<0>(readout) == nullptr or std::get<1>(readout) == nullptr or std::get<2>(readout) == nullptr or std::get<3>(readout) == nullptr) {
        RequestSuspension();
        return;
    }

    // Put the output objects in the frame and push the frame
    frame->Put("CCMDigitalReadout", std::get<0>(readout), I3Frame::DAQ);
    frame->Put("CCMTriggers", std::get<1>(readout), I3Frame::DAQ);
    frame->Put("TriggerTimes", std::get<2>(readout), I3Frame::DAQ);
    frame->Put("CCMEventHeader", std::get<3>(readout), I3Frame::DAQ);
    PushFrame(frame);
    ++counter;
}

void MergedSource::Finish() {
    log_notice_stream(
        "Merged " << offsets.size() << " DAQ streams into " << counter << " separate triggers. Encountered "
        << counter-incomplete_counter << " complete triggers and " << incomplete_counter << " incomplete triggers");
}
