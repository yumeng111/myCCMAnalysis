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
    std::vector<std::vector<int64_t>> time_cache;
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
    TimeReader(std::vector<std::string> const & file_names, size_t n_skip=0) :
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

    }

    std::vector<int64_t> GetTimes(size_t N, size_t board_idx) {
        while(time_cache[board_idx].size() < N) {
            PopTimes();
        }
        std::vector<int64_t> result(N);
        std::copy(std::begin(time_cache[board_idx]), std::begin(time_cache[board_idx]) + N, std::begin(result));
        return result;
    }
};



