#ifndef CCMAnalysis_BinaryUtilities_H
#define CCMAnalysis_BinaryUtilities_H

#include <string>
#include <vector>
#include <iostream>

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

namespace CCMAnalysis {
namespace Binary {

inline std::ostream & write_binary(std::ostream & os, uint32_t const & x);
inline std::istream & read_binary(std::istream & is, uint32_t & s);
inline std::ostream & write_binary(std::ostream & os, uint16_t const & x);
inline std::istream & read_binary(std::istream & is, uint16_t & s);
inline std::ostream & write_binary(std::ostream & os, float const & x);
inline std::istream & read_binary(std::istream & is, float & s);
inline std::ostream & write_binary(std::ostream & os, double const & x);
inline std::istream & read_binary(std::istream & is, double & s);
inline std::ostream & write_binary(std::ostream & os, long double const & x);
inline std::istream & read_binary(std::istream & is, long double & s);
inline std::ostream & write_binary(std::ostream & os, char const * x, size_t n);
inline std::istream & read_binary(std::istream & is, char * x, size_t n);
inline std::ostream & write_binary(std::ostream & os, struct timespec const & x);
inline std::istream & read_binary(std::istream & is, struct timespec & s);
inline std::ostream & write_binary(std::ostream & os, std::string const & s);
inline std::istream & read_binary(std::istream & is, std::string & s);
inline std::ostream & write_binary(std::ostream & os, ChannelHeader const & header);
inline std::istream & read_binary(std::istream & is, ChannelHeader & header);
inline std::ostream & write_binary(std::ostream & os, DigitizerBoard const & board);
inline std::istream & read_binary(std::istream & is, DigitizerBoard & board);
inline std::ostream & write_binary(std::ostream & os, CCMDAQMachineConfig const & config);
inline std::istream & read_binary(std::istream & is, CCMDAQMachineConfig & config);
inline std::ostream & write_binary(std::ostream & os, CCMDAQConfig const & config);
inline std::istream & read_binary(std::istream & is, CCMDAQConfig & config);
inline std::ostream & write_binary(std::ostream & os, CCMTrigger const & trigger);
inline std::istream & read_binary(std::istream & is, CCMTrigger const & trigger);
inline std::ostream & write_binary(std::ostream & os, CCMTriggerReadout const & trigger_readout);
inline std::istream & read_binary(std::istream & is, CCMTriggerReadout & trigger_readout);
inline std::ostream & write_magic_number(std::ostream & os);
inline std::istream & read_magic_number(std::istream & is, bool & result);
inline std::ostream & write_binary(std::ostream & os, CCMData const & data);
inline std::istream & read_config(std::istream & is, CCMData & data);
inline std::istream & read_size(std::istream & is, size_t & size);
inline std::istream & read_binary(std::istream & is, CCMData & data);

//bool IsVersionedBinary(std::)

template<typename T>
inline std::ostream & write_binary(std::ostream & os, std::vector<T> const & v);

template<typename T>
inline std::istream & read_binary(std::istream & is, std::vector<T> & v);

template<typename T>
inline std::ostream & write_binary_contiguous_vector(std::ostream & os, std::vector<T> const & v);

template<typename T>
inline std::istream & read_binary_contiguous_vector(std::istream & is, std::vector<T> & v);

template<>
inline std::ostream & write_binary(std::ostream & os, std::vector<uint32_t> const & v);
template<>
inline std::istream & read_binary(std::istream & is, std::vector<uint32_t> & v);

template<>
inline std::ostream & write_binary(std::ostream & os, std::vector<uint16_t> const & v);
template<>
inline std::istream & read_binary(std::istream & is, std::vector<uint16_t> & v);

template<>
inline std::ostream & write_binary(std::ostream & os, std::vector<float> const & v);
template<>
inline std::istream & read_binary(std::istream & is, std::vector<float> & v);

inline void merge_triggers(
        CCMAnalysis::Binary::CCMTriggerReadout & output_trigger,
        CCMAnalysis::Binary::CCMTriggerReadout && source_trigger,
        size_t board_idx,
        size_t first_idx,
        size_t last_idx,
        bool fill_computer_time = false) {
    CCMAnalysis::Binary::CCMTrigger & o = output_trigger.triggers[0];
    CCMAnalysis::Binary::CCMTrigger & s = source_trigger.triggers[0];
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
    if(source_trigger.samples.size() == s.channel_masks.size()) {
        assert(source_trigger.samples.size() > first_idx and source_trigger.samples.size() >= last_idx);
        std::copy(source_trigger.samples.begin() + first_idx, source_trigger.samples.begin() + last_idx, std::back_inserter(output_trigger.samples));
    } else {
        size_t new_first_idx = 0;
        for(size_t i=0; i<first_idx; ++i) {
            if(s.channel_masks[i])
                ++new_first_idx;
        }
        size_t new_last_idx = new_first_idx + last_idx - first_idx;
        assert(source_trigger.samples.size() > new_first_idx and source_trigger.samples.size() >= new_last_idx);
        std::copy(source_trigger.samples.begin() + new_first_idx, source_trigger.samples.begin() + new_last_idx, std::back_inserter(output_trigger.samples));
    }
}

inline void merge_triggers(
        CCMAnalysis::Binary::CCMTriggerReadout & output_trigger,
        CCMAnalysis::Binary::CCMTriggerReadout const & source_trigger,
        size_t board_idx,
        size_t first_idx,
        size_t last_idx,
        bool fill_computer_time = false) {
    CCMAnalysis::Binary::CCMTriggerReadout trigger = source_trigger;
    return merge_triggers(output_trigger, std::move(trigger), board_idx, first_idx, last_idx, fill_computer_time);
}

inline void merge_empty_trigger(
        CCMAnalysis::Binary::CCMTriggerReadout & output_trigger,
        size_t n_channels,
        bool fill_computer_time = false) {
    CCMAnalysis::Binary::CCMTrigger & o = output_trigger.triggers[0];
    std::fill_n(std::back_inserter(o.channel_sizes), n_channels, 0);
    std::fill_n(std::back_inserter(o.channel_masks), n_channels, 0);
    std::fill_n(std::back_inserter(o.channel_temperatures), n_channels, 0);
    output_trigger.samples.reserve(output_trigger.samples.size() + n_channels);
    for(size_t i=0; i<n_channels; ++i) {
        output_trigger.samples.emplace_back();
    }
    //std::fill_n(std::back_inserter(output_trigger.samples), n_channels, std::vector<uint16_t>());
    std::fill_n(std::back_inserter(o.board_event_numbers), 1, 0);
    std::fill_n(std::back_inserter(o.board_times), 1, 0);
    if(fill_computer_time) {
        struct timespec t; t.tv_sec = 0; t.tv_nsec = 0;
        std::fill_n(std::back_inserter(o.board_computer_times), 1, t);
    }
}

int32_t subtract_times(uint32_t t0, uint32_t t1);

inline bool operator <(timespec const & lhs, timespec const & rhs) {
    if (lhs.tv_sec == rhs.tv_sec)
        return lhs.tv_nsec < rhs.tv_nsec;
    else
        return lhs.tv_sec < rhs.tv_sec;
}



} // namespace Binary
} // namespace CCMAnalsysis

#include "CCMAnalysis/CCMBinary/BinaryUtilities.cxx"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.tcc"

#endif // CCMAnalysis_BinaryUtilities_H
