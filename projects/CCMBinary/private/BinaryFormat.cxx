#ifndef CCMAnalysis_BinaryFormat_CXX
#define CCMAnalysis_BinaryFormat_CXX

#include <ctime>
#include <cstdio>
#include <string>
#include <time.h>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>

#include <icetray/I3Logging.h>
#include <icetray/serialization.h>

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

bool operator==(timespec const & a, timespec const & b) {
    return (a.tv_sec == b.tv_sec) and (a.tv_nsec == b.tv_nsec);
}

namespace CCMAnalysis {
namespace Binary {

template <class Archive>
void
ChannelHeader::serialize(Archive& ar, unsigned version) {
    if(version > channelheader_version_)
        log_fatal("Attempting to read version %u from file but running version %u of ChannelHeader class.", version, channelheader_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("physical_board_id", physical_board_id);
    ar & make_nvp("board_serial_number", board_serial_number);
    ar & make_nvp("physical_channel_type", physical_channel_type);
    ar & make_nvp("physical_channel_id", physical_channel_id);
    ar & make_nvp("caen_optical_link_number", caen_optical_link_number);
    ar & make_nvp("caen_optical_link_board_number", caen_optical_link_board_number);
    ar & make_nvp("caen_channel_number", caen_channel_number);
}

bool ChannelHeader::operator==(ChannelHeader const & other) const {
    return std::tie(
    physical_board_id,
    board_serial_number,
    physical_channel_type,
    physical_channel_id,
    caen_optical_link_number,
    caen_optical_link_board_number,
    caen_channel_number)
      == std::tie(
    other.physical_board_id,
    other.board_serial_number,
    other.physical_channel_type,
    other.physical_channel_id,
    other.caen_optical_link_number,
    other.caen_optical_link_board_number,
    other.caen_channel_number);
}

bool ChannelHeader::operator!=(ChannelHeader const & other) const {
    return not this->operator==(other);
}

std::ostream & ChannelHeader::Print(std::ostream & os) const {
    os << '[';
    os << physical_board_id;
    os << '(' << caen_optical_link_number << ',' << caen_optical_link_board_number << ',' << caen_channel_number << "), ";
    os << physical_channel_type << ", ";
    os << physical_channel_id << ']';
    return os;
}

template <class Archive>
void
DigitizerBoard::serialize(Archive& ar, unsigned version) {
    if(version > digitizerboard_version_)
        log_fatal("Attempting to read version %u from file but running version %u of DigitizerBoard class.", version, digitizerboard_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("physical_board_id", physical_board_id);
    ar & make_nvp("board_serial_number", board_serial_number);
    ar & make_nvp("caen_optical_link_number", caen_optical_link_number);
    ar & make_nvp("caen_optical_link_board_number", caen_optical_link_board_number);
    ar & make_nvp("channels", channels);
    ar & make_nvp("trigger_out_type", trigger_out_type);
}

bool DigitizerBoard::operator==(DigitizerBoard const & other) const {
    return std::tie(
    physical_board_id,
    board_serial_number,
    caen_optical_link_number,
    caen_optical_link_board_number,
    channels,
    trigger_out_type)
        == std::tie(
    other.physical_board_id,
    other.board_serial_number,
    other.caen_optical_link_number,
    other.caen_optical_link_board_number,
    other.channels,
    other.trigger_out_type);
}

bool DigitizerBoard::operator!=(DigitizerBoard const & other) const {
    return not this->operator==(other);
}

std::ostream & DigitizerBoard::Print(std::ostream & os) const {
    os << "[Board:";
    os << physical_board_id;
    os << '(' << caen_optical_link_number << ',' << caen_optical_link_board_number << "),\n";
    os << "TriggerOutType:" << trigger_out_type << ", \n";
    if(channels.size() > 0) {
        os << "Channels:\n";
        bool first = true;
        for(size_t i=0; i<channels.size(); ++i) {
            if(first)
                first = false;
            else {
                os << ',';
                os << '\n';
            }
            os << '[' << channels[i].physical_channel_type << ": ";
            os << channels[i].physical_channel_id << ']';
        }
    } else {
        os << "Channels: []";
    }
    os << "\n]";
    return os;
}

template <class Archive>
void
CCMDAQMachineConfig::serialize(Archive& ar, unsigned version) {
    if(version > ccmdaqmachineconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMDAQMachineConfig class.", version, ccmdaqmachineconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("machine_identifier", machine_identifier);
    ar & make_nvp("num_digitizer_boards", num_digitizer_boards);
    ar & make_nvp("num_channels", num_channels);
    ar & make_nvp("num_samples", num_samples);
    ar & make_nvp("trigger_percent_after", trigger_percent_after);
    ar & make_nvp("trigger_time_tolerance", trigger_time_tolerance);
    ar & make_nvp("missed_trigger_tolerance", missed_trigger_tolerance);
    ar & make_nvp("offset_estimate_min_triggers", offset_estimate_min_triggers);
    ar & make_nvp("offset_estimate_abs_error_threshold", offset_estimate_abs_error_threshold);
    ar & make_nvp("offset_estimate_rel_error_threshold", offset_estimate_rel_error_threshold);
    ar & make_nvp("offset_estimate_tau", offset_estimate_tau);
}

bool CCMDAQMachineConfig::operator==(CCMDAQMachineConfig const & other) const {
    return std::tie(
    machine_identifier,
    num_digitizer_boards,
    num_channels,
    num_samples,
    trigger_percent_after,
    trigger_time_tolerance,
    missed_trigger_tolerance,
    offset_estimate_min_triggers,
    offset_estimate_abs_error_threshold,
    offset_estimate_rel_error_threshold,
    offset_estimate_tau)
        == std::tie(
    other.machine_identifier,
    other.num_digitizer_boards,
    other.num_channels,
    other.num_samples,
    other.trigger_percent_after,
    other.trigger_time_tolerance,
    other.missed_trigger_tolerance,
    other.offset_estimate_min_triggers,
    other.offset_estimate_abs_error_threshold,
    other.offset_estimate_rel_error_threshold,
    other.offset_estimate_tau);
}

bool CCMDAQMachineConfig::operator!=(CCMDAQMachineConfig const & other) const {
    return not this->operator==(other);
}

std::ostream & CCMDAQMachineConfig::Print(std::ostream & os) const {
    os << "[Machine:" << machine_identifier << ",\n";
    os << "NBoards:" << num_digitizer_boards << ",\n";
    os << "NChannels:" << num_channels << ",\n";
    os << "NSamples:" << num_samples << ",\n";
    os << "TriggerPercentAfter:" << trigger_percent_after << ",\n";
    os << "TriggerTimeTolerance:" << trigger_time_tolerance << ",\n";
    os << "MissedTriggerTolerance:" << missed_trigger_tolerance << ",\n";
    os << "OffsetEstimateMinTriggers:" << offset_estimate_min_triggers << ",\n";
    os << "OffsetEstimateAbsErrorThreshold:" << offset_estimate_abs_error_threshold << ",\n";
    os << "OffsetEstimateRelErrorThreshold:" << offset_estimate_rel_error_threshold << ",\n";
    os << "OffsetEstimateTau:" << offset_estimate_tau << '\n';
    os << ']';
    return os;
}

template <class Archive>
void
CCMDAQConfig::serialize(Archive& ar, unsigned version) {
    if(version > ccmdaqconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMDAQConfig class.", version, ccmdaqconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("machine_configurations", machine_configurations);
    ar & make_nvp("digitizer_boards", digitizer_boards);
}

bool CCMDAQConfig::operator==(CCMDAQConfig const & other) const {
    return std::tie(machine_configurations, digitizer_boards)
        == std::tie(other.machine_configurations, other.digitizer_boards);
}

bool CCMDAQConfig::operator!=(CCMDAQConfig const & other) const {
    return not this->operator==(other);
}

std::ostream & CCMDAQConfig::Print(std::ostream & os) const {
    os << "[Machines:[\n";
    bool first = true;
    size_t board_idx = 0;
    for(size_t i=0; i<machine_configurations.size(); ++i) {
        if(first)
            first = false;
        else
            os << ",\n";
        os << machine_configurations[i];
        os << "\n[Boards:";
        size_t max_idx = board_idx + machine_configurations[i].num_digitizer_boards;
        bool bfirst = true;
        for(; board_idx<max_idx; ++board_idx) {
            if(bfirst)
                bfirst = false;
            else
                os << ", ";
            os << digitizer_boards[board_idx].physical_board_id;
        }
    }
    os << "],\n";
    os << "[BoardConfigurations:\n";
    first = true;
    for(size_t i=0; i<digitizer_boards.size(); ++i) {
        if(first)
            first = false;
        else
            os << ",\n";
        os << digitizer_boards[i];
    }
    os << "]\n]";
    return os;
}

template <class Archive>
void
CCMTrigger::serialize(Archive& ar, unsigned version) {
    if(version > ccmtrigger_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMTrigger class.", version, ccmtrigger_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("channel_sizes", channel_sizes);
    ar & make_nvp("channel_masks", channel_masks);
    ar & make_nvp("channel_temperatures", channel_temperatures);
    ar & make_nvp("board_event_numbers", board_event_numbers);
    ar & make_nvp("board_times", board_times);
    ar & make_nvp("board_computer_times", board_computer_times);
}

namespace {
std::string GetTimeStringToSecond(time_t now) {
    char buf[sizeof "2011-10-08T07:07:09Z"];
    std::strftime(buf, sizeof buf, "%FT%TZ", gmtime(&now));
    return std::string(buf);
}
}

bool CCMTrigger::operator==(CCMTrigger const & other) const {
    return std::tie(
    channel_sizes,
    channel_masks,
    channel_temperatures,
    board_event_numbers,
    board_times,
    board_computer_times)
      == std::tie(
    other.channel_sizes,
    other.channel_masks,
    other.channel_temperatures,
    other.board_event_numbers,
    other.board_times,
    other.board_computer_times);
}

bool CCMTrigger::operator!=(CCMTrigger const & other) const {
    return not this->operator==(other);
}

std::ostream & CCMTrigger::Print(std::ostream & os) const {
    os << "[CCMTrigger\n";
    os << "ChannelSizes:[";
    bool first = true;
    for(size_t i=0; i< channel_sizes.size(); ++i) {
        if(first) first = false; else os << ",";
        os << channel_sizes[i];
    }
    os << "],\n";
    os << "ChannelMasks:[";
    first = true;
    for(size_t i=0; i< channel_masks.size(); ++i) {
        if(first) first = false; else os << ",";
        os << channel_masks[i];
    }
    os << "],\n";
    os << "ChannelTemperatures:[";
    first = true;
    for(size_t i=0; i< channel_temperatures.size(); ++i) {
        if(first) first = false; else os << ",";
        os << channel_temperatures[i];
    }
    os << "],\n";
    os << "BoardEventNumbers:[";
    first = true;
    for(size_t i=0; i<board_event_numbers.size(); ++i) {
        if(first) first = false; else os << ",";
        os << board_event_numbers[i];
    }
    os << "],\n";
    os << "BoardTimes:[";
    first = true;
    for(size_t i=0; i<board_times.size(); ++i) {
        if(first) first = false; else os << ",";
        os << "(" << (bool((0x1 << 31) & board_times[i])) << ")" << ((~(0x1 << 31)) & board_times[i]);
    }
    os << "],\n";
    first = true;
    os << "BoardComputerTimes:[";
    for(size_t i=0; i<board_computer_times.size(); ++i) {
        if(first) first = false; else os << ",";
        os << GetTimeStringToSecond(board_computer_times[i].tv_sec) << "." << std::setfill('0') << std::setw(9) << board_computer_times[i].tv_nsec % 1000000000l;
    }
    os << "]\n";
    os << "]";
    return os;
}


template <class Archive>
void
CCMTriggerReadout::serialize(Archive& ar, unsigned version) {
    if(version > ccmtriggerreadout_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMTriggerReadout class.", version, ccmtriggerreadout_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("event_number", event_number);
    ar & make_nvp("computer_time", computer_time);
    ar & make_nvp("triggers", triggers);
    ar & make_nvp("samples", samples);
}

bool CCMTriggerReadout::operator==(CCMTriggerReadout const & other) const {
    return std::tie(
        this->event_number,
        this->computer_time,
        this->triggers,
        this->samples) == std::tie(
        other.event_number,
        other.computer_time,
        other.triggers,
        other.samples);
}

bool CCMTriggerReadout::operator!=(CCMTriggerReadout const & other) const {
    return not this->operator==(other);
}

std::ostream & CCMTriggerReadout::Print(std::ostream & os) const {
    os << "[CCMTriggerReadout\n";
    os << "ComputerTime:";
    os << GetTimeStringToSecond(computer_time.tv_sec) << "." << std::setfill('0') << std::setw(9) << computer_time.tv_nsec % 1000000000l;
    os << "\n";
    os << "Triggers:[";
    bool first = true;
    for(size_t i=0; i<triggers.size(); ++i) {
        if(first)
            first = false;
        else
            os << ",";
        os << "\n";
        os << triggers[i];
    }
    os << "],\n";
    os << "Samples:[";
    first = true;
    for(size_t i=0; i< samples.size(); ++i) {
        if(first)
            first = false;
        else
            os << ",";
        os << samples[i].size();
    }
    os << "]\n]";
    return os;
}

template <class Archive>
void
CCMData::serialize(Archive& ar, unsigned version) {
    if(version > ccmdata_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMData class.", version, ccmdata_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("daq_config", daq_config);
    ar & make_nvp("trigger_readout", trigger_readout);
}

bool CCMData::operator==(CCMData const & other) const {
    return std::tie(
        this->daq_config,
        this->trigger_readout) == std::tie(
        other.daq_config,
        other.trigger_readout);
}

bool CCMData::operator!=(CCMData const & other) const {
    return not this->operator==(other);
}

} // namespace Binary
} // namespace CCMAnalsysis

std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::ChannelHeader c) {
    return(c.Print(os));
}
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::DigitizerBoard c) {
    return(c.Print(os));
}
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMDAQMachineConfig c) {
    return(c.Print(os));
}
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMDAQConfig c) {
    return(c.Print(os));
}
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMTrigger c) {
    return(c.Print(os));
}
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMTriggerReadout c) {
    return(c.Print(os));
}

I3_SERIALIZABLE(CCMAnalysis::Binary::ChannelHeader);
I3_SERIALIZABLE(CCMAnalysis::Binary::DigitizerBoard);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMDAQMachineConfig);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMDAQConfig);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMTrigger);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMTriggerReadout);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMData);

I3_SERIALIZABLE(CCMAnalysis::Binary::I3VectorCCMTrigger);

#endif // CCMAnalysis_BinaryFormat_CXX
