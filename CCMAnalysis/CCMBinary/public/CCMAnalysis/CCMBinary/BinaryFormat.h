#ifndef CCMAnalysis_BinaryFormat_H
#define CCMAnalysis_BinaryFormat_H

#include <string>
#include <time.h>
#include <vector>

namespace CCMAnalysis {
namespace Binary {

struct ChannelHeader {
    uint32_t version = 0;
    std::string physical_board_id;
    std::string board_serial_number;
    std::string physical_channel_type;
    std::string physical_channel_id;
    uint32_t caen_optical_link_number;
    uint32_t caen_optical_link_board_number;
    uint32_t caen_channel_number;
};

struct DigitizerBoard {
    uint32_t version = 0;
    std::string physical_board_id;
    std::string board_serial_number;
    uint32_t caen_optical_link_number;
    uint32_t caen_optical_link_board_number;
    std::vector<ChannelHeader> channels;
    std::string trigger_out_type;
};

struct CCMDAQMachineConfig {
    uint32_t version = 0;
    std::string machine_identifier;
    uint32_t num_digitizer_boards;
    uint32_t num_channels;
    uint32_t num_samples;
    uint32_t trigger_percent_after;
    uint32_t trigger_time_tolerance;
    uint32_t missed_trigger_tolerance;
};

struct CCMDAQConfig {
    uint32_t version = 0;
    std::vector<CCMDAQMachineConfig> machine_configurations;
    std::vector<DigitizerBoard> digitizer_boards;
};

struct CCMTrigger {
    uint32_t version = 0;
    /// The number of samples for each channel (should be equal to #NSAMPLES)
    std::vector<uint16_t> channel_sizes;
    /// Lets you know if the channel was masked in the hardware
    std::vector<uint16_t> channel_masks;
    /// The temperatures of each channel (updated periodically)
    std::vector<uint32_t> channel_temperatures;
    /// The event number on each board
    std::vector<uint32_t> board_event_numbers;
    /// The internal clock time for each board
    std::vector<uint32_t> board_times;
};

struct CCMTriggerReadout {
    uint32_t version = 0;
    uint32_t event_number;
    /// The computer time of the event
    struct timespec computer_time;  // needed to explicitly declare as a struct - stackoverflow 11153334
    /// The triggers within the readout window
    std::vector<CCMTrigger> triggers;
    /// The ADC value for each sample on a given channel
    ///  samples has one entry for each sample
    ///  each channel entry is a vector of ADC counts ordered in time
    std::vector<std::vector<uint16_t>> samples;
};

struct CCMData {
    uint32_t version = 0;
    CCMDAQConfig daq_config;
    std::vector<CCMTriggerReadout> trigger_readout;
};

static const char lexical_magic_number[] = "CCMVersionedBinary";
static constexpr size_t magic_size = sizeof(lexical_magic_number) - 1;

} // namespace Binary
} // namespace CCMAnalsysis

#endif // CCMAnalysis_BinaryFormat_H
