#ifndef CCMAnalysis_BinaryFormat_H
#define CCMAnalysis_BinaryFormat_H

namespace CCMAnalysis {
namespace Binary {

struct ChannelHeader {
    uint32_t version;
    std::string physical_board_id;
    std::string board_serial_number;
    std::string physical_channel_id;
    uint32_t caen_optical_link_number;
    uint32_t caen_optical_link_board_number;
    uint32_t caen_channel_number;
};

struct DigitizerBoard {
    uint32_t version;
    std::string physical_board_id;
    std::string board_serial_number;
    uint32_t caen_optical_link_number;
    uint32_t caen_optical_link_board_number;
    std::vector<ChannelHeader> channels;
};

struct CCMDAQConfig {
    uint32_t version;
    std::string machine_identifier;
    uint32_t num_digitizer_boards;
    uint32_t num_channels;
    uint32_t num_samples;
    uint32_t ring_buffer_size;
    float trigger_percent_after;
    std::vector<DigitizerBoard> digitizer_boards;
};

struct CCMDigitizerReadout {
    uint32_t version;
    /// The number of samples for each channel (should be equal to #NSAMPLES)
    std::vector<uint16_t> channel_sizes;
    /// Lets you know if the channel was masked in the hardware
    std::vector<uint16_t> channel_masks;
    /// The temperatures of each channel (updated periodically)
    std::vector<uint16_t> board_temperatures;
    /// The event number on each board
    std::vector<uint32_t> board_event_numbers;
    /// The internal clock time for each board
    std::vector<uint32_t> board_times;
    /// The ADC value for each sample on a given channel
    std::vector<std::vector<uint16_t>> samples;
};

struct CCMTriggerReadout {
    uint32_t version;
    uint32_t evNum;
    /// The digitizer information
    CCMDigitizerReadout digitizer_readout;
    /// The computer time of the event
    struct timespec computer_time;  // needed to explicitly declare as a struct - stackoverflow 11153334
};

struct CCMData {
    uint32_t version;
    CCMDAQConfig daq_config;
    std::vector<CCMTriggerReadout> trigger_readout;
};

static const char lexical_magic_number[] = "CCMVersionedBinary";
static constexpr size_t magic_size = sizeof(lexical_magic_number) - 1;

} // namespace Binary
} // namespace CCMAnalsysis

#endif // CCMAnalysis_BinaryFormat_H
