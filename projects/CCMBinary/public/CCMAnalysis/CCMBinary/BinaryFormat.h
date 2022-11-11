#ifndef CCMAnalysis_BinaryFormat_H
#define CCMAnalysis_BinaryFormat_H

#include <string>
#include <time.h>
#include <vector>
#include <iostream>

#include <icetray/I3FrameObject.h>
#include <icetray/serialization.h>

namespace CCMAnalysis {
namespace Binary {

static const unsigned channelheader_version_ = 0;
struct ChannelHeader : public I3FrameObject {
    uint32_t version = channelheader_version_;
    std::string physical_board_id;
    std::string board_serial_number;
    std::string physical_channel_type;
    std::string physical_channel_id;
    uint32_t caen_optical_link_number;
    uint32_t caen_optical_link_board_number;
    uint32_t caen_channel_number;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
    bool operator==(ChannelHeader const & other) const;
    bool operator!=(ChannelHeader const & other) const;
    std::ostream & Print(std::ostream & os) const override;
};

static const unsigned digitizerboard_version_ = 0;
struct DigitizerBoard : public I3FrameObject{
    uint32_t version = digitizerboard_version_;
    std::string physical_board_id;
    std::string board_serial_number;
    uint32_t caen_optical_link_number;
    uint32_t caen_optical_link_board_number;
    std::vector<ChannelHeader> channels;
    std::string trigger_out_type;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
    bool operator==(DigitizerBoard const & other) const;
    bool operator!=(DigitizerBoard const & other) const;
    std::ostream & Print(std::ostream & os) const override;
};

static const unsigned ccmdaqmachineconfig_version_ = 0;
struct CCMDAQMachineConfig : public I3FrameObject{
    uint32_t version = ccmdaqmachineconfig_version_;
    std::string machine_identifier;
    uint32_t num_digitizer_boards;
    uint32_t num_channels;
    uint32_t num_samples;
    uint32_t trigger_percent_after;
    uint32_t trigger_time_tolerance;
    uint32_t missed_trigger_tolerance;
    uint32_t offset_estimate_min_triggers;
    long double offset_estimate_abs_error_threshold;
    long double offset_estimate_rel_error_threshold;
    long double offset_estimate_tau;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
    bool operator==(CCMDAQMachineConfig const & other) const;
    bool operator!=(CCMDAQMachineConfig const & other) const;
    std::ostream & Print(std::ostream & os) const override;
};

static const unsigned ccmdaqconfig_version_ = 0;
struct CCMDAQConfig : public I3FrameObject {
    uint32_t version = ccmdaqconfig_version_;
    std::vector<CCMDAQMachineConfig> machine_configurations;
    std::vector<DigitizerBoard> digitizer_boards;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
    bool operator==(CCMDAQConfig const & other) const;
    bool operator!=(CCMDAQConfig const & other) const;
    std::ostream & Print(std::ostream & os) const override;
};

static const unsigned ccmtrigger_version_ = 0;
struct CCMTrigger : public I3FrameObject {
    uint32_t version = ccmtrigger_version_;
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
    /// The computer time the event was read out
    std::vector<struct timespec> board_computer_times;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
    std::ostream & Print(std::ostream & os) const override;
    bool operator==(CCMTrigger const & other) const;
    bool operator!=(CCMTrigger const & other) const;
};

static const unsigned ccmtriggerreadout_version_ = 0;
struct CCMTriggerReadout : public I3FrameObject {
    uint32_t version = ccmtriggerreadout_version_;
    uint32_t event_number;
    /// The computer time of the event
    struct timespec computer_time;  // needed to explicitly declare as a struct - stackoverflow 11153334
    /// The triggers within the readout window
    std::vector<CCMTrigger> triggers;
    /// The ADC value for each sample on a given channel
    ///  samples has one entry for each sample
    ///  each channel entry is a vector of ADC counts ordered in time
    std::vector<std::vector<uint16_t>> samples;

    /// Not saved
    /// Used only for event alignment
    bool stop = false;
    CCMTriggerReadout() = default;
    CCMTriggerReadout(CCMTriggerReadout const & other) :
        version(other.version),
        event_number(other.event_number),
        computer_time(other.computer_time),
        triggers(other.triggers),
        samples(other.samples),
        stop(other.stop) {}
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
    std::ostream & Print(std::ostream & os) const override;
};

static const unsigned ccmdata_version_ = 0;
struct CCMData : public I3FrameObject {
    uint32_t version = ccmdata_version_;
    CCMDAQConfig daq_config;
    std::vector<CCMTriggerReadout> trigger_readout;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

static const char lexical_magic_number[] = "CCMVersionedBinary";
static constexpr size_t magic_size = sizeof(lexical_magic_number) - 1;

I3_POINTER_TYPEDEFS(ChannelHeader);
I3_POINTER_TYPEDEFS(DigitizerBoard);
I3_POINTER_TYPEDEFS(CCMDAQMachineConfig);
I3_POINTER_TYPEDEFS(CCMDAQConfig);
I3_POINTER_TYPEDEFS(CCMTrigger);
I3_POINTER_TYPEDEFS(CCMTriggerReadout);
I3_POINTER_TYPEDEFS(CCMData);
} // namespace Binary
} // namespace CCMAnalsysis

std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::ChannelHeader c);
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::DigitizerBoard c);
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMDAQMachineConfig c);
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMDAQConfig c);
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMTrigger c);
std::ostream& operator<<(std::ostream& os, const CCMAnalysis::Binary::CCMTriggerReadout c);

I3_CLASS_VERSION(CCMAnalysis::Binary::ChannelHeader, CCMAnalysis::Binary::channelheader_version_);
I3_CLASS_VERSION(CCMAnalysis::Binary::DigitizerBoard, CCMAnalysis::Binary::digitizerboard_version_);
I3_CLASS_VERSION(CCMAnalysis::Binary::CCMDAQMachineConfig, CCMAnalysis::Binary::ccmdaqmachineconfig_version_);
I3_CLASS_VERSION(CCMAnalysis::Binary::CCMDAQConfig, CCMAnalysis::Binary::ccmdaqconfig_version_);
I3_CLASS_VERSION(CCMAnalysis::Binary::CCMTrigger, CCMAnalysis::Binary::ccmtrigger_version_);
I3_CLASS_VERSION(CCMAnalysis::Binary::CCMTriggerReadout, CCMAnalysis::Binary::ccmtriggerreadout_version_);
I3_CLASS_VERSION(CCMAnalysis::Binary::CCMData, CCMAnalysis::Binary::ccmdata_version_);

#endif // CCMAnalysis_BinaryFormat_H
