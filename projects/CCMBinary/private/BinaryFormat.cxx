#ifndef CCMAnalysis_BinaryFormat_CXX
#define CCMAnalysis_BinaryFormat_CXX

#include <string>
#include <vector>
#include <cstring>
#include <sstream>
#include <iostream>

#include <icetray/I3Logging.h>
#include <icetray/serialization.h>

#include "CCMAnalysis/CCMBinary/BinaryFormat.h"

namespace CCMAnalysis {
namespace Binary {

template <class Archive>
void
ChannelHeader::serialize(Archive& ar, unsigned version) {
    if(version > channelheader_version_)
        log_fatal("Attempting to read version %u from file but running version %u of ChannelHeader class.", version, channelheader_version_);

    ar & make_nvp("physical_board_id", physical_board_id);
    ar & make_nvp("board_serial_number", board_serial_number);
    ar & make_nvp("physical_channel_type", physical_channel_type);
    ar & make_nvp("physical_channel_id", physical_channel_id);
    ar & make_nvp("caen_optical_link_number", caen_optical_link_number);
    ar & make_nvp("caen_optical_link_board_number", caen_optical_link_board_number);
    ar & make_nvp("caen_channel_number", caen_channel_number);
}

template <class Archive>
void
DigitizerBoard::serialize(Archive& ar, unsigned version) {
    if(version > digitizerboard_version_)
        log_fatal("Attempting to read version %u from file but running version %u of DigitizerBoard class.", version, digitizerboard_version_);
    ar & make_nvp("physical_board_id", physical_board_id);
    ar & make_nvp("board_serial_number", board_serial_number);
    ar & make_nvp("caen_optical_link_number", caen_optical_link_number);
    ar & make_nvp("caen_optical_link_board_number", caen_optical_link_board_number);
    ar & make_nvp("channels", channels);
    ar & make_nvp("trigger_out_type", trigger_out_type);
}

template <class Archive>
void
CCMDAQMachineConfig::serialize(Archive& ar, unsigned version) {
    if(version > ccmdaqmachineconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMDAQMachineConfig class.", version, ccmdaqmachineconfig_version_);
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

template <class Archive>
void
CCMDAQConfig::serialize(Archive& ar, unsigned version) {
    if(version > ccmdaqconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMDAQConfig class.", version, ccmdaqconfig_version_);
        ar & make_nvp("machine_configurations", machine_configurations);
        ar & make_nvp("digitizer_boards", digitizer_boards);
}

template <class Archive>
void
CCMTrigger::serialize(Archive& ar, unsigned version) {
    if(version > ccmtrigger_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMTrigger class.", version, ccmtrigger_version_);
    ar & make_nvp("channel_sizes", channel_sizes);
    ar & make_nvp("channel_masks", channel_masks);
    ar & make_nvp("channel_temperatures", channel_temperatures);
    ar & make_nvp("board_event_numbers", board_event_numbers);
    ar & make_nvp("board_times", board_times);
    ar & make_nvp("board_computer_times", board_computer_times);
}

template <class Archive>
void
CCMTriggerReadout::serialize(Archive& ar, unsigned version) {
    if(version > ccmtriggerreadout_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMTriggerReadout class.", version, ccmtriggerreadout_version_);
    ar & make_nvp("event_number", event_number);
    ar & make_nvp("computer_time", computer_time);
    ar & make_nvp("triggers", triggers);
    ar & make_nvp("samples", samples);
}

template <class Archive>
void
CCMData::serialize(Archive& ar, unsigned version) {
    if(version > ccmdata_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMData class.", version, ccmdata_version_);
    ar & make_nvp("daq_config", daq_config);
    ar & make_nvp("trigger_readout", trigger_readout);
}

} // namespace Binary
} // namespace CCMAnalsysis

I3_SERIALIZABLE(CCMAnalysis::Binary::ChannelHeader);
I3_SERIALIZABLE(CCMAnalysis::Binary::DigitizerBoard);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMDAQMachineConfig);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMDAQConfig);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMTrigger);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMTriggerReadout);
I3_SERIALIZABLE(CCMAnalysis::Binary::CCMData);

#endif // CCMAnalysis_BinaryFormat_CXX
