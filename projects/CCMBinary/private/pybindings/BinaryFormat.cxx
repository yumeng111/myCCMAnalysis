#include <CCMBinary/BinaryFormat.h>
#include <CCMBinary/load_project.h>
#include <boost/python.hpp>

using namespace boost::python

I3_PYTHON_MODULE(CCMBinary)
{
  load_project("CCMBinary", false);
  import("icecube.icetray");

  scope().attr("channelheader_version_")=channelheader_version_;
  class_<ChannelHeader,bases<I3FrameObject> >("ChannelHeader",
    "Container class for channel header information")
    .def_readonly("version",&ChannelHeader::version)
    .def_readwrite("physical_board_id",&ChannelHeader::physical_board_id)
    .def_readwrite("board_serial_number",&ChannelHeader::board_serial_number)
    .def_readwrite("physical_channel_type",&ChannelHeader::physical_channel_type)
    .def_readwrite("physical_channel_id",&ChannelHeader::physical_channel_id)
    .def_readwrite("caen_optical_link_number",&ChannelHeader::caen_optical_link_number)
    .def_readwrite("caen_optical_link_board_number",&ChannelHeader::caen_optical_link_board_number)
    .def_readwrite("caen_channel_number",&ChannelHeader::caen_channel_number);
  
  scope().attr("digitizerboard_version_")=digitizerboard_version_;
  class_<DigitizerBoard,bases<I3FrameObject> >("DigitizerBoard",
    "Container class for digitizer board information")
    .def_readonly("version",&DigitizerBoard::version)
    .def_readwrite("physical_board_id",&DigitizerBoard::physical_board_id)
    .def_readwrite("board_serial_number",&DigitizerBoard::board_serial_number)
    .def_readwrite("caen_optical_link_number",&DigitizerBoard::caen_optical_link_number)
    .def_readwrite("caen_optical_link_board_number",&DigitizerBoard::caen_optical_link_board_number)
    .def_readwrite("channels",&DigitizerBoard::channels);
  
  scope().attr("ccmdaqmachineconfig_version_")=ccmdaqmachineconfig_version_;
  class_<CCMDAQMachineConfig,bases<I3FrameObject> >("CCMDAQMachineConfig",
    "Container class for DAQ machine information")
    .def_readonly("version",&CCMDAQMachineConfig::version)
    .def_readwrite("num_digitizer_boards",&CCMDAQMachineConfig::num_digitizer_boards)
    .def_readwrite("num_channels",&CCMDAQMachineConfig::num_channels)
    .def_readwrite("num_samples",&CCMDAQMachineConfig::num_samples)
    .def_readwrite("trigger_percent_after",&CCMDAQMachineConfig::trigger_percent_after)
    .def_readwrite("trigger_time_tolerance",&CCMDAQMachineConfig::trigger_time_tolerance)
    .def_readwrite("missed_trigger_tolerance",&CCMDAQMachineConfig::missed_trigger_tolerance)
    .def_readwrite("offset_estimate_min_triggers",&CCMDAQMachineConfig::offset_estimate_min_triggers)
    .def_readwrite("offset_estimate_abs_error_threshold",&CCMDAQMachineConfig::offset_estimate_abs_error_threshold)
    .def_readwrite("offset_estimate_rel_error_threshold",&CCMDAQMachineConfig::offset_estimate_rel_error_threshold)
    .def_readwrite("offset_estimate_tau",&CCMDAQMachineConfig::offset_estimate_tau);
  
  scope().attr("ccmdaqconfig_version_")=ccmdaqconfig_version_;
  class_<CCMDAQConfig,bases<I3FrameObject> >("CCMDAQConfig",
    "Container class for DAQ config information")
    .def_readonly("version",&CCMDAQConfig::version)
    .def_readwrite("machine_configurations",&CCMDAQConfig::machine_configurations)
    .def_readwrite("digitizer_boards",&CCMDAQConfig::digitizer_boards);
  
  scope().attr("ccmtrigger_version_")=ccmtrigger_version_;
  class_<CCMTrigger,bases<I3FrameObject> >("CCMTrigger",
    "Container class for trigger information")
    .def_readonly("version",&CCMTrigger::version)
    .def_readwrite("channel_sizes",&CCMTrigger::channel_sizes)
    .def_readwrite("channel_masks",&CCMTrigger::channel_masks)
    .def_readwrite("channel_temperatures",&CCMTrigger::channel_temperatures)
    .def_readwrite("board_event_numbers",&CCMTrigger::board_event_numbers)
    .def_readwrite("board_times",&CCMTrigger::board_times)
    .def_readwrite("board_computer_times",&CCMTrigger::board_computer_times);
  
  scope().attr("ccmdata_version_")=ccmdata_version_;
  class_<CCMData,bases<I3FrameObject> >("CCMData",
    "Container class for CCM data information")
    .def_readonly("version",&CCMData::version)
    .def_readonly("daq_config",&CCMData::daq_config)
    .def_readonly("trigger_readout",&CCMData::trigger_readout);

}
