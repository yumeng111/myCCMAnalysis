#include <time.h>

#include <CCMAnalysis/CCMBinary/BinaryFormat.h>
#include <CCMAnalysis/CCMBinary/BinaryUtilities.h>
#include <icetray/load_project.h>
#include <boost/python.hpp>
// #include <icetray/scratch.h>
#include "I3Vectors.h"

using namespace boost::python;

bool operator==(timespec const & a, timespec const & b) {
    return (a.tv_sec == b.tv_sec) and (a.tv_nsec == b.tv_nsec);
}

I3_PYTHON_MODULE(CCMBinary)
{
  load_project("CCMBinary", false);
  import("icecube.icetray");

  scope().attr("channelheader_version_") = CCMAnalysis::Binary::channelheader_version_;
  class_<CCMAnalysis::Binary::ChannelHeader, boost::shared_ptr<CCMAnalysis::Binary::ChannelHeader>, bases<I3FrameObject> >("ChannelHeader",
    "Container class for channel header information")
    .def_readonly("version",&CCMAnalysis::Binary::ChannelHeader::version)
    .def_readwrite("physical_board_id",&CCMAnalysis::Binary::ChannelHeader::physical_board_id)
    .def_readwrite("board_serial_number",&CCMAnalysis::Binary::ChannelHeader::board_serial_number)
    .def_readwrite("physical_channel_type",&CCMAnalysis::Binary::ChannelHeader::physical_channel_type)
    .def_readwrite("physical_channel_id",&CCMAnalysis::Binary::ChannelHeader::physical_channel_id)
    .def_readwrite("caen_optical_link_number",&CCMAnalysis::Binary::ChannelHeader::caen_optical_link_number)
    .def_readwrite("caen_optical_link_board_number",&CCMAnalysis::Binary::ChannelHeader::caen_optical_link_board_number)
    .def_readwrite("caen_channel_number",&CCMAnalysis::Binary::ChannelHeader::caen_channel_number);

  class_<std::vector<CCMAnalysis::Binary::ChannelHeader> >("vector_ChannelHeader")
  .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::ChannelHeader> >());

  scope().attr("digitizerboard_version_") = CCMAnalysis::Binary::digitizerboard_version_;
  class_<CCMAnalysis::Binary::DigitizerBoard,bases<I3FrameObject> >("DigitizerBoard",
    "Container class for digitizer board information")
    .def_readonly("version",&CCMAnalysis::Binary::DigitizerBoard::version)
    .def_readwrite("physical_board_id",&CCMAnalysis::Binary::DigitizerBoard::physical_board_id)
    .def_readwrite("board_serial_number",&CCMAnalysis::Binary::DigitizerBoard::board_serial_number)
    .def_readwrite("caen_optical_link_number",&CCMAnalysis::Binary::DigitizerBoard::caen_optical_link_number)
    .def_readwrite("caen_optical_link_board_number",&CCMAnalysis::Binary::DigitizerBoard::caen_optical_link_board_number)
    .def_readwrite("channels",&CCMAnalysis::Binary::DigitizerBoard::channels);

  class_<std::vector<CCMAnalysis::Binary::DigitizerBoard> >("vector_DigitizerBoard")
    .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::DigitizerBoard>, true>())
    ;
  from_python_sequence<std::vector<CCMAnalysis::Binary::DigitizerBoard>, variable_capacity_policy>();

  scope().attr("ccmdaqmachineconfig_version_") = CCMAnalysis::Binary::ccmdaqmachineconfig_version_;
  class_<CCMAnalysis::Binary::CCMDAQMachineConfig,bases<I3FrameObject> >("CCMDAQMachineConfig",
    "Container class for DAQ machine information")
    .def_readonly("version",&CCMAnalysis::Binary::CCMDAQMachineConfig::version)
    .def_readwrite("num_digitizer_boards",&CCMAnalysis::Binary::CCMDAQMachineConfig::num_digitizer_boards)
    .def_readwrite("num_channels",&CCMAnalysis::Binary::CCMDAQMachineConfig::num_channels)
    .def_readwrite("num_samples",&CCMAnalysis::Binary::CCMDAQMachineConfig::num_samples)
    .def_readwrite("trigger_percent_after",&CCMAnalysis::Binary::CCMDAQMachineConfig::trigger_percent_after)
    .def_readwrite("trigger_time_tolerance",&CCMAnalysis::Binary::CCMDAQMachineConfig::trigger_time_tolerance)
    .def_readwrite("missed_trigger_tolerance",&CCMAnalysis::Binary::CCMDAQMachineConfig::missed_trigger_tolerance)
    .def_readwrite("offset_estimate_min_triggers",&CCMAnalysis::Binary::CCMDAQMachineConfig::offset_estimate_min_triggers)
    .def_readwrite("offset_estimate_abs_error_threshold",&CCMAnalysis::Binary::CCMDAQMachineConfig::offset_estimate_abs_error_threshold)
    .def_readwrite("offset_estimate_rel_error_threshold",&CCMAnalysis::Binary::CCMDAQMachineConfig::offset_estimate_rel_error_threshold)
    .def_readwrite("offset_estimate_tau",&CCMAnalysis::Binary::CCMDAQMachineConfig::offset_estimate_tau);

  class_<std::vector<CCMAnalysis::Binary::CCMDAQMachineConfig> >("vector_CCMDAQMachineConfig")
  .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::CCMDAQMachineConfig> >());

  scope().attr("ccmdaqconfig_version_") = CCMAnalysis::Binary::ccmdaqconfig_version_;
  class_<CCMAnalysis::Binary::CCMDAQConfig,bases<I3FrameObject> >("CCMDAQConfig",
    "Container class for DAQ config information")
    .def_readonly("version",&CCMAnalysis::Binary::CCMDAQConfig::version)
    .def_readwrite("machine_configurations",&CCMAnalysis::Binary::CCMDAQConfig::machine_configurations)
    .def_readwrite("digitizer_boards",&CCMAnalysis::Binary::CCMDAQConfig::digitizer_boards);

  class_<std::vector<CCMAnalysis::Binary::CCMDAQConfig> >("vector_CCMDAQConfig")
      .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::CCMDAQConfig>, true>())
      ;
  from_python_sequence<std::vector<CCMAnalysis::Binary::CCMDAQConfig>, variable_capacity_policy>();

  scope().attr("ccmtrigger_version_") = CCMAnalysis::Binary::ccmtrigger_version_;
  class_<CCMAnalysis::Binary::CCMTrigger,bases<I3FrameObject> >("CCMTrigger",
    "Container class for trigger information")
    .def_readonly("version",&CCMAnalysis::Binary::CCMTrigger::version)
    .def_readwrite("channel_sizes",&CCMAnalysis::Binary::CCMTrigger::channel_sizes)
    .def_readwrite("channel_masks",&CCMAnalysis::Binary::CCMTrigger::channel_masks)
    .def_readwrite("channel_temperatures",&CCMAnalysis::Binary::CCMTrigger::channel_temperatures)
    .def_readwrite("board_event_numbers",&CCMAnalysis::Binary::CCMTrigger::board_event_numbers)
    .def_readwrite("board_times",&CCMAnalysis::Binary::CCMTrigger::board_times)
    .def_readwrite("board_computer_times",&CCMAnalysis::Binary::CCMTrigger::board_computer_times);

  class_<std::vector<CCMAnalysis::Binary::CCMTrigger> >("vector_CCMTrigger")
  .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::CCMTrigger> >());

  scope().attr("ccmtriggerreadout_version_") = CCMAnalysis::Binary::ccmtriggerreadout_version_;
  class_<CCMAnalysis::Binary::CCMTriggerReadout, boost::shared_ptr<CCMAnalysis::Binary::CCMTriggerReadout>, bases<I3FrameObject> >("CCMTriggerReadout",
    "Container class for CCM trigger readout information")
    .def(init<>())
    .def(init<CCMAnalysis::Binary::CCMTriggerReadout>())
    .def_readonly("version",&CCMAnalysis::Binary::CCMTriggerReadout::version)
    .def_readonly("event_number",&CCMAnalysis::Binary::CCMTriggerReadout::event_number)
    .def_readonly("triggers",&CCMAnalysis::Binary::CCMTriggerReadout::triggers)
    .def_readonly("samples",&CCMAnalysis::Binary::CCMTriggerReadout::samples);

  class_<std::vector<CCMAnalysis::Binary::CCMTriggerReadout> >("vector_CCMTriggerReadout")
  .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::CCMTriggerReadout> >());

  scope().attr("ccmdata_version_") = CCMAnalysis::Binary::ccmdata_version_;
  class_<CCMAnalysis::Binary::CCMData,bases<I3FrameObject> >("CCMData",
    "Container class for CCM data information")
    .def_readonly("version",&CCMAnalysis::Binary::CCMData::version)
    .def_readonly("daq_config",&CCMAnalysis::Binary::CCMData::daq_config)
    .def_readonly("trigger_readout",&CCMAnalysis::Binary::CCMData::trigger_readout);

  class_<std::vector<CCMAnalysis::Binary::CCMData> >("vector_CCMData")
  .def(vector_indexing_suite<std::vector<CCMAnalysis::Binary::CCMData> >());

  // Define vectors of types
  class_<std::vector<uint16_t> >("vector_uint16_t")
  .def(vector_indexing_suite<std::vector<uint16_t> >());

  class_<std::vector<std::vector<uint16_t>> >("vector_vector_uint16_t")
  .def(vector_indexing_suite<std::vector<std::vector<uint16_t>> >());

  class_<std::vector<uint32_t> >("vector_uint32_t")
  .def(vector_indexing_suite<std::vector<uint32_t> >());

  class_<timespec>("timespec",
    "std::timespec")
    .def_readwrite("tv_sec",&timespec::tv_sec)
    .def_readwrite("tv_nsec",&timespec::tv_nsec);

  class_<std::vector<timespec> >("vector_timespec")
    .def(vector_indexing_suite<std::vector<timespec>, true>())
    ;
  from_python_sequence<std::vector<timespec>, variable_capacity_policy>();
    
  void (*merge1)(CCMAnalysis::Binary::CCMTriggerReadout &, CCMAnalysis::Binary::CCMTriggerReadout const &, size_t, size_t, size_t, bool) = &CCMAnalysis::Binary::merge_triggers;
  // CCMAnalysis::Binary::CCMTriggerReadout&(*merge2)(CCMAnalysis::Binary::CCMTriggerReadout &, CCMAnalysis::Binary::CCMTriggerReadout &&, size_t, size_t, size_t, bool) = &CCMAnalysis::Binary::merge_triggers;
  def("merge_triggers", merge1);

  void (*merge_empty)(CCMAnalysis::Binary::CCMTriggerReadout &, size_t, bool) = &CCMAnalysis::Binary::merge_empty_trigger;
  def("merge_empty_trigger", merge_empty);

  register_i3vector_of<I3Vector<uint16_t>>("I3VectorUInt16");
  register_i3vector_of<CCMAnalysis::Binary::CCMTrigger>("CCMTrigger");

}
