
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusBeamTime.h>
#include <icetray/python/dataclass_suite.hpp>

namespace bp = boost::python;

CCMRecoPulseSeriesMapPtr
underhanded_apply(const CCMRecoPulseSeriesMapApplySPECalPlusBeamTime& shifter, const I3Frame &frame)
{
  return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(shifter.Apply(frame));
}

void register_CCMRecoPulseSeriesMapApplySPECalPlusBeamTime()
{
  bp::class_<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime, bp::bases<I3FrameObject>,
      CCMRecoPulseSeriesMapApplySPECalPlusBeamTimePtr>("CCMRecoPulseSeriesMapApplySPECalPlusBeamTime",
      bp::init<const std::string &, const std::string &, const std::string &, const std::string &, const std::string &>(bp::args("pulses_key", "calibration_key", "nim_pulses_key", "geometry_key", "bcm_summary_key")))
    .def("apply", &underhanded_apply, "Apply the correction to an I3Frame, returning an CCMRecoPulseSeries.")
    .add_property("pulses_source", &CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::GetPulsesSource)
    .add_property("calibration_source", &CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::GetCalibrationSource)
    .add_property("nim_pulses_source", &CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::GetNIMPulsesSource)
    .add_property("geometry_source", &CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::GetGeometrySource)
    .add_property("bcm_summary_source", &CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::GetBCMSummarySource)
    .def(bp::dataclass_suite<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime>())
    ;
  ;

  register_pointer_conversions<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime>();
}
