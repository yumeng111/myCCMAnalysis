
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime.h>
#include <icetray/python/dataclass_suite.hpp>

namespace bp = boost::python;

CCMRecoPulseSeriesMapPtr
underhanded_apply(const CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime& shifter, const I3Frame &frame)
{
  return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(shifter.Apply(frame));
}

void register_CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime()
{
  bp::class_<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime, bp::bases<I3FrameObject>,
      CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimePtr>("CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime",
      bp::init<const std::string &, const std::string &, const std::string &, const std::string &>(bp::args("pulses_key", "calibration_key", "nim_pulses_key", "geometry_key")))
    .def("apply", &underhanded_apply, "Apply the correction to an I3Frame, returning an CCMRecoPulseSeries.")
    .add_property("pulses_source", &CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime::GetPulsesSource)
    .add_property("calibration_source", &CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime::GetCalibrationSource)
    .add_property("nim_pulses_source", &CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime::GetNIMPulsesSource)
    .add_property("geometry_source", &CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime::GetGeometrySource)
    .def(bp::dataclass_suite<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>())
    ;
  ;

  register_pointer_conversions<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>();
}
