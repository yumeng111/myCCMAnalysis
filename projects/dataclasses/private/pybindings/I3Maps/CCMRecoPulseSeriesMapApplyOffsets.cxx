
#include <dataclasses/CCMRecoPulseSeriesMapApplyOffsets.h>
#include <icetray/python/dataclass_suite.hpp>

namespace bp = boost::python;

CCMRecoPulseSeriesMapPtr
underhanded_apply(const CCMRecoPulseSeriesMapApplyOffsets& shifter, const I3Frame &frame)
{
  return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(shifter.Apply(frame));
}

void register_CCMRecoPulseSeriesMapApplyOffsets()
{
  bp::class_<CCMRecoPulseSeriesMapApplyOffsets, bp::bases<I3FrameObject>,
      CCMRecoPulseSeriesMapApplyOffsetsPtr>("CCMRecoPulseSeriesMapApplyOffsets",
      bp::init<const std::string &, const std::string &>(bp::args("pulses_key", "offsets_key")))
    .def("apply", &underhanded_apply, "Apply the correction to an I3Frame, returning an CCMRecoPulseSeries.")
    .add_property("pulses_source", &CCMRecoPulseSeriesMapApplyOffsets::GetPulsesSource)
    .add_property("offsets_source", &CCMRecoPulseSeriesMapApplyOffsets::GetOffsetsSource)
    .def(bp::dataclass_suite<CCMRecoPulseSeriesMapApplyOffsets>())
    ;
  ;

  register_pointer_conversions<CCMRecoPulseSeriesMapApplyOffsets>();
}
