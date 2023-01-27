
#include <dataclasses/I3MapCCMPMTKeyUnion.h>
#include <icetray/python/dataclass_suite.hpp>

namespace bp = boost::python;

CCMRecoPulseSeriesMapPtr
underhanded_apply(const CCMRecoPulseSeriesMapUnion& mask, const I3Frame &frame)
{
	return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(mask.Apply(frame));
}

void register_CCMRecoPulseSeriesMapUnion()
{
	bp::class_<CCMRecoPulseSeriesMapUnion, bp::bases<I3FrameObject>,
	    CCMRecoPulseSeriesMapUnionPtr>("CCMRecoPulseSeriesMapUnion",
	    bp::init<const I3Frame&, const std::vector<std::string> &>(bp::args("frame", "keys")))
		.def("apply", &underhanded_apply, "Apply the union to an I3Frame, returning an CCMRecoPulseSeries.")
		.add_property("sources", &CCMRecoPulseSeriesMapUnion::GetSources)
		.def(bp::dataclass_suite<CCMRecoPulseSeriesMapUnion>())
	    ;
	;
	
	register_pointer_conversions<CCMRecoPulseSeriesMapUnion>();
}
