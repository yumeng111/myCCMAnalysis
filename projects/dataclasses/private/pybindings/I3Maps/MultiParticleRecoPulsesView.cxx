
#include "dataclasses/MultiParticleRecoPulsesView.h"
#include "icetray/python/dataclass_suite.hpp"

namespace bp = boost::python;

CCMRecoPulseSeriesMapPtr
underhanded_apply(const MultiParticleRecoPulsesView& mask, const I3Frame &frame)
{
	return boost::const_pointer_cast<CCMRecoPulseSeriesMap>(mask.Apply(frame));
}

void register_MultiParticleRecoPulsesView()
{
	bp::class_<MultiParticleRecoPulsesView, bp::bases<I3FrameObject>,
	    MultiParticleRecoPulsesViewPtr>("MultiParticleRecoPulsesView",
	    bp::init<const I3Frame&, const std::string &>(bp::args("frame", "key")))
		.def("apply", &underhanded_apply, "Apply the union to an I3Frame, returning an CCMRecoPulseSeries.")
		.add_property("source", &MultiParticleRecoPulsesView::GetSource)
		.def(bp::dataclass_suite<MultiParticleRecoPulsesView>())
	    ;
	;

	register_pointer_conversions<MultiParticleRecoPulsesView>();
}
