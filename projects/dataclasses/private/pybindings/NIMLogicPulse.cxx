#include <vector>
#include <icetray/I3Units.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/I3Map.h>
#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

namespace bp = boost::python;

void register_NIMLogicPulse() {
    bp::class_<NIMLogicPulse, NIMLogicPulsePtr>("NIMLogicPulse")
    PROPERTY(NIMLogicPulse, time, NIMPulseTime)
    PROPERTY(NIMLogicPulse, length, NIMPulseLength)
    .def( bp::self == bp::self )
    .def( bp::self != bp::self )
    .def(bp::dataclass_suite<NIMLogicPulse>())
    ;

    bp::class_<NIMLogicPulseSeries, bp::bases<I3FrameObject>,
        boost::shared_ptr<NIMLogicPulseSeries>>("NIMLogicPulseSeries")
        .def(bp::dataclass_suite<NIMLogicPulseSeries>())
        ;

    register_pointer_conversions<NIMLogicPulseSeries>();

    bp::class_<NIMLogicPulseSeriesMap, bp::bases<I3FrameObject>,
        boost::shared_ptr<NIMLogicPulseSeriesMap>>("NIMLogicPulseSeriesMap")
        .def(bp::dataclass_suite<NIMLogicPulseSeriesMap>())
        ;

    register_pointer_conversions<NIMLogicPulseSeriesMap>();
}
