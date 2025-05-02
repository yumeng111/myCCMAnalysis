//
//   Copyright (c) 2013   Alex Olivas
//

#include <simclasses/SampledRecoPulse.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_SampledRecoPulse() {
    {
    scope mcpe_scope =
        class_<SampledRecoPulse, boost::shared_ptr<SampledRecoPulse> >("SampledRecoPulse")
        .def(dataclass_suite<SampledRecoPulse>())
        .def(init<float, float, float, float, bool>())
        .def_readwrite("light_time",&SampledRecoPulse::light_time)
        .def_readwrite("afterpulse_time",&SampledRecoPulse::afterpulse_time)
        .def_readwrite("reco_time",&SampledRecoPulse::reco_time)
        .def_readwrite("amplitude",&SampledRecoPulse::amplitude)
        .def_readwrite("in_gev",&SampledRecoPulse::in_gev)
        ;
    }


    class_<SampledRecoPulseSeries, SampledRecoPulseSeriesPtr, bases<I3FrameObject>>("SampledRecoPulseSeries")
        .def(dataclass_suite<SampledRecoPulseSeries>())
        ;
    register_pointer_conversions<SampledRecoPulseSeries>();

    class_<SampledRecoPulseSeriesMap, SampledRecoPulseSeriesMapPtr, bases<I3FrameObject> >("SampledRecoPulseSeriesMap")
        .def(dataclass_suite<SampledRecoPulseSeriesMap>())
        ;
    register_pointer_conversions<SampledRecoPulseSeriesMap>();
}
