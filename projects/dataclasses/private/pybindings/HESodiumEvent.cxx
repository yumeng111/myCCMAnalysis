#include <vector>
#include <iostream>

#include <dataclasses/physics/HESodiumEvent.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>


namespace bp = boost::python;

std::string to_str(HESodiumEvent const & theSummary) {
    std::ostringstream oss;
    oss << theSummary << std::flush;
    return oss.str();
}


void register_HESodiumEvent() {
    std::string(*convert_to_str)(HESodiumEvent const &) = to_str;
    bp::class_<HESodiumEvent, bp::bases<I3FrameObject>, boost::shared_ptr<HESodiumEvent> >("HESodiumEvent")

    .def_readwrite("photon_vertex",&HESodiumEvent::photon_vertex)
    .def_readwrite("electron_vertex",&HESodiumEvent::electron_vertex)
    .def_readwrite("positron_vertex",&HESodiumEvent::positron_vertex)
    .def(bp::dataclass_suite<HESodiumEvent>())
    ;
    register_pointer_conversions<HESodiumEvent>();

    bp::class_<HESodiumEventSeries, HESodiumEventSeriesPtr, bp::bases<I3FrameObject>>("HESodiumEventSeries")
        .def(bp::dataclass_suite<HESodiumEventSeries>())
        ;
    register_pointer_conversions<HESodiumEventSeries>();
}
