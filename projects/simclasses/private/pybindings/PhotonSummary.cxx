#include <vector>
#include <iostream>

#include <simclasses/PhotonSummary.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

void register_PhotonSummary() {
    {
    scope mcpe_scope =
        class_<PhotonSummary, boost::shared_ptr<PhotonSummary> >("PhotonSummary")
        .def(dataclass_suite<PhotonSummary>())
        .def(init<float, float, float, float, float, float, size_t >())
        .def_readwrite("g4_distance_uv",&PhotonSummary::g4_distance_uv)
        .def_readwrite("g4_distance_visible",&PhotonSummary::g4_distance_visible)
        .def_readwrite("calculated_distance_uv",&PhotonSummary::calculated_distance_uv)
        .def_readwrite("calculated_distance_visible",&PhotonSummary::calculated_distance_visible)
        .def_readwrite("g4_time",&PhotonSummary::g4_time)
        .def_readwrite("calculated_time",&PhotonSummary::calculated_time)
        .def_readwrite("n_wls",&PhotonSummary::n_wls)
        ;
    }


    class_<PhotonSummarySeries, PhotonSummarySeriesPtr, bases<I3FrameObject>>("PhotonSummarySeries")
        .def(dataclass_suite<PhotonSummarySeries>())
        ;
    register_pointer_conversions<PhotonSummarySeries>();

}

