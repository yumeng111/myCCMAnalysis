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
        .def(init<float, float, float, float, float, float, size_t, std::vector<size_t>, WLSLocationSeries, PhotonSummary::PhotonSource, PhotonSummary::PhotonSource, PhotonSummary::PhotonSource >())
        .def_readwrite("g4_distance_uv",&PhotonSummary::g4_distance_uv)
        .def_readwrite("g4_distance_visible",&PhotonSummary::g4_distance_visible)
        .def_readwrite("calculated_distance_uv",&PhotonSummary::calculated_distance_uv)
        .def_readwrite("calculated_distance_visible",&PhotonSummary::calculated_distance_visible)
        .def_readwrite("g4_time",&PhotonSummary::g4_time)
        .def_readwrite("calculated_time",&PhotonSummary::calculated_time)
        .def_readwrite("n_wls",&PhotonSummary::n_wls)
        .def_readwrite("n_photons_per_wls",&PhotonSummary::n_photons_per_wls)
        .def_readwrite("wls_loc",&PhotonSummary::wls_loc)
        .def_readwrite("photon_source",&PhotonSummary::photon_source)
        .def_readwrite("temp_parent",&PhotonSummary::temp_parent)
        .def_readwrite("current_process",&PhotonSummary::current_process)
        ;
    enum_<PhotonSummary::PhotonSource>("PhotonSource")
      .value("Unknown", PhotonSummary::PhotonSource::Unknown)
      .value("Scintillation", PhotonSummary::PhotonSource::Scintillation)
      .value("Cerenkov", PhotonSummary::PhotonSource::Cerenkov)
      .value("OpWLS", PhotonSummary::PhotonSource::OpWLS)
      .export_values()
      ;
    //enum_<PhotonSummary::WLSLocation>("WLSLocation")
    //  .value("Unknown", PhotonSummary::WLSLocation::Unknown)
    //  .value("PMT", PhotonSummary::WLSLocation::PMT)
    //  .value("FoilTop", PhotonSummary::WLSLocation::FoilTop)
    //  .value("FoilBottom", PhotonSummary::WLSLocation::FoilBottom)
    //  .value("FoilSides", PhotonSummary::WLSLocation::FoilSides)
    //  .export_values()
    //  ;
    }


    class_<PhotonSummarySeries, PhotonSummarySeriesPtr, bases<I3FrameObject>>("PhotonSummarySeries")
        .def(dataclass_suite<PhotonSummarySeries>())
        ;
    register_pointer_conversions<PhotonSummarySeries>();

}

