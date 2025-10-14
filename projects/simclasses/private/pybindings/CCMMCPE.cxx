//
//   Copyright (c) 2013   Alex Olivas
//

#include <simclasses/CCMMCPE.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_CCMMCPE() {
    {
    scope mcpe_scope =
        class_<CCMMCPE, boost::shared_ptr<CCMMCPE> >("CCMMCPE")
        .def(dataclass_suite<CCMMCPE>())
        //.def(init<size_t, size_t, size_t, std::vector<size_t>, WLSLocationSeries, float, float, float, float, float, float, float, I3Position, I3Direction, CCMMCPE::PhotonSource >())
        .def(init<size_t, size_t, std::vector<size_t>, WLSLocationSeries, float, float, float, float, float, CCMMCPE::PhotonSource >())
        .def_readwrite("parent_id",&CCMMCPE::parent_id)
        .def_readwrite("track_id",&CCMMCPE::track_id)
        .def_readwrite("n_photons_per_wls",&CCMMCPE::n_photons_per_wls)
        .def_readwrite("wls_loc",&CCMMCPE::wls_loc)
        .def_readwrite("time",&CCMMCPE::time)
        //.def_readwrite("calculated_time",&CCMMCPE::calculated_time)
        .def_readwrite("wavelength",&CCMMCPE::wavelength)
        .def_readwrite("distance_uv",&CCMMCPE::distance_uv)
        .def_readwrite("original_wavelength",&CCMMCPE::original_wavelength)
        .def_readwrite("distance_visible",&CCMMCPE::distance_visible)
        //.def_readwrite("calculated_distance_uv",&CCMMCPE::calculated_distance_uv)
        //.def_readwrite("calculated_distance_visible",&CCMMCPE::calculated_distance_visible)
        //.def_readwrite("position",&CCMMCPE::position)
        //.def_readwrite("direction",&CCMMCPE::direction)
        .def_readwrite("photon_source",&CCMMCPE::photon_source)
        ;
    enum_<CCMMCPE::PhotonSource>("PhotonSource")
      .value("Unknown", CCMMCPE::PhotonSource::Unknown)
      .value("Scintillation", CCMMCPE::PhotonSource::Scintillation)
      .value("Cherenkov", CCMMCPE::PhotonSource::Cherenkov)
      .value("OpWLS", CCMMCPE::PhotonSource::OpWLS)
      .export_values()
      ;

    }


    class_<CCMMCPESeries, CCMMCPESeriesPtr, bases<I3FrameObject>>("CCMMCPESeries")
        .def(dataclass_suite<CCMMCPESeries>())
        ;
    register_pointer_conversions<CCMMCPESeries>();

    class_<CCMMCPESeriesMap, CCMMCPESeriesMapPtr, bases<I3FrameObject> >("CCMMCPESeriesMap")
        .def(dataclass_suite<CCMMCPESeriesMap>())
        ;
    register_pointer_conversions<CCMMCPESeriesMap>();

    class_<CCMMCPESeriesMapByID, CCMMCPESeriesMapByIDPtr, bases<I3FrameObject> >("CCMMCPESeriesMapByID")
        .def(dataclass_suite<CCMMCPESeriesMapByID>())
        ;
    register_pointer_conversions<CCMMCPESeriesMapByID>();
}
