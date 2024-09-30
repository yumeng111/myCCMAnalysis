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
        .def(init<size_t, size_t, size_t, std::vector<size_t>, float, float, float, float, float, float, float, I3Position, I3Direction, CCMMCPE::PhotonSource >())
        .def_readwrite("parent_id",&CCMMCPE::parent_id)
        .def_readwrite("track_id",&CCMMCPE::track_id)
        .def_readwrite("n_wls",&CCMMCPE::n_wls)
        .def_readwrite("n_photons_per_wls",&CCMMCPE::n_photons_per_wls)
        .def_readwrite("g4_time",&CCMMCPE::g4_time)
        .def_readwrite("calculated_time",&CCMMCPE::calculated_time)
        .def_readwrite("wavelength",&CCMMCPE::wavelength)
        .def_readwrite("g4_distance_uv",&CCMMCPE::g4_distance_uv)
        .def_readwrite("g4_distance_visible",&CCMMCPE::g4_distance_visible)
        .def_readwrite("calculated_distance_uv",&CCMMCPE::calculated_distance_uv)
        .def_readwrite("calculated_distance_visible",&CCMMCPE::calculated_distance_visible)
        .def_readwrite("position",&CCMMCPE::position)
        .def_readwrite("direction",&CCMMCPE::direction)
        .def_readwrite("photon_source",&CCMMCPE::photon_source)
        ;
    enum_<CCMMCPE::PhotonSource>("PhotonSource")
      .value("Unknown", CCMMCPE::PhotonSource::Unknown)
      .value("Scintillation", CCMMCPE::PhotonSource::Scintillation)
      .value("Cerenkov", CCMMCPE::PhotonSource::Cerenkov)
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
}
