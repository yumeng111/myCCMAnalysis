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
        .def(init<size_t, size_t, float, float, float, float, float, I3Position, I3Direction, CCMMCPE::PhotonSource >())
        .def_readwrite("parent_id",&CCMMCPE::parent_id)
        .def_readwrite("track_id",&CCMMCPE::track_id)
        .def_readwrite("global_time",&CCMMCPE::global_time)
        .def_readwrite("local_time",&CCMMCPE::local_time)
        .def_readwrite("wavelength",&CCMMCPE::wavelength)
        .def_readwrite("distance_uv",&CCMMCPE::distance_uv)
        .def_readwrite("distance_visible",&CCMMCPE::distance_visible)
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
