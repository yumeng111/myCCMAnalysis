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
        .def(init<float, float, I3Position, I3Direction, CCMMCPE::PhotonSource >())
        .def_readwrite("time",&CCMMCPE::time)
        .def_readwrite("wavelength",&CCMMCPE::wavelength)
        .def_readwrite("position",&CCMMCPE::position)
        .def_readwrite("direction",&CCMMCPE::direction)
        .def_readwrite("photon_source",&CCMMCPE::photon_source)
        ;
    enum_<CCMMCPE::PhotonSource>("PhotonSource")
      .value("Unknown", CCMMCPE::PhotonSource::Unknown)
      .value("Scintillation", CCMMCPE::PhotonSource::Scintillation)
      .value("Cerenkov", CCMMCPE::PhotonSource::Cerenkov)
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
