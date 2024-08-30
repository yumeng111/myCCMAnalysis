//
//   Copyright (c) 2013   Alex Olivas
//

#include <simclasses/LightweightCCMMCPE.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

void register_LightweightCCMMCPE() {
    {
    scope mcpe_scope =
        class_<LightweightCCMMCPE, boost::shared_ptr<LightweightCCMMCPE> >("LightweightCCMMCPE")
        .def(dataclass_suite<LightweightCCMMCPE>())
        .def(init<float, float, float, LightweightCCMMCPE::PhotonSource >())
        .def_readwrite("g4_time",&LightweightCCMMCPE::g4_time)
        .def_readwrite("wavelength",&LightweightCCMMCPE::wavelength)
        .def_readwrite("g4_distance_uv",&LightweightCCMMCPE::g4_distance_uv)
        .def_readwrite("photon_source",&LightweightCCMMCPE::photon_source)
        ;
    enum_<LightweightCCMMCPE::PhotonSource>("PhotonSource")
      .value("Unknown", LightweightCCMMCPE::PhotonSource::Unknown)
      .value("Scintillation", LightweightCCMMCPE::PhotonSource::Scintillation)
      .value("Cerenkov", LightweightCCMMCPE::PhotonSource::Cerenkov)
      .value("OpWLS", LightweightCCMMCPE::PhotonSource::OpWLS)
      .export_values()
      ;
        
    }


    class_<LightweightCCMMCPESeries, LightweightCCMMCPESeriesPtr, bases<I3FrameObject>>("LightweightCCMMCPESeries")
        .def(dataclass_suite<LightweightCCMMCPESeries>())
        ;
    register_pointer_conversions<LightweightCCMMCPESeries>();

    class_<LightweightCCMMCPESeriesMap, LightweightCCMMCPESeriesMapPtr, bases<I3FrameObject> >("LightweightCCMMCPESeriesMap")
        .def(dataclass_suite<LightweightCCMMCPESeriesMap>())
        ;
    register_pointer_conversions<LightweightCCMMCPESeriesMap>();
}
