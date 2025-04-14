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
        .def(init<std::vector<size_t>, WLSLocationSeries, float, float, float>())
        .def_readwrite("n_photons_per_wls",&LightweightCCMMCPE::n_photons_per_wls)
        .def_readwrite("wls_loc",&LightweightCCMMCPE::wls_loc)
        .def_readwrite("time",&LightweightCCMMCPE::time)
        .def_readwrite("wavelength",&LightweightCCMMCPE::wavelength)
        .def_readwrite("distance_uv",&LightweightCCMMCPE::distance_uv)
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
