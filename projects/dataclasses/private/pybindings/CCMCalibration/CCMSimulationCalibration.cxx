//
//   Copyright (c) 2013   Alex Olivas
//

#include <dataclasses/calibration/CCMSimulationCalibration.h>
#include <icetray/python/dataclass_suite.hpp>

using namespace boost::python;

template<typename theType>
std::string to_str(const theType& theStr){
    std::ostringstream oss;
    oss << theStr << std::flush;
    return oss.str();
}

std::string to_str_CCMSimulationCalibration(const CCMSimulationCalibration& self){
    return to_str<CCMSimulationCalibration>(self); };

void register_CCMSimulationCalibration() {
{
    class_<CCMSimulationCalibration, bases<I3FrameObject>, boost::shared_ptr<CCMSimulationCalibration> >("CCMSimulationCalibration")
        .def(copy_suite<CCMSimulationCalibration>())
        .def(dataclass_suite<CCMSimulationCalibration>())
#define I3DOMCALPROPS (PMTEfficiencies)(LatePulseMu)(LatePulseSigma)(LatePulseScale)(PMTSPEMu)(PMTSPESigma)(Rs)(Rt)(TauS)(TauT)(TauOther)
        BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, CCMSimulationCalibration, I3DOMCALPROPS)
#undef I3DOMCALPROPS
        .def_readwrite("PMTEfficiencies",&CCMSimulationCalibration::PMTEfficiencies)
        .def_readwrite("LatePulseMu",&CCMSimulationCalibration::LatePulseMu)
        .def_readwrite("LatePulseSigma",&CCMSimulationCalibration::LatePulseSigma)
        .def_readwrite("LatePulseScale",&CCMSimulationCalibration::LatePulseScale)
        .def_readwrite("PMTSPEMu",&CCMSimulationCalibration::PMTSPEMu)
        .def_readwrite("PMTSPESigma",&CCMSimulationCalibration::PMTSPESigma)
        .def_readwrite("Rs",&CCMSimulationCalibration::Rs)
        .def_readwrite("Rt",&CCMSimulationCalibration::Rt)
        .def_readwrite("tau_s",&CCMSimulationCalibration::tau_s)
        .def_readwrite("tau_t",&CCMSimulationCalibration::tau_t)
        .def_readwrite("tau_other",&CCMSimulationCalibration::tau_other)
        ;
}
    register_pointer_conversions<CCMSimulationCalibration>();
}
