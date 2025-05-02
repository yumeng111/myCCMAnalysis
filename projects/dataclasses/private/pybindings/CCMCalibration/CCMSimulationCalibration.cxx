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

void register_CCMPulseTimeDistributionParameters() {
    class_<CCMPulseTimeDistributionParameters, boost::shared_ptr<CCMPulseTimeDistributionParameters> >("CCMPulseTimeDistributionParameters")
        .def(copy_suite<CCMPulseTimeDistributionParameters>())
        .def(dataclass_suite<CCMPulseTimeDistributionParameters>())
        .def(init<double, double, double>())
        .def(init<std::array<double, 3> >())
        .def(init<std::initializer_list<double> >())
        .def_readwrite("mu", &CCMPulseTimeDistributionParameters::mu)
        .def_readwrite("sigma", &CCMPulseTimeDistributionParameters::sigma)
        .def_readwrite("fraction", &CCMPulseTimeDistributionParameters::fraction)
        .def("__eq__", &CCMPulseTimeDistributionParameters::operator==)
        .def("__str__", to_str<CCMPulseTimeDistributionParameters>)
        ;
    class_<CCMPulseTimeDistributionParametersSeries, boost::shared_ptr<CCMPulseTimeDistributionParametersSeries>>("CCMPulseTimeDistributionParametersSeries")
        .def(dataclass_suite<CCMPulseTimeDistributionParametersSeries>())
        ;
}

void register_CCMSimulationPMTCalibration() {
    class_<CCMSimulationPMTCalibration, boost::shared_ptr<CCMSimulationPMTCalibration> >("CCMSimulationPMTCalibration")
        .def(copy_suite<CCMSimulationPMTCalibration>())
        .def(dataclass_suite<CCMSimulationPMTCalibration>())
        .def_readwrite("pmt_efficiency",&CCMSimulationPMTCalibration::pmt_efficiency)
        .def_readwrite("pmt_spe_mu",&CCMSimulationPMTCalibration::pmt_spe_mu)
        .def_readwrite("pmt_spe_sigma",&CCMSimulationPMTCalibration::pmt_spe_sigma)
        .def_readwrite("pmt_spe_threshold",&CCMSimulationPMTCalibration::pmt_spe_threshold)
        .def_readwrite("main_pulse_mu",&CCMSimulationPMTCalibration::main_pulse_mu)
        .def_readwrite("main_pulse_sigma",&CCMSimulationPMTCalibration::main_pulse_sigma)
        .def_readwrite("late_pulses",&CCMSimulationPMTCalibration::late_pulses)
        ;
    class_<CCMSimulationPMTCalibrationMap, boost::shared_ptr<CCMSimulationPMTCalibrationMap>>("CCMSimulationPMTCalibrationMap")
        .def(dataclass_suite<CCMSimulationPMTCalibrationMap>())
        ;
}

void register_CCMSimulationCalibration() {
{
    class_<CCMSimulationCalibration, bases<I3FrameObject>, boost::shared_ptr<CCMSimulationCalibration> >("CCMSimulationCalibration")
        .def(copy_suite<CCMSimulationCalibration>())
        .def(dataclass_suite<CCMSimulationCalibration>())
        .def_readwrite("pmt_calibration",&CCMSimulationCalibration::pmt_calibration)
        .def_readwrite("Rs",&CCMSimulationCalibration::Rs)
        .def_readwrite("Rt",&CCMSimulationCalibration::Rt)
        .def_readwrite("tau_s",&CCMSimulationCalibration::tau_s)
        .def_readwrite("tau_t",&CCMSimulationCalibration::tau_t)
        .def_readwrite("tau_other",&CCMSimulationCalibration::tau_other)
        .def_readwrite("uv_absorption_a",&CCMSimulationCalibration::uv_absorption_a)
        .def_readwrite("uv_absorption_b",&CCMSimulationCalibration::uv_absorption_b)
        .def_readwrite("uv_absorption_d",&CCMSimulationCalibration::uv_absorption_d)
        .def_readwrite("uv_absorption_scaling",&CCMSimulationCalibration::uv_absorption_scaling)
        ;
}
    register_pointer_conversions<CCMSimulationCalibration>();
}
