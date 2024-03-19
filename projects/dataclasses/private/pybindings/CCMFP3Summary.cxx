#include <vector>
#include <iostream>

#include <dataclasses/physics/CCMFP3Summary.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

namespace bp = boost::python;

std::string to_str(CCMFP3Gamma const & theGamma) {
    std::ostringstream oss;
    oss << theGamma << std::flush;
    return oss.str();
}

std::string to_str(CCMFP3Summary const & theSummary) {
    std::ostringstream oss;
    oss << theSummary << std::flush;
    return oss.str();
}

void register_CCMFP3Gamma() {

    std::string(*convert_to_str)(CCMFP3Gamma const &) = to_str;
    bp::class_<CCMFP3Gamma, boost::shared_ptr<CCMFP3Gamma> >("CCMFP3Gamma")

    .def_readwrite("start_time", &CCMFP3Gamma::gamma_start_time)
    .def_readwrite("end_time", &CCMFP3Gamma::gamma_end_time)
    .def_readwrite("peak_time", &CCMFP3Gamma::gamma_peak_time)
    .def_readwrite("peak_value", &CCMFP3Gamma::gamma_peak_value)
    .def_readwrite("integral", &CCMFP3Gamma::gamma_integral)
    .def_readwrite("derivative", &CCMFP3Gamma::gamma_derivative)
    .def_readwrite("second_derivative", &CCMFP3Gamma::gamma_second_derivative)
    .def_readwrite("local_average", &CCMFP3Gamma::gamma_local_average)
    .def(bp::dataclass_suite<CCMFP3Gamma>())
    ;
}

void register_CCMFP3Summary() {

    std::string(*convert_to_str)(CCMFP3Summary const &) = to_str;
    bp::class_<CCMFP3Summary, bp::bases<I3FrameObject>, boost::shared_ptr<CCMFP3Summary> >("CCMFP3Summary")

    .def_readwrite("waveform_length", &CCMFP3Summary::fp3_waveform_length)
    .def_readwrite("baseline", &CCMFP3Summary::fp3_baseline)
    .def_readwrite("baseline_stddev", &CCMFP3Summary::fp3_baseline_stddev)
    .def_readwrite("num_noise_peaks", &CCMFP3Summary::fp3_num_noise_peaks)
    .def_readwrite("average_noise_level", &CCMFP3Summary::fp3_average_noise_level)
    .def_readwrite("max_noise_level", &CCMFP3Summary::fp3_max_noise_level)
    .def_readwrite("neutron_start_time", &CCMFP3Summary::fp3_neutron_start_time)
    .def_readwrite("neutron_end_time", &CCMFP3Summary::fp3_neutron_end_time)
    .def_readwrite("neutron_derivative", &CCMFP3Summary::fp3_neutron_derivative)
    .def_readwrite("neutron_second_derivative", &CCMFP3Summary::fp3_neutron_second_derivative)
    .def_readwrite("neutron_local_average", &CCMFP3Summary::fp3_neutron_local_average)
    .def_readwrite("neutron_peak_time", &CCMFP3Summary::fp3_neutron_peak_time)
    .def_readwrite("neutron_peak_value", &CCMFP3Summary::fp3_neutron_peak_value)
    .def_readwrite("neutron_integral", &CCMFP3Summary::fp3_neutron_integral)
    .def_readwrite("gammas", &CCMFP3Summary::fp3_gammas)
    .def(bp::dataclass_suite<CCMFP3Summary>())
    ;

    register_pointer_conversions<CCMFP3Summary>();
}

