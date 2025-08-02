#include <iostream>

#include "dataclasses/physics/NuclearCascadeStep.h"

#include "icetray/python/dataclass_suite.hpp"
#include "dataclasses/ostream_overloads.hpp"


namespace bp = boost::python;

std::string to_str(NuclearCascadeStep const & theSummary) {
    std::ostringstream oss;
    oss << theSummary << std::flush;
    return oss.str();
}

void register_NuclearCascadeStep() {
    std::string(*convert_to_str)(NuclearCascadeStep const &) = to_str;
    bp::class_<NuclearCascadeStep, bp::bases<I3FrameObject>, boost::shared_ptr<NuclearCascadeStep> >("NuclearCascadeStep")
    .def_readwrite("initial_level_index", &NuclearCascadeStep::initial_level_index)
    .def_readwrite("initial_level_energy", &NuclearCascadeStep::initial_level_energy)
    .def_readwrite("initial_level_spin", &NuclearCascadeStep::initial_level_spin)
    .def_readwrite("initial_level_parity", &NuclearCascadeStep::initial_level_parity)
    .def_readwrite("final_level_index", &NuclearCascadeStep::final_level_index)
    .def_readwrite("final_level_energy", &NuclearCascadeStep::final_level_energy)
    .def_readwrite("final_level_spin", &NuclearCascadeStep::final_level_spin)
    .def_readwrite("final_level_parity", &NuclearCascadeStep::final_level_parity)
    .def_readwrite("gamma_energy", &NuclearCascadeStep::gamma_energy)
    .def_readwrite("T12_ns", &NuclearCascadeStep::T12_ns)
    .def_readwrite("tau_ns", &NuclearCascadeStep::tau_ns)
    .def_readwrite("sampled_delay_ns", &NuclearCascadeStep::sampled_delay_ns)
    .def_readwrite("cumulative_time_ns", &NuclearCascadeStep::cumulative_time_ns)
    .def(bp::dataclass_suite<NuclearCascadeStep>())
    ;
    register_pointer_conversions<NuclearCascadeStep>();
}
