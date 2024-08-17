#include <vector>
#include <iostream>

#include <dataclasses/physics/AnalyticLightYieldGenerator.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>


namespace bp = boost::python;

std::string to_str(AnalyticLightYieldGenerator const & theSummary) {
    std::ostringstream oss;
    oss << theSummary << std::flush;
    return oss.str();
}


void register_AnalyticLightYieldGenerator() {
    std::string(*convert_to_str)(AnalyticLightYieldGenerator const &) = to_str;
    bp::class_<AnalyticLightYieldGenerator, bp::bases<I3FrameObject>, boost::shared_ptr<AnalyticLightYieldGenerator> >("AnalyticLightYieldGenerator")

    .def_readwrite("Rs",&AnalyticLightYieldGenerator::Rs)
    .def_readwrite("Rt",&AnalyticLightYieldGenerator::Rt)
    .def_readwrite("tau_s",&AnalyticLightYieldGenerator::tau_s)
    .def_readwrite("tau_t",&AnalyticLightYieldGenerator::tau_t)
    .def_readwrite("tau_rec",&AnalyticLightYieldGenerator::tau_rec)
    .def_readwrite("tau_TPB",&AnalyticLightYieldGenerator::tau_TPB)
    .def_readwrite("uv_absorption",&AnalyticLightYieldGenerator::uv_absorption)
    .def_readwrite("normalization",&AnalyticLightYieldGenerator::normalization)
    .def_readwrite("time_offset",&AnalyticLightYieldGenerator::time_offset)
    .def_readwrite("const_offset",&AnalyticLightYieldGenerator::const_offset)
    .def_readwrite("n_sodium_events",&AnalyticLightYieldGenerator::n_sodium_events)
    .def_readwrite("z_offset",&AnalyticLightYieldGenerator::z_offset)
    .def_readwrite("light_profile_type",&AnalyticLightYieldGenerator::light_profile_type)
    .def(bp::dataclass_suite<AnalyticLightYieldGenerator>())
    ;
    boost::python::enum_<AnalyticLightYieldGenerator::LArLightProfileType>("LArLightProfileType")
      .value("Unknown", AnalyticLightYieldGenerator::LArLightProfileType::Unknown)
      .value("Full", AnalyticLightYieldGenerator::LArLightProfileType::Full)
      .value("Simplified", AnalyticLightYieldGenerator::LArLightProfileType::Simplified)
      .export_values()
    ;
    register_pointer_conversions<AnalyticLightYieldGenerator>();

}
