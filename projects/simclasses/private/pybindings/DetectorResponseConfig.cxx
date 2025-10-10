#include <simclasses/DetectorResponseConfig.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

void register_DetectorResponseConfig() {
    class_<DetectorResponseConfig>("DetectorResponseConfig")
        .def(dataclass_suite<DetectorResponseConfig>())
        .def(init<>())
        .def(self == self)
        .def_readwrite("rayleigh_scattering_length", &DetectorResponseConfig::rayleigh_scattering_length_)
        .def_readwrite("enable_uv_absorption", &DetectorResponseConfig::enable_uv_absorption_)
        .def_readwrite("uv_absorption_a", &DetectorResponseConfig::uv_absorption_a_)
        .def_readwrite("uv_absorption_b", &DetectorResponseConfig::uv_absorption_b_)
        .def_readwrite("uv_absorption_d", &DetectorResponseConfig::uv_absorption_d_)
        .def_readwrite("uv_absorption_scaling", &DetectorResponseConfig::uv_absorption_scaling_)
        .def_readwrite("pmt_tpb_qe", &DetectorResponseConfig::pmt_tpb_qe_)
        .def_readwrite("endcap_tpb_qe", &DetectorResponseConfig::endcap_tpb_qe_)
        .def_readwrite("side_tpb_qe", &DetectorResponseConfig::side_tpb_qe_)
        .def_readwrite("pmt_tpb_thickness", &DetectorResponseConfig::pmt_tpb_thickness_)
        .def_readwrite("endcap_tpb_thickness", &DetectorResponseConfig::endcap_tpb_thickness_)
        .def_readwrite("side_tpb_thickness", &DetectorResponseConfig::side_tpb_thickness_)
        .def_readwrite("tpb_abs_tau", &DetectorResponseConfig::tpb_abs_tau_)
        .def_readwrite("tpb_abs_norm", &DetectorResponseConfig::tpb_abs_norm_)
        .def_readwrite("tpb_abs_scale", &DetectorResponseConfig::tpb_abs_scale_)
        .def_readwrite("mie_gg", &DetectorResponseConfig::mie_gg_)
        .def_readwrite("mie_ratio", &DetectorResponseConfig::mie_ratio_)
        .def_readwrite("normalization", &DetectorResponseConfig::normalization_)
        .def_readwrite("photon_sampling_factor", &DetectorResponseConfig::photon_sampling_factor_)
        .def_readwrite("mie_scattering_length_200nm", &DetectorResponseConfig::mie_scattering_length_200nm_)
        .def_readwrite("mie_scattering_cutoff", &DetectorResponseConfig::mie_scattering_cutoff_)
        .def_readwrite("refractive_index_a0", &DetectorResponseConfig::refractive_index_a0_)
        .def_readwrite("refractive_index_aUV", &DetectorResponseConfig::refractive_index_aUV_)
        .def_readwrite("refractive_index_gamma_UV", &DetectorResponseConfig::refractive_index_gamma_UV_)
        .def_readwrite("refractive_index_wavelength_UV", &DetectorResponseConfig::refractive_index_wavelength_UV_)
        .def_readwrite("rayleigh_scattering_length_128nm", &DetectorResponseConfig::rayleigh_scattering_length_128nm_)
        .def_readwrite("birks_constant", &DetectorResponseConfig::birks_constant_)
        .def_readwrite("mean_excitation_energy", &DetectorResponseConfig::mean_excitation_energy_)
        .def_readwrite("tpb_wls_time_constant", &DetectorResponseConfig::tpb_wls_time_constant_)
        .add_property("__str__", make_function(&DetectorResponseConfig::Print, return_internal_reference<>()))
    ;
}
