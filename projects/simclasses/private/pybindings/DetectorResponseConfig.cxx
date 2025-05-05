
//class DetectorResponseConfig : public I3FrameObject {
//public:
//    double rayleigh_scattering_length_;
//    bool enable_uv_absorption_;
//    double uv_absorption_a_;
//    double uv_absorption_b_;
//    double uv_absorption_d_;
//    double uv_absorption_scaling_;
//    double pmt_tpb_qe_;
//    double endcap_tpb_qe_;
//    double side_tpb_qe_;
//    double pmt_tpb_thickness_;
//    double endcap_tpb_thickness_;
//    double side_tpb_thickness_;
//    double tpb_abs_tau_;
//    double tpb_abs_norm_;
//    double tpb_abs_scale_;
//    double mie_gg_;
//    double mie_ratio_;
//    double normalization_;
//    double photon_sampling_factor_;
//
//    DetectorResponseConfig() = default;
//    virtual ~DetectorResponseConfig() override = default;
//
//    std::ostream& Print(std::ostream&) const;
//
//    bool operator==(const DetectorResponseConfig& rhs) const;
//private:
//    friend class icecube::serialization::access;
//    template<class Archive> void save(Archive& ar, unsigned version) const;
//    template<class Archive> void load(Archive& ar, unsigned version);
//    I3_SERIALIZATION_SPLIT_MEMBER();
//};
//
//std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm);
//std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm);

#include <vector>
#include <iostream>

#include <simclasses/DetectorResponseConfig.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

void register_DetectorResponseConfig() {
    class_<DetectorResponseConfig>("DetectorResponseConfig")
        .def(dataclass_suite<DetectorResponseConfig>())
        .def(init<>())
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
        .add_property("__str__", make_function(&DetectorResponseConfig::Print, return_internal_reference<>()))
    ;
}
