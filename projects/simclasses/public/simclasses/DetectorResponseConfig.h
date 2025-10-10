#ifndef DetectorResponseConfig_H_INCLUDED
#define DetectorResponseConfig_H_INCLUDED

#include <icetray/I3DefaultName.h>
#include "icetray/I3FrameObject.h"
#include "icetray/I3Units.h"
#include "dataclasses/Utility.h"

#include <iostream>

static const unsigned g4ccmconfig_version_ = 8;
class DetectorResponseConfig : public I3FrameObject {
public:
    double rayleigh_scattering_length_;
    bool enable_uv_absorption_;
    double uv_absorption_a_;
    double uv_absorption_b_;
    double uv_absorption_d_;
    double uv_absorption_scaling_;
    double pmt_tpb_qe_;
    double endcap_tpb_qe_;
    double side_tpb_qe_;
    double pmt_tpb_thickness_;
    double endcap_tpb_thickness_;
    double side_tpb_thickness_;
    double tpb_abs_tau_;
    double tpb_abs_norm_;
    double tpb_abs_scale_;

    double mie_gg_;
    double mie_ratio_;
    double normalization_;
    double photon_sampling_factor_;

    // Mie scattering length at 200nm
    double mie_scattering_length_200nm_;
    double mie_scattering_cutoff_;

    // Refractive index parameters in harmonic oscillator model
    double refractive_index_a0_;
    double refractive_index_aUV_;
    double refractive_index_gamma_UV_;
    double refractive_index_wavelength_UV_;

    // Mie scattering parameters
    double rayleigh_scattering_length_128nm_;

    DetectorResponseConfig() = default;
    virtual ~DetectorResponseConfig() override = default;

    std::ostream& Print(std::ostream&) const override;

    bool operator==(const DetectorResponseConfig& rhs) const;
private:
    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm);
std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm);

template <class Archive>
void DetectorResponseConfig::save(Archive& ar, unsigned version) const {
    if(version>g4ccmconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of DetectorResponseConfig class.",version,g4ccmconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("rayleigh_scattering_length", rayleigh_scattering_length_);
    ar & make_nvp("enable_uv_absorption", enable_uv_absorption_);
    ar & make_nvp("uv_absorption_a", uv_absorption_a_);
    ar & make_nvp("uv_absorption_b", uv_absorption_b_);
    ar & make_nvp("uv_absorption_d", uv_absorption_d_);
    ar & make_nvp("uv_absorption_scaling", uv_absorption_scaling_);
    ar & make_nvp("pmt_tpb_qe", pmt_tpb_qe_);
    ar & make_nvp("endcap_tpb_qe", endcap_tpb_qe_);
    ar & make_nvp("side_tpb_qe", side_tpb_qe_);
    ar & make_nvp("pmt_tpb_thickness", pmt_tpb_thickness_);
    ar & make_nvp("endcap_tpb_thickness", endcap_tpb_thickness_);
    ar & make_nvp("side_tpb_thickness", side_tpb_thickness_);
    ar & make_nvp("tpb_abs_tau", tpb_abs_tau_);
    ar & make_nvp("tpb_abs_norm", tpb_abs_norm_);
    ar & make_nvp("tpb_abs_scale", tpb_abs_scale_);
    ar & make_nvp("mie_gg", mie_gg_);
    ar & make_nvp("mie_ratio", mie_ratio_);
    ar & make_nvp("normalization", normalization_);
    ar & make_nvp("photon_sampling_factor", photon_sampling_factor_);
    ar & make_nvp("mie_scattering_length_200nm", mie_scattering_length_200nm_);
    ar & make_nvp("mie_scattering_cutoff", mie_scattering_cutoff_);
    ar & make_nvp("refractive_index_a0", refractive_index_a0_);
    ar & make_nvp("refractive_index_aUV", refractive_index_aUV_);
    ar & make_nvp("refractive_index_gamma_UV", refractive_index_gamma_UV_);
    ar & make_nvp("refractive_index_wavelength_UV", refractive_index_wavelength_UV_);
    ar & make_nvp("rayleigh_scattering_length_128nm", rayleigh_scattering_length_128nm_);
}

template <class Archive>
void DetectorResponseConfig::load(Archive& ar, unsigned version) {
    if(version>g4ccmconfig_version_)
        log_fatal("Attempting to read version %u from file but running version %u of DetectorResponseConfig class.",version,g4ccmconfig_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("rayleigh_scattering_length", rayleigh_scattering_length_);
    if(version >= 7)
        ar & make_nvp("enable_uv_absorption", enable_uv_absorption_);
    else
        enable_uv_absorption_ = false;

    double uv_absorption_length_;
    if(version == 5)
        ar & make_nvp("uv_absorption_length", uv_absorption_length_);
    if(version >= 6) {
        ar & make_nvp("uv_absorption_a", uv_absorption_a_);
        ar & make_nvp("uv_absorption_b", uv_absorption_b_);
        ar & make_nvp("uv_absorption_d", uv_absorption_d_);
        ar & make_nvp("uv_absorption_scaling", uv_absorption_scaling_);
    }
    ar & make_nvp("pmt_tpb_qe", pmt_tpb_qe_);
    ar & make_nvp("endcap_tpb_qe", endcap_tpb_qe_);
    ar & make_nvp("side_tpb_qe", side_tpb_qe_);
    ar & make_nvp("pmt_tpb_thickness", pmt_tpb_thickness_);
    ar & make_nvp("endcap_tpb_thickness", endcap_tpb_thickness_);
    ar & make_nvp("side_tpb_thickness", side_tpb_thickness_);

    if(version == 0) {
        tpb_abs_tau_ = 0.0;
        tpb_abs_norm_ = 0.0;
        tpb_abs_scale_ = 0.0;
        mie_gg_ = 0.0;
        mie_ratio_ = 0.0;
        uv_absorption_length_ = 0.0;
    } else if(version == 1) {
        ar & make_nvp("tpb_abs_tau", tpb_abs_tau_);
        ar & make_nvp("tpb_abs_norm", tpb_abs_norm_);
        tpb_abs_scale_ = 0.0;
        mie_gg_ = 0.0;
        mie_ratio_ = 0.0;
        uv_absorption_length_ = 0.0;
    } else if(version == 2) {
        ar & make_nvp("tpb_abs_tau", tpb_abs_tau_);
        ar & make_nvp("tpb_abs_norm", tpb_abs_norm_);
        ar & make_nvp("tpb_abs_scale", tpb_abs_scale_);
        mie_gg_ = 0.0;
        mie_ratio_ = 0.0;
        uv_absorption_length_ = 0.0;
    }

    if(version >= 3) {
        ar & make_nvp("tpb_abs_tau", tpb_abs_tau_);
        ar & make_nvp("tpb_abs_norm", tpb_abs_norm_);
        ar & make_nvp("tpb_abs_scale", tpb_abs_scale_);
        ar & make_nvp("mie_gg", mie_gg_);
        ar & make_nvp("mie_ratio", mie_ratio_);
    }

    if(version == 3) {
        uv_absorption_length_ = 0.0;
    } else if(version == 4) {
        ar & make_nvp("uv_absorption_length", uv_absorption_length_);
    } else if(version == 5) {
        ar & make_nvp("normalization", normalization_);
    } else {
        ar & make_nvp("normalization", normalization_);
        ar & make_nvp("photon_sampling_factor", photon_sampling_factor_);
    }

    if(version < 6) {
        // Attempt to make something close to a flat absorption length
        uv_absorption_a_ = 1e-8;
        uv_absorption_b_ = 1e8 * log(1.0 - exp(-1.0));
        uv_absorption_d_ = uv_absorption_length_;
        uv_absorption_scaling_ = 1.0;
    }

    if(normalization_ < 0.0 or normalization_ > 10.0 or normalization_ < 1e-8) {
        log_warn("Normalization value %f is out of range [1e-8,10.0]. Resetting to 1.0 since fits of this era have the normalization on the 1.0 boundary", normalization_);
        normalization_ = 1.0;
    }

    if(version < 8) {
        // Set Mie scattering length at 200nm
        mie_scattering_length_200nm_ = 9.37 * I3Units::cm;
        mie_scattering_cutoff_ = 200.0 * (I3Units::m * 1e-9); // nm

        // Set refractive index parameters to some reasonable values
        refractive_index_a0_ = 1.10232;
        refractive_index_aUV_ = 0.00001058;
        refractive_index_gamma_UV_ = 0.0018280705;
        refractive_index_wavelength_UV_ = 106.6 * (I3Units::m * 1e-9); // nm

        // Set Rayleigh scattering length at 128nm
        rayleigh_scattering_length_128nm_ = 99.9769 * I3Units::cm;
    } else {
        ar & make_nvp("mie_scattering_length_200nm", mie_scattering_length_200nm_);
        ar & make_nvp("mie_scattering_cutoff", mie_scattering_cutoff_);
        ar & make_nvp("refractive_index_a0", refractive_index_a0_);
        ar & make_nvp("refractive_index_aUV", refractive_index_aUV_);
        ar & make_nvp("refractive_index_gamma_UV", refractive_index_gamma_UV_);
        ar & make_nvp("refractive_index_wavelength_UV", refractive_index_wavelength_UV_);
        ar & make_nvp("rayleigh_scattering_length_128nm", rayleigh_scattering_length_128nm_);
    }
}

I3_CLASS_VERSION(DetectorResponseConfig, g4ccmconfig_version_);
I3_POINTER_TYPEDEFS(DetectorResponseConfig);
I3_DEFAULT_NAME(DetectorResponseConfig);

#endif // DetectorResponseConfig_H_INCLUDED

