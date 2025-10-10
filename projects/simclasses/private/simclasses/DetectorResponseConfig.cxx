#include "simclasses/DetectorResponseConfig.h"

#include "icetray/I3FrameObject.h"

std::ostream & DetectorResponseConfig::Print(std::ostream& oss) const {
    oss << "DetectorResponseConfig: "
        << "\n  Rayleigh Scattering Length :" << rayleigh_scattering_length_
        << "\n  Enable UV Absorption :" << enable_uv_absorption_
        << "\n  UV Absorption A :" << uv_absorption_a_
        << "\n  UV Absorption B :" << uv_absorption_b_
        << "\n  UV Absorption D :" << uv_absorption_d_
        << "\n  UV Absorption Scaling :" << uv_absorption_scaling_
        << "\n  PMT TPB QE :" << pmt_tpb_qe_
        << "\n  Endcap TPB QE :" << endcap_tpb_qe_
        << "\n  Side TPB QE :" << side_tpb_qe_
        << "\n  PMT TPB Thickness :" << pmt_tpb_thickness_
        << "\n  Endcap TPB Thickness :" << endcap_tpb_thickness_
        << "\n  Side TPB Thickness :" << side_tpb_thickness_
        << "\n  TPB Absorption Tau :" << tpb_abs_tau_
        << "\n  TPB Absorption Norm :" << tpb_abs_norm_
        << "\n  TPB Absorption Scale :" << tpb_abs_scale_
        << "\n  Mie GG :" << mie_gg_
        << "\n  Mie Ratio :" << mie_ratio_
        << "\n  Normalization :" << normalization_
        << "\n  Photon Sampling Factor :" << photon_sampling_factor_
        << "\n  Mie Scattering Length at 200nm :" << mie_scattering_length_200nm_
        << "\n  Mie Scattering Cutoff :" << mie_scattering_cutoff_
        << "\n  Refractive Index a0 :" << refractive_index_a0_
        << "\n  Refractive Index aUV :" << refractive_index_aUV_
        << "\n  Refractive Index gamma_UV :" << refractive_index_gamma_UV_
        << "\n  Refractive Index wavelength_UV :" << refractive_index_wavelength_UV_
        << "\n  Rayleigh Scattering Length at 128nm :" << rayleigh_scattering_length_128nm_
        << std::endl;
    return oss;
}

bool DetectorResponseConfig::operator==(const DetectorResponseConfig& rhs) const {
    return (rayleigh_scattering_length_ == rhs.rayleigh_scattering_length_ &&
            enable_uv_absorption_ == rhs.enable_uv_absorption_ &&
            uv_absorption_a_ == rhs.uv_absorption_a_ &&
            uv_absorption_b_ == rhs.uv_absorption_b_ &&
            uv_absorption_d_ == rhs.uv_absorption_d_ &&
            uv_absorption_scaling_ == rhs.uv_absorption_scaling_ &&
            pmt_tpb_qe_ == rhs.pmt_tpb_qe_ &&
            endcap_tpb_qe_ == rhs.endcap_tpb_qe_ &&
            side_tpb_qe_ == rhs.side_tpb_qe_ &&
            pmt_tpb_thickness_ == rhs.pmt_tpb_thickness_ &&
            endcap_tpb_thickness_ == rhs.endcap_tpb_thickness_ &&
            side_tpb_thickness_ == rhs.side_tpb_thickness_ &&
            tpb_abs_tau_ == rhs.tpb_abs_tau_ &&
            tpb_abs_norm_ == rhs.tpb_abs_norm_ &&
            tpb_abs_scale_ == rhs.tpb_abs_scale_ &&
            mie_gg_ == rhs.mie_gg_ &&
            mie_ratio_ == rhs.mie_ratio_ &&
            normalization_ == rhs.normalization_ &&
            photon_sampling_factor_ == rhs.photon_sampling_factor_ &&
            mie_scattering_length_200nm_ == rhs.mie_scattering_length_200nm_ &&
            mie_scattering_cutoff_ == rhs.mie_scattering_cutoff_ &&
            refractive_index_a0_ == rhs.refractive_index_a0_ &&
            refractive_index_aUV_ == rhs.refractive_index_aUV_ &&
            refractive_index_gamma_UV_ == rhs.refractive_index_gamma_UV_ &&
            refractive_index_wavelength_UV_ == rhs.refractive_index_wavelength_UV_ &&
            rayleigh_scattering_length_128nm_ == rhs.rayleigh_scattering_length_128nm_);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm) {
    return bcm.Print(oss);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm) {
    return bcm.Print(oss);
}

I3_SERIALIZABLE(DetectorResponseConfig);
I3_EXPORT_EXTRA_KEY(DetectorResponseConfig, DetectorResponseConfig)
