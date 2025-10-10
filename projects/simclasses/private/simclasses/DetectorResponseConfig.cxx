#include "simclasses/DetectorResponseConfig.h"

#include "icetray/I3FrameObject.h"

#if __has_include("G4Version.hh")
    #define CCM_HAS_G4 1
    #include "G4SystemOfUnits.hh"
#endif

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
        << "\n  Birks Constant :" << birks_constant_
        << "\n  Mean Excitation Energy :" << mean_excitation_energy_
        << "\n  TPB WLS Time Constant :" << tpb_wls_time_constant_
#ifdef CCM_HAS_G4
        << "\n  Use G4 Units :" << use_g4_units_
#endif
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
            rayleigh_scattering_length_128nm_ == rhs.rayleigh_scattering_length_128nm_ &&
            birks_constant_ == rhs.birks_constant_ &&
            mean_excitation_energy_ == rhs.mean_excitation_energy_ &&
            tpb_wls_time_constant_ == rhs.tpb_wls_time_constant_
        );
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig const & bcm) {
    return bcm.Print(oss);
}

std::ostream& operator<<(std::ostream& oss, DetectorResponseConfig & bcm) {
    return bcm.Print(oss);
}

#ifdef CCM_HAS_G4
void DetectorResponseConfig::to_g4_units() {
    if(use_g4_units_) return;
    // Convert units to Geant4 internal units
    rayleigh_scattering_length_ = rayleigh_scattering_length_ / I3Units::cm * CLHEP::cm;
    uv_absorption_a_ = uv_absorption_a_ / (1./I3Units::nanometer) * (1./CLHEP::nanometer);
    uv_absorption_b_ = uv_absorption_b_ / I3Units::nanometer * CLHEP::nanometer;
    uv_absorption_d_ = uv_absorption_d_ / I3Units::cm * CLHEP::cm;
    pmt_tpb_thickness_ = pmt_tpb_thickness_ / I3Units::micrometer * CLHEP::micrometer;
    endcap_tpb_thickness_ = endcap_tpb_thickness_ / I3Units::micrometer * CLHEP::micrometer;
    side_tpb_thickness_ = side_tpb_thickness_ / I3Units::micrometer * CLHEP::micrometer;
    tpb_abs_tau_ = tpb_abs_tau_ / (1./I3Units::nanometer) * (1./CLHEP::nanometer);
    tpb_abs_norm_ = tpb_abs_norm_ / I3Units::nanometer * CLHEP::nanometer;
    mie_scattering_length_200nm_ = mie_scattering_length_200nm_ / I3Units::cm * CLHEP::cm;
    mie_scattering_cutoff_ = mie_scattering_cutoff_ / I3Units::nanometer * CLHEP::nanometer;
    refractive_index_wavelength_UV_ = refractive_index_wavelength_UV_ / I3Units::nanometer * CLHEP::nanometer;
    rayleigh_scattering_length_128nm_ = rayleigh_scattering_length_128nm_ / I3Units::cm * CLHEP::cm;
    birks_constant_ = birks_constant_ / (I3Units::cm/I3Units::MeV) * (CLHEP::cm/CLHEP::MeV);
    mean_excitation_energy_ = mean_excitation_energy_ / I3Units::eV * CLHEP::eV;
    tpb_wls_time_constant_ = tpb_wls_time_constant_ / I3Units::ns * CLHEP::ns;
    use_g4_units_ = true;
}

void DetectorResponseConfig::to_i3_units() {
    if(!use_g4_units_) return;
    // Convert units to IceTray internal units
    rayleigh_scattering_length_ = rayleigh_scattering_length_ / CLHEP::cm * I3Units::cm;
    uv_absorption_a_ = uv_absorption_a_ / (1./CLHEP::nanometer) * (1./I3Units::nanometer);
    uv_absorption_b_ = uv_absorption_b_ / CLHEP::nanometer * I3Units::nanometer;
    uv_absorption_d_ = uv_absorption_d_ / CLHEP::cm * I3Units::cm;
    pmt_tpb_thickness_ = pmt_tpb_thickness_ / CLHEP::micrometer * I3Units::micrometer;
    endcap_tpb_thickness_ = endcap_tpb_thickness_ / CLHEP::micrometer * I3Units::micrometer;
    side_tpb_thickness_ = side_tpb_thickness_ / CLHEP::micrometer * I3Units::micrometer;
    tpb_abs_tau_ = tpb_abs_tau_ / (1./CLHEP::nanometer) * (1./I3Units::nanometer);
    tpb_abs_norm_ = tpb_abs_norm_ / CLHEP::nanometer * I3Units::nanometer;
    mie_scattering_length_200nm_ = mie_scattering_length_200nm_ / CLHEP::cm * I3Units::cm;
    mie_scattering_cutoff_ = mie_scattering_cutoff_ / CLHEP::nanometer * I3Units::nanometer;
    refractive_index_wavelength_UV_ = refractive_index_wavelength_UV_ / CLHEP::nanometer * I3Units::nanometer;
    rayleigh_scattering_length_128nm_ = rayleigh_scattering_length_128nm_ / CLHEP::cm * I3Units::cm;
    birks_constant_ = birks_constant_ / (CLHEP::cm/CLHEP::MeV) * (I3Units::cm/I3Units::MeV);
    mean_excitation_energy_ = mean_excitation_energy_ / CLHEP::eV * I3Units::eV;
    tpb_wls_time_constant_ = tpb_wls_time_constant_ / CLHEP::ns * I3Units::ns;
    use_g4_units_ = false;
}

#endif

I3_SERIALIZABLE(DetectorResponseConfig);
I3_EXPORT_EXTRA_KEY(DetectorResponseConfig, DetectorResponseConfig)

#undef CCM_HAS_G4
