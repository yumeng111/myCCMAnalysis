#include <ios>
#include <ostream>
#include "simclasses/DetectorResponseConfig.h"

#include "icetray/I3FrameObject.h"

#if __has_include("G4Version.hh")
    #define CCM_HAS_G4 1
    #include "G4SystemOfUnits.hh"
#endif

std::ostream & DetectorResponseConfig::Print(std::ostream& oss) const {
#ifdef CCM_HAS_G4
    DetectorResponseConfig x = *this;
    if(use_g4_units_) {
        x.to_i3_units();
    }
#endif
    oss << "DetectorResponseConfig: "
        << "\n  Rayleigh Scattering Length: " << x.rayleigh_scattering_length_ / I3Units::cm << " cm"
        << "\n  Enable UV Absorption: " << x.enable_uv_absorption_
        << "\n  UV Absorption A: " << x.uv_absorption_a_ / (1.0/I3Units::nanometer) << " 1/nm"
        << "\n  UV Absorption B: " << x.uv_absorption_b_ / I3Units::nanometer << " nm"
        << "\n  UV Absorption D: " << x.uv_absorption_d_ / I3Units::cm << " cm"
        << "\n  UV Absorption Scaling: " << x.uv_absorption_scaling_
        << "\n  PMT TPB QE: " << x.pmt_tpb_qe_
        << "\n  Endcap TPB QE: " << x.endcap_tpb_qe_
        << "\n  Side TPB QE: " << x.side_tpb_qe_
        << "\n  PMT TPB Thickness: " << x.pmt_tpb_thickness_ / I3Units::micrometer << " um"
        << "\n  Endcap TPB Thickness: " << x.endcap_tpb_thickness_ / I3Units::micrometer << " um"
        << "\n  Side TPB Thickness: " << x.side_tpb_thickness_ / I3Units::micrometer << " um"
        << "\n  TPB Absorption Tau: " << x.tpb_abs_tau_ / (1.0 / I3Units::nanometer) << " 1/nm"
        << "\n  TPB Absorption Norm: " << x.tpb_abs_norm_ / I3Units::nanometer << " nm"
        << "\n  TPB Absorption Scale: " << x.tpb_abs_scale_
        << "\n  Mie GG: " << x.mie_gg_
        << "\n  Mie Ratio: " << x.mie_ratio_
        << "\n  Normalization: " << x.normalization_
        << "\n  Photon Sampling Factor: " << x.photon_sampling_factor_
        << "\n  Mie Scattering Length at 200nm: " << x.mie_scattering_length_200nm_ / I3Units::cm << " cm"
        << "\n  Mie Scattering Cutoff: " << x.mie_scattering_cutoff_ / I3Units::nanometer << " nm"
        << "\n  Refractive Index a0: " << x.refractive_index_a0_
        << "\n  Refractive Index aUV: " << x.refractive_index_aUV_
        << "\n  Refractive Index gamma_UV: " << x.refractive_index_gamma_UV_
        << "\n  Refractive Index wavelength_UV: " << x.refractive_index_wavelength_UV_ / I3Units::nanometer << " nm"
        << "\n  Rayleigh Scattering Length at 128nm: " << x.rayleigh_scattering_length_128nm_ / I3Units::cm << " cm"
        << "\n  Birks Constant: " << x.birks_constant_ / (I3Units::cm/I3Units::MeV) << " cm/MeV"
        << "\n  Mean Excitation Energy: " << x.mean_excitation_energy_ / I3Units::eV << " eV"
        << "\n  TPB WLS Time Constant: " << x.tpb_wls_time_constant_ / I3Units::ns << " ns"
#ifdef CCM_HAS_G4
        << std::boolalpha
        << "\n  Use G4 Units: " << use_g4_units_
        << std::noboolalpha
#else
        << "\n  Use G4 Units: False (not compiled with Geant4)"
#endif
        << std::endl;
    return oss;
}

bool DetectorResponseConfig::operator==(DetectorResponseConfig const & r) const {
#if CCM_HAS_G4
    DetectorResponseConfig lhs = *this;
    if(use_g4_units_) {
        lhs.to_i3_units();
    }
    DetectorResponseConfig rhs = r;
    if(rhs.use_g4_units_) {
        rhs.to_i3_units();
    }
#else
    DetectorResponseConfig const & lhs = *this;
    DetectorResponseConfig const & rhs = r;
#endif
    return (lhs.rayleigh_scattering_length_ == rhs.rayleigh_scattering_length_ &&
            lhs.enable_uv_absorption_ == rhs.enable_uv_absorption_ &&
            lhs.uv_absorption_a_ == rhs.uv_absorption_a_ &&
            lhs.uv_absorption_b_ == rhs.uv_absorption_b_ &&
            lhs.uv_absorption_d_ == rhs.uv_absorption_d_ &&
            lhs.uv_absorption_scaling_ == rhs.uv_absorption_scaling_ &&
            lhs.pmt_tpb_qe_ == rhs.pmt_tpb_qe_ &&
            lhs.endcap_tpb_qe_ == rhs.endcap_tpb_qe_ &&
            lhs.side_tpb_qe_ == rhs.side_tpb_qe_ &&
            lhs.pmt_tpb_thickness_ == rhs.pmt_tpb_thickness_ &&
            lhs.endcap_tpb_thickness_ == rhs.endcap_tpb_thickness_ &&
            lhs.side_tpb_thickness_ == rhs.side_tpb_thickness_ &&
            lhs.tpb_abs_tau_ == rhs.tpb_abs_tau_ &&
            lhs.tpb_abs_norm_ == rhs.tpb_abs_norm_ &&
            lhs.tpb_abs_scale_ == rhs.tpb_abs_scale_ &&
            lhs.mie_gg_ == rhs.mie_gg_ &&
            lhs.mie_ratio_ == rhs.mie_ratio_ &&
            lhs.normalization_ == rhs.normalization_ &&
            lhs.photon_sampling_factor_ == rhs.photon_sampling_factor_ &&
            lhs.mie_scattering_length_200nm_ == rhs.mie_scattering_length_200nm_ &&
            lhs.mie_scattering_cutoff_ == rhs.mie_scattering_cutoff_ &&
            lhs.refractive_index_a0_ == rhs.refractive_index_a0_ &&
            lhs.refractive_index_aUV_ == rhs.refractive_index_aUV_ &&
            lhs.refractive_index_gamma_UV_ == rhs.refractive_index_gamma_UV_ &&
            lhs.refractive_index_wavelength_UV_ == rhs.refractive_index_wavelength_UV_ &&
            lhs.rayleigh_scattering_length_128nm_ == rhs.rayleigh_scattering_length_128nm_ &&
            lhs.birks_constant_ == rhs.birks_constant_ &&
            lhs.mean_excitation_energy_ == rhs.mean_excitation_energy_ &&
            lhs.tpb_wls_time_constant_ == rhs.tpb_wls_time_constant_
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
