#include <ios>
#include <ostream>
#include "simclasses/CCMSimulationSettings.h"

#include "icetray/I3FrameObject.h"

#if __has_include("G4Version.hh")
    #define CCM_HAS_G4 1
    #include "G4SystemOfUnits.hh"
#endif

std::ostream & CCMSimulationSettings::Print(std::ostream& oss) const {
#ifdef CCM_HAS_G4
    CCMSimulationSettings x = *this;
    if(use_g4_units_) {
        x.to_i3_units();
    }
#endif
    oss << "CCMSimulationSettings: "
        << std::boolalpha
        << "\n  Save All Energy Losses Tree: " << x.save_all_energy_losses_tree_
        << "\n  Veto SD Save Energy Losses Vector: " << x.veto_sd_save_energy_losses_vector_
        << "\n  Veto SD Save Energy Losses Tree: " << x.veto_sd_save_energy_losses_tree_
        << "\n  Veto SD Prune Tree: " << x.veto_sd_prune_tree_
        << "\n  Interior SD Save Energy Losses Vector: " << x.interior_sd_save_energy_losses_vector_
        << "\n  Interior SD Save Energy Losses Tree: " << x.interior_sd_save_energy_losses_tree_
        << "\n  Interior SD Prune Tree: " << x.interior_sd_prune_tree_
        << "\n  Kill Neutrinos: " << x.kill_neutrinos_
        << "\n  Kill Photons: " << x.kill_photons_
        << "\n  Kill Scintillation: " << x.kill_scintillation_
        << "\n  Kill Cherenkov: " << x.kill_cherenkov_
        << "\n  Do Time Cut: " << x.do_time_cut_
        << "\n  Time Cut: " << x.time_cut_ << " ns"
        << "\n  Detailed Photon Tracking: " << x.detailed_photon_tracking_
        << "\n  Track Particles: " << x.track_particles_
        << "\n  Track Energy Losses: " << x.track_energy_losses_
        << "\n  Simulate Nuclear Recoils: " << x.simulate_nuclear_recoils_
        << "\n  G4 Range Cut: " << x.g4_range_cut_ / I3Units::mm << " mm"
        << "\n  G4 Edep Min: " << x.g4_edep_min_ / I3Units::keV << " keV"
        << "\n  G4 E Tracking Min: " << x.g4_e_tracking_min_ / I3Units::keV << " keV"
        << "\n  Record Hits: " << x.record_hits_
        << "\n  Reset Time for Radioactivation: " << x.reset_time_for_radioactivation_
        << "\n  Source Rod In: " << x.source_rod_in_
        << "\n  Source Rod Location: " << x.source_rod_location_ / I3Units::cm << " cm"
        << "\n  Cobalt Source Run: " << x.cobalt_source_run_
        << "\n  Sodium Source Run: " << x.sodium_source_run_
        << "\n  Training Source: " << x.training_source_
        << "\n  Decay X: " << x.decay_x_ / I3Units::cm << " cm"
        << "\n  Decay Y: " << x.decay_y_ / I3Units::cm << " cm"
        << "\n  Decay Z: " << x.decay_z_ / I3Units::cm << " cm"
        << "\n  Random Seed: " << x.random_seed_
#ifdef CCM_HAS_G4
        << "\n  Use G4 Units: " << use_g4_units_
#else
        << "\n  Use G4 Units: False (not compiled with Geant4)"
#endif
        << std::noboolalpha
        << std::endl;
    return oss;
}

bool CCMSimulationSettings::operator==(CCMSimulationSettings const & r) const {
#ifdef CCM_HAS_G4
    CCMSimulationSettings lhs = *this;
    if(use_g4_units_) {
        lhs.to_i3_units();
    }
    CCMSimulationSettings rhs = r;
    if(rhs.use_g4_units_) {
        rhs.to_i3_units();
    }
#else
    CCMSimulationSettings const & lhs = *this;
    CCMSimulationSettings const & rhs = r;
#endif
    return (
        lhs.save_all_energy_losses_tree_ == rhs.save_all_energy_losses_tree_ &&
        lhs.veto_sd_save_energy_losses_vector_ == rhs.veto_sd_save_energy_losses_vector_ &&
        lhs.veto_sd_save_energy_losses_tree_ == rhs.veto_sd_save_energy_losses_tree_ &&
        lhs.veto_sd_prune_tree_ == rhs.veto_sd_prune_tree_ &&
        lhs.interior_sd_save_energy_losses_vector_ == rhs.interior_sd_save_energy_losses_vector_ &&
        lhs.interior_sd_save_energy_losses_tree_ == rhs.interior_sd_save_energy_losses_tree_ &&
        lhs.interior_sd_prune_tree_ == rhs.interior_sd_prune_tree_ &&
        lhs.kill_neutrinos_ == rhs.kill_neutrinos_ &&
        lhs.kill_photons_ == rhs.kill_photons_ &&
        lhs.kill_scintillation_ == rhs.kill_scintillation_ &&
        lhs.kill_cherenkov_ == rhs.kill_cherenkov_ &&
        lhs.do_time_cut_ == rhs.do_time_cut_ &&
        lhs.time_cut_ == rhs.time_cut_ &&
        lhs.detailed_photon_tracking_ == rhs.detailed_photon_tracking_ &&
        lhs.track_particles_ == rhs.track_particles_ &&
        lhs.track_energy_losses_ == rhs.track_energy_losses_ &&
        lhs.simulate_nuclear_recoils_ == rhs.simulate_nuclear_recoils_ &&
        lhs.g4_range_cut_ == rhs.g4_range_cut_ &&
        lhs.g4_edep_min_ == rhs.g4_edep_min_ &&
        lhs.g4_e_tracking_min_ == rhs.g4_e_tracking_min_ &&
        lhs.record_hits_ == rhs.record_hits_ &&
        lhs.reset_time_for_radioactivation_ == rhs.reset_time_for_radioactivation_ &&
        lhs.source_rod_in_ == rhs.source_rod_in_ &&
        lhs.source_rod_location_ == rhs.source_rod_location_ &&
        lhs.cobalt_source_run_ == rhs.cobalt_source_run_ &&
        lhs.sodium_source_run_ == rhs.sodium_source_run_ &&
        lhs.training_source_ == rhs.training_source_ &&
        lhs.decay_x_ == rhs.decay_x_ &&
        lhs.decay_y_ == rhs.decay_y_ &&
        lhs.decay_z_ == rhs.decay_z_ &&
        lhs.random_seed_ == rhs.random_seed_
    );
}

std::ostream& operator<<(std::ostream& oss, CCMSimulationSettings const & bcm) {
    return bcm.Print(oss);
}

std::ostream& operator<<(std::ostream& oss, CCMSimulationSettings & bcm) {
    return bcm.Print(oss);
}

#ifdef CCM_HAS_G4
void CCMSimulationSettings::to_g4_units() {
    if(use_g4_units_) return;
    // Convert units to Geant4 internal units
    g4_range_cut_ = g4_range_cut_ / I3Units::mm * CLHEP::mm;
    g4_edep_min_ = g4_edep_min_ / I3Units::MeV * CLHEP::MeV;
    g4_e_tracking_min_ = g4_e_tracking_min_ / I3Units::MeV * CLHEP::MeV;
    time_cut_ = time_cut_ / I3Units::ns * CLHEP::ns;
    source_rod_location_ = source_rod_location_ / I3Units::cm * CLHEP::cm;
    decay_x_ = decay_x_ / I3Units::cm * CLHEP::cm;
    decay_y_ = decay_y_ / I3Units::cm * CLHEP::cm;
    decay_z_ = decay_z_ / I3Units::cm * CLHEP::cm;
    use_g4_units_ = true;
}

void CCMSimulationSettings::to_i3_units() {
    if(!use_g4_units_) return;
    // Convert units to IceTray internal units
    g4_range_cut_ = g4_range_cut_ / CLHEP::mm * I3Units::mm;
    g4_edep_min_ = g4_edep_min_ / CLHEP::MeV * I3Units::MeV;
    g4_e_tracking_min_ = g4_e_tracking_min_ / CLHEP::MeV * I3Units::MeV;
    time_cut_ = time_cut_ / CLHEP::ns * I3Units::ns;
    source_rod_location_ = source_rod_location_ / CLHEP::cm * I3Units::cm;
    decay_x_ = decay_x_ / CLHEP::cm * I3Units::cm;
    decay_y_ = decay_y_ / CLHEP::cm * I3Units::cm;
    decay_z_ = decay_z_ / CLHEP::cm * I3Units::cm;
    use_g4_units_ = false;
}
#endif

I3_SERIALIZABLE(CCMSimulationSettings);
I3_EXPORT_EXTRA_KEY(CCMSimulationSettings, CCMSimulationSettings)

#undef CCM_HAS_G4
