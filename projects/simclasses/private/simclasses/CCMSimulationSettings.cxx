#include <ios>
#include <ostream>
#include "simclasses/CCMSimulationSettings.h"

#include "icetray/I3FrameObject.h"

#if __has_include("G4Version.hh")
    #define CCM_HAS_G4 1
    #include "G4SystemOfUnits.hh"
#endif

std::ostream & CCMSimulationSettings::Print(std::ostream& oss) const {
    oss << "CCMSimulationSettings: "
        //<< "\n  Rayleigh Scattering Length :" << rayleigh_scattering_length_
        << std::boolalpha
        << "\n  Save All Energy Losses Tree :" << save_all_energy_losses_tree_
        << "\n  Veto SD Save Energy Losses Vector :" << veto_sd_save_energy_losses_vector_
        << "\n  Veto SD Save Energy Losses Tree :" << veto_sd_save_energy_losses_tree_
        << "\n  Veto SD Prune Tree :" << veto_sd_prune_tree_
        << "\n  Interior SD Save Energy Losses Vector :" << interior_sd_save_energy_losses_vector_
        << "\n  Interior SD Save Energy Losses Tree :" << interior_sd_save_energy_losses_tree_
        << "\n  Interior SD Prune Tree :" << interior_sd_prune_tree_
        << "\n  Kill Neutrinos :" << kill_neutrinos_
        << "\n  Kill Photons :" << kill_photons_
        << "\n  Kill Scintillation :" << kill_scintillation_
        << "\n  Kill Cherenkov :" << kill_cherenkov_
        << "\n  Do Time Cut :" << do_time_cut_
        << "\n  Time Cut :" << time_cut_
        << "\n  Detailed Photon Tracking :" << detailed_photon_tracking_
        << "\n  Track Particles :" << track_particles_
        << "\n  Track Energy Losses :" << track_energy_losses_
        << "\n  Simulate Nuclear Recoils :" << simulate_nuclear_recoils_
        << "\n  G4 Range Cut :" << g4_range_cut_
        << "\n  G4 Edep Min :" << g4_edep_min_
        << "\n  G4 E Tracking Min :" << g4_e_tracking_min_
        << "\n  Record Hits :" << record_hits_
        << "\n  Source Rod In :" << source_rod_in_
        << "\n  Source Rod Location :" << source_rod_location_
        << "\n  Cobalt Source Run :" << cobalt_source_run_
        << "\n  Sodium Source Run :" << sodium_source_run_
        << "\n  Training Source :" << training_source_
        << "\n  Decay X :" << decay_x_
        << "\n  Decay Y :" << decay_y_
        << "\n  Decay Z :" << decay_z_
        << "\n  Random Seed :" << random_seed_
#ifdef CCM_HAS_G4
        << "\n  Use G4 Units :" << use_g4_units_
#endif
        << std::noboolalpha
        << std::endl;
    return oss;
}

bool CCMSimulationSettings::operator==(const CCMSimulationSettings& rhs) const {
    return (
        save_all_energy_losses_tree_ == rhs.save_all_energy_losses_tree_ &&
        veto_sd_save_energy_losses_vector_ == rhs.veto_sd_save_energy_losses_vector_ &&
        veto_sd_save_energy_losses_tree_ == rhs.veto_sd_save_energy_losses_tree_ &&
        veto_sd_prune_tree_ == rhs.veto_sd_prune_tree_ &&
        interior_sd_save_energy_losses_vector_ == rhs.interior_sd_save_energy_losses_vector_ &&
        interior_sd_save_energy_losses_tree_ == rhs.interior_sd_save_energy_losses_tree_ &&
        interior_sd_prune_tree_ == rhs.interior_sd_prune_tree_ &&
        kill_neutrinos_ == rhs.kill_neutrinos_ &&
        kill_photons_ == rhs.kill_photons_ &&
        kill_scintillation_ == rhs.kill_scintillation_ &&
        kill_cherenkov_ == rhs.kill_cherenkov_ &&
        do_time_cut_ == rhs.do_time_cut_ &&
        time_cut_ == rhs.time_cut_ &&
        detailed_photon_tracking_ == rhs.detailed_photon_tracking_ &&
        track_particles_ == rhs.track_particles_ &&
        track_energy_losses_ == rhs.track_energy_losses_ &&
        simulate_nuclear_recoils_ == rhs.simulate_nuclear_recoils_ &&
        g4_range_cut_ == rhs.g4_range_cut_ &&
        g4_edep_min_ == rhs.g4_edep_min_ &&
        g4_e_tracking_min_ == rhs.g4_e_tracking_min_ &&
        record_hits_ == rhs.record_hits_ &&
        source_rod_in_ == rhs.source_rod_in_ &&
        source_rod_location_ == rhs.source_rod_location_ &&
        cobalt_source_run_ == rhs.cobalt_source_run_ &&
        sodium_source_run_ == rhs.sodium_source_run_ &&
        training_source_ == rhs.training_source_ &&
        decay_x_ == rhs.decay_x_ &&
        decay_y_ == rhs.decay_y_ &&
        decay_z_ == rhs.decay_z_ &&
        random_seed_ == rhs.random_seed_
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
