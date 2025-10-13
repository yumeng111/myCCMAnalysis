#ifndef CCMSimulationSettings_H_INCLUDED
#define CCMSimulationSettings_H_INCLUDED

#include <icetray/I3DefaultName.h>
#include "icetray/I3FrameObject.h"
#include "icetray/I3Units.h"
#include "dataclasses/Utility.h"

#include <string>
#include <iostream>
#include <sstream>

#if __has_include("G4Version.hh")
    #define CCM_HAS_G4 1
#endif

static const unsigned ccm_simulation_settings_version_ = 1;
class CCMSimulationSettings : public I3FrameObject {
public:
    bool save_all_energy_losses_tree_ = false;
    bool veto_sd_save_energy_losses_vector_ = false;
    bool veto_sd_save_energy_losses_tree_ = false;
    bool veto_sd_prune_tree_ = false;
    bool interior_sd_save_energy_losses_vector_ = false;
    bool interior_sd_save_energy_losses_tree_ = false;
    bool interior_sd_prune_tree_ = false;
    bool kill_neutrinos_ = true;
    bool kill_photons_ = false;
    bool kill_scintillation_ = false;
    bool kill_cherenkov_ = false;
    bool do_time_cut_ = false;
    double time_cut_ = 200.0 * I3Units::ns;
    bool detailed_photon_tracking_ = false;
    bool track_particles_ = true;
    bool track_energy_losses_ = false;
    bool simulate_nuclear_recoils_ = false;
    double g4_range_cut_ = 0.7 * I3Units::mm;
    double g4_edep_min_ = 1.0 * I3Units::keV;
    double g4_e_tracking_min_ = 1.0 * I3Units::keV;
    bool record_hits_ = true;

    bool reset_time_for_radioactivation_ = false;
    bool source_rod_in_ = false;
    double source_rod_location_ = 0.0*I3Units::cm;
    bool cobalt_source_run_ = false;
    bool sodium_source_run_ = false;
    bool training_source_ = false;
    double decay_x_ = 0.0*I3Units::cm;
    double decay_y_ = 0.0*I3Units::cm;
    double decay_z_ = 0.0*I3Units::cm;
    long random_seed_ = 0;

#ifdef CCM_HAS_G4
    bool use_g4_units_ = false;
    void to_g4_units();
    void to_i3_units();
#endif

    CCMSimulationSettings() = default;
    virtual ~CCMSimulationSettings() override = default;

    std::ostream& Print(std::ostream&) const override;

    bool operator==(const CCMSimulationSettings& rhs) const;
private:
    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

std::ostream& operator<<(std::ostream& oss, CCMSimulationSettings const & bcm);
std::ostream& operator<<(std::ostream& oss, CCMSimulationSettings & bcm);

template <class Archive>
void CCMSimulationSettings::save(Archive& ar, unsigned version) const {
    if(version > ccm_simulation_settings_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMSimulationSettings class.",version,ccm_simulation_settings_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("save_all_energy_losses_tree", save_all_energy_losses_tree_);
    ar & make_nvp("veto_sd_save_energy_losses_vector", veto_sd_save_energy_losses_vector_);
    ar & make_nvp("veto_sd_save_energy_losses_tree", veto_sd_save_energy_losses_tree_);
    ar & make_nvp("veto_sd_prune_tree", veto_sd_prune_tree_);
    ar & make_nvp("interior_sd_save_energy_losses_vector", interior_sd_save_energy_losses_vector_);
    ar & make_nvp("interior_sd_save_energy_losses_tree", interior_sd_save_energy_losses_tree_);
    ar & make_nvp("interior_sd_prune_tree", interior_sd_prune_tree_);
    ar & make_nvp("kill_neutrinos", kill_neutrinos_);
    ar & make_nvp("kill_photons", kill_photons_);
    ar & make_nvp("kill_scintillation", kill_scintillation_);
    ar & make_nvp("kill_cherenkov", kill_cherenkov_);
    ar & make_nvp("do_time_cut", do_time_cut_);
    ar & make_nvp("time_cut", time_cut_);
    ar & make_nvp("detailed_photon_tracking", detailed_photon_tracking_);
    ar & make_nvp("track_particles", track_particles_);
    ar & make_nvp("track_energy_losses", track_energy_losses_);
    ar & make_nvp("simulate_nuclear_recoils", simulate_nuclear_recoils_);
    ar & make_nvp("g4_range_cut", g4_range_cut_);
    ar & make_nvp("g4_edep_min", g4_edep_min_);
    ar & make_nvp("g4_e_tracking_min", g4_e_tracking_min_);
    ar & make_nvp("record_hits", record_hits_);
    ar & make_nvp("reset_time_for_radioactivation", reset_time_for_radioactivation_);
    ar & make_nvp("source_rod_in", source_rod_in_);
    ar & make_nvp("source_rod_location", source_rod_location_);
    ar & make_nvp("cobalt_source_run", cobalt_source_run_);
    ar & make_nvp("sodium_source_run", sodium_source_run_);
    ar & make_nvp("training_source", training_source_);
    ar & make_nvp("decay_x", decay_x_);
    ar & make_nvp("decay_y", decay_y_);
    ar & make_nvp("decay_z", decay_z_);
    ar & make_nvp("random_seed", random_seed_);
}

template <class Archive>
void CCMSimulationSettings::load(Archive& ar, unsigned version) {
    if(version > ccm_simulation_settings_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMSimulationSettings class.",version,ccm_simulation_settings_version_);
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("save_all_energy_losses_tree", save_all_energy_losses_tree_);
    ar & make_nvp("veto_sd_save_energy_losses_vector", veto_sd_save_energy_losses_vector_);
    ar & make_nvp("veto_sd_save_energy_losses_tree", veto_sd_save_energy_losses_tree_);
    ar & make_nvp("veto_sd_prune_tree", veto_sd_prune_tree_);
    ar & make_nvp("interior_sd_save_energy_losses_vector", interior_sd_save_energy_losses_vector_);
    ar & make_nvp("interior_sd_save_energy_losses_tree", interior_sd_save_energy_losses_tree_);
    ar & make_nvp("interior_sd_prune_tree", interior_sd_prune_tree_);
    ar & make_nvp("kill_neutrinos", kill_neutrinos_);
    ar & make_nvp("kill_photons", kill_photons_);
    ar & make_nvp("kill_scintillation", kill_scintillation_);
    ar & make_nvp("kill_cherenkov", kill_cherenkov_);
    ar & make_nvp("do_time_cut", do_time_cut_);
    ar & make_nvp("time_cut", time_cut_);
    ar & make_nvp("detailed_photon_tracking", detailed_photon_tracking_);
    ar & make_nvp("track_particles", track_particles_);
    ar & make_nvp("track_energy_losses", track_energy_losses_);
    ar & make_nvp("simulate_nuclear_recoils", simulate_nuclear_recoils_);
    ar & make_nvp("g4_range_cut", g4_range_cut_);
    ar & make_nvp("g4_edep_min", g4_edep_min_);
    ar & make_nvp("g4_e_tracking_min", g4_e_tracking_min_);
    ar & make_nvp("record_hits", record_hits_);
    if(ccm_simulation_settings_version_ >= 1) {
        ar & make_nvp("reset_time_for_radioactivation", reset_time_for_radioactivation_);
    } else {
        reset_time_for_radioactivation_ = true;
    }
    ar & make_nvp("source_rod_in", source_rod_in_);
    ar & make_nvp("source_rod_location", source_rod_location_);
    ar & make_nvp("cobalt_source_run", cobalt_source_run_);
    ar & make_nvp("sodium_source_run", sodium_source_run_);
    ar & make_nvp("training_source", training_source_);
    ar & make_nvp("decay_x", decay_x_);
    ar & make_nvp("decay_y", decay_y_);
    ar & make_nvp("decay_z", decay_z_);
    ar & make_nvp("random_seed", random_seed_);
}

I3_CLASS_VERSION(CCMSimulationSettings, ccm_simulation_settings_version_);
I3_POINTER_TYPEDEFS(CCMSimulationSettings);
I3_DEFAULT_NAME(CCMSimulationSettings);

#undef CCM_HAS_G4

#endif // CCMSimulationSettings_H_INCLUDED

