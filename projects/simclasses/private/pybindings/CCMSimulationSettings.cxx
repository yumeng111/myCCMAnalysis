#include <vector>
#include <iostream>

#include <simclasses/CCMSimulationSettings.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

void register_CCMSimulationSettings() {
    class_<CCMSimulationSettings>("CCMSimulationSettings")
        .def(dataclass_suite<CCMSimulationSettings>())
        .def(init<>())
        .def(self == self)
        .def_readwrite("save_all_energy_losses_tree", &CCMSimulationSettings::save_all_energy_losses_tree_)
        .def_readwrite("veto_sd_save_energy_losses_vector", &CCMSimulationSettings::veto_sd_save_energy_losses_vector_)
        .def_readwrite("veto_sd_save_energy_losses_tree", &CCMSimulationSettings::veto_sd_save_energy_losses_tree_)
        .def_readwrite("veto_sd_prune_tree", &CCMSimulationSettings::veto_sd_prune_tree_)
        .def_readwrite("interior_sd_save_energy_losses_vector", &CCMSimulationSettings::interior_sd_save_energy_losses_vector_)
        .def_readwrite("interior_sd_save_energy_losses_tree", &CCMSimulationSettings::interior_sd_save_energy_losses_tree_)
        .def_readwrite("interior_sd_prune_tree", &CCMSimulationSettings::interior_sd_prune_tree_)
        .def_readwrite("kill_neutrinos", &CCMSimulationSettings::kill_neutrinos_)
        .def_readwrite("kill_photons", &CCMSimulationSettings::kill_photons_)
        .def_readwrite("kill_scintillation", &CCMSimulationSettings::kill_scintillation_)
        .def_readwrite("kill_cherenkov", &CCMSimulationSettings::kill_cherenkov_)
        .def_readwrite("do_time_cut", &CCMSimulationSettings::do_time_cut_)
        .def_readwrite("time_cut", &CCMSimulationSettings::time_cut_)
        .def_readwrite("detailed_photon_tracking", &CCMSimulationSettings::detailed_photon_tracking_)
        .def_readwrite("track_particles", &CCMSimulationSettings::track_particles_)
        .def_readwrite("track_energy_losses", &CCMSimulationSettings::track_energy_losses_)
        .def_readwrite("simulate_nuclear_recoils", &CCMSimulationSettings::simulate_nuclear_recoils_)
        .def_readwrite("g4_range_cut", &CCMSimulationSettings::g4_range_cut_)
        .def_readwrite("g4_edep_min", &CCMSimulationSettings::g4_edep_min_)
        .def_readwrite("g4_e_tracking_min", &CCMSimulationSettings::g4_e_tracking_min_)
        .def_readwrite("record_hits", &CCMSimulationSettings::record_hits_)
        .def_readwrite("reset_time_for_radioactivation", &CCMSimulationSettings::reset_time_for_radioactivation_)
        .def_readwrite("source_rod_in", &CCMSimulationSettings::source_rod_in_)
        .def_readwrite("source_rod_location", &CCMSimulationSettings::source_rod_location_)
        .def_readwrite("cobalt_source_run", &CCMSimulationSettings::cobalt_source_run_)
        .def_readwrite("sodium_source_run", &CCMSimulationSettings::sodium_source_run_)
        .def_readwrite("training_source", &CCMSimulationSettings::training_source_)
        .def_readwrite("decay_x", &CCMSimulationSettings::decay_x_)
        .def_readwrite("decay_y", &CCMSimulationSettings::decay_y_)
        .def_readwrite("decay_z", &CCMSimulationSettings::decay_z_)
        .def_readwrite("random_seed", &CCMSimulationSettings::random_seed_)
        .add_property("__str__", make_function(&CCMSimulationSettings::Print, return_internal_reference<>()))
    ;
}
