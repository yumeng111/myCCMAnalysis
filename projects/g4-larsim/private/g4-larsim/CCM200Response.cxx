// standard library stuff

#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"

#include "g4-larsim/CCM200Response.h"
#include "g4-larsim/CCMDetectorResponse.h"
#include "g4-larsim/g4classes/G4Interface.h"

#include "icetray/I3Units.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/IcetrayFwd.h"
#include "icetray/I3FrameObject.h"
#include "icetray/I3ServiceBase.h"
#include "icetray/I3SingleServiceFactory.h"

#include "phys-services/I3RandomService.h"
#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>

CCM200Response::CCM200Response(const I3Context& context) :
    CCMDetectorResponse(context) {

    // Defaults for simulation settings and detector config
    defaultSimulationSettings_ = CCMSimulationSettings();
    defaultDetectorConfig_ = DetectorResponseConfig();
    CCMSimulationSettings & s = simulationSettings_;
    DetectorResponseConfig & c = detectorConfig_;

    AddParameter("SaveAllEnergyLossesTree", "save all energy losses to tree", s.save_all_energy_losses_tree_);
    AddParameter("VetoSDSaveEnergyLossesVector", "save energy losses in veto sensitive detector to vector", s.veto_sd_save_energy_losses_vector_);
    AddParameter("VetoSDSaveEnergyLossesTree", "save energy losses in veto sensitive detector to tree", s.veto_sd_save_energy_losses_tree_);
    AddParameter("VetoSDPruneTree", "prune tree in veto sensitive detector", s.veto_sd_prune_tree_);
    AddParameter("InteriorSDSaveEnergyLossesVector", "save energy losses in interior sensitive detector to vector", s.interior_sd_save_energy_losses_vector_);
    AddParameter("InteriorSDSaveEnergyLossesTree", "save energy losses in interior sensitive detector to tree", s.interior_sd_save_energy_losses_tree_);
    AddParameter("InteriorSDPruneTree", "prune tree in interior sensitive detector", s.interior_sd_prune_tree_);

    AddParameter("KillNeutrinos", "kill neutrinos", s.kill_neutrinos_);
    AddParameter("KillPhotons", "kill photons", s.kill_photons_);
    AddParameter("KillScintillation", "kill scintillation", s.kill_scintillation_);
    AddParameter("KillCherenkov", "kill cherenkov", s.kill_cherenkov_);
    AddParameter("TimeCut", "only track events up to N nsec", s.do_time_cut_);
    AddParameter("TimeCutValue", "only track events up to this time if \"TimeCut\" is set to True [ns]", s.time_cut_ / I3Units::ns);
    AddParameter("DetailedPhotonTracking", "track all optical photons to get distance travelled", s.detailed_photon_tracking_);
    AddParameter("TrackParticles", "track particles", s.track_particles_);
    AddParameter("TrackEnergyLosses", "track energy losses", s.track_energy_losses_);

    AddParameter("SimulateNuclearRecoils", "simulate nuclear recoils", s.simulate_nuclear_recoils_);
    AddParameter("G4RangeCut", "production cut for secondary particles (I3Units) [m]", s.g4_range_cut_);
    AddParameter("G4EDepMin", "minimum energy deposit for secondary particles (I3Units) [GeV]", s.g4_edep_min_);
    AddParameter("G4ETrackingMin", "minimum energy for tracking (I3Units) [GeV]", s.g4_e_tracking_min_);

    AddParameter("RecordHits", "record hits", s.record_hits_);
    AddParameter("ResetTimeForRadioActivation", "If true, then reset the global geant4 time to zero when a radioactive decay occurs", s.reset_time_for_radioactivation_);
    AddParameter("SourceRodIn", "true if we want to simulate the sodium source rod", s.source_rod_in_);
    AddParameter("SourceRodLocation", "z location of the end of the sodium source rod (I3Units) [m]", s.source_rod_location_);
    AddParameter("CobaltSourceRun", "true if we want to simulate cobalt source pellet", s.cobalt_source_run_);
    AddParameter("SodiumSourceRun", "true if we want to simulate sodium source pellet", s.sodium_source_run_);
    AddParameter("TrainingSource", "true if we want to simulate training source events", s.training_source_);
    AddParameter("DecayX", "if generating training source data, provid X position (I3Units) [m]", s.decay_x_);
    AddParameter("DecayY", "if generating training source data, provid Y position (I3Units) [m]", s.decay_y_);
    AddParameter("DecayZ", "if generating training source data, provid Z position (I3Units) [m]", s.decay_z_);
    AddParameter("Rayleigh128Length", "Rayleigh scattering length for 128nm light (I3Units) [m]", c.rayleigh_scattering_length_);
    AddParameter("EnableUVAbsorption", "Enable UV absorption", c.enable_uv_absorption_);
    AddParameter("UVAbsA", "Set UV absorption slope [1/nm]", c.uv_absorption_a_ / (1.0/I3Units::nanometer));
    AddParameter("UVAbsB", "Set UV absorption offset [nm]", c.uv_absorption_b_ / I3Units::nanometer);
    AddParameter("UVAbsD", "Set UV absorption reference distance [cm]", c.uv_absorption_d_ / I3Units::cm);
    AddParameter("UVAbsScaling", "Set UV absorption scale [dimensionless]", c.uv_absorption_scaling_);
    AddParameter("WLSNPhotonsEndCapFoil", "mean number of photons produced per WLS for TPB foils on the end caps of the detector", c.endcap_tpb_qe_);
    AddParameter("WLSNPhotonsSideFoil", "mean number of photons produced per WLS for TPB foils on the sides of the detector", c.side_tpb_qe_);
    AddParameter("WLSNPhotonsPMT", "mean number of photons produced per WLS for TPB on PMTs", c.pmt_tpb_qe_);
    AddParameter("EndCapFoilTPBThickness", "thickness of TPB on the endcap foil (I3Units) [m]", c.endcap_tpb_thickness_);
    AddParameter("SideFoilTPBThickness", "thickness of TPB on the side foil (I3Units) [m]", c.side_tpb_thickness_);
    AddParameter("PMTTPBThickness", "thickness of TPB on the PMTs (I3Units) [m]", c.pmt_tpb_thickness_);
    AddParameter("TPBAbsorptionTau", "factor in exponential for visible part of tpb absorption spectrum [1/nm]", c.tpb_abs_tau_ / (1.0/I3Units::nanometer));
    AddParameter("TPBAbsorptionNorm", "normalization for exponential for visible part of tpb absorption spectrum [nm]", c.tpb_abs_norm_ / I3Units::nanometer);
    AddParameter("TPBAbsorptionScale", "overall scaling for entire tpb absorption curve [dimensionless]", c.tpb_abs_scale_);
    AddParameter("MieGG", "used to calculate angular spread of outgoing photons from mie scattering in tpb", c.mie_gg_);
    AddParameter("MieRatio", "ratio of forward : backwards scattering from mie scattering in tpb", c.mie_ratio_);
    AddParameter("Normalization", "normalization of number of photons produced in liquid argon", c.normalization_);
    AddParameter("PhotonSampling", "scaling of number of photons produced in liquid argon", c.photon_sampling_factor_);
    AddParameter("RindexGamma", "gamma for index of refraction fit (changes cherenkov light production)", c.refractive_index_gamma_UV_);
    AddParameter("BirksConstant", "Birks constant for scintillation quenching [cm/MeV]", c.birks_constant_ / (I3Units::cm/I3Units::MeV));
    AddParameter("MeanExcitationEnergy", "Mean excitation energy of liquid argon [eV]", c.mean_excitation_energy_ / I3Units::eV);
    AddParameter("TPBWLSTimeConstant", "Time constant for re-emission of wavelength-shifted photons in TPB [ns]", c.tpb_wls_time_constant_ / I3Units::ns);
    AddParameter("RandomSeed", "seed for geant4 random generator", s.random_seed_);
}

void CCM200Response::Configure() {
    CCMSimulationSettings & s = simulationSettings_;
    DetectorResponseConfig & c = detectorConfig_;
    GetParameter("SaveAllEnergyLossesTree", s.save_all_energy_losses_tree_);
    GetParameter("VetoSDSaveEnergyLossesVector", s.veto_sd_save_energy_losses_vector_);
    GetParameter("VetoSDSaveEnergyLossesTree", s.veto_sd_save_energy_losses_tree_);
    GetParameter("VetoSDPruneTree", s.veto_sd_prune_tree_);
    GetParameter("InteriorSDSaveEnergyLossesVector", s.interior_sd_save_energy_losses_vector_);
    GetParameter("InteriorSDSaveEnergyLossesTree", s.interior_sd_save_energy_losses_tree_);
    GetParameter("InteriorSDPruneTree", s.interior_sd_prune_tree_);

    GetParameter("KillNeutrinos", s.kill_neutrinos_);
    GetParameter("KillPhotons", s.kill_photons_);
    GetParameter("KillScintillation", s.kill_scintillation_);
    GetParameter("KillCherenkov", s.kill_cherenkov_);
    GetParameter("TimeCut", s.do_time_cut_);
    GetParameter("TimeCutValue", s.time_cut_);
    GetParameter("DetailedPhotonTracking", s.detailed_photon_tracking_);
    GetParameter("TrackParticles", s.track_particles_);
    GetParameter("TrackEnergyLosses", s.track_energy_losses_);

    GetParameter("SimulateNuclearRecoils", s.simulate_nuclear_recoils_);
    GetParameter("G4RangeCut", s.g4_range_cut_);
    GetParameter("G4EDepMin", s.g4_edep_min_);
    GetParameter("G4ETrackingMin", s.g4_e_tracking_min_);

    GetParameter("RecordHits", s.record_hits_);
    GetParameter("ResetTimeForRadioActivation", s.reset_time_for_radioactivation_);
    GetParameter("SourceRodIn", s.source_rod_in_);
    GetParameter("SourceRodLocation", s.source_rod_location_);
    GetParameter("CobaltSourceRun", s.cobalt_source_run_);
    GetParameter("SodiumSourceRun", s.sodium_source_run_);
    GetParameter("TrainingSource", s.training_source_);
    GetParameter("DecayX", s.decay_x_);
    GetParameter("DecayY", s.decay_y_);
    GetParameter("DecayZ", s.decay_z_);
    GetParameter("Rayleigh128Length", c.rayleigh_scattering_length_);
    GetParameter("EnableUVAbsorption", c.enable_uv_absorption_);
    GetParameter("UVAbsA", c.uv_absorption_a_);
    GetParameter("UVAbsB", c.uv_absorption_b_);
    GetParameter("UVAbsD", c.uv_absorption_d_);
    GetParameter("UVAbsScaling", c.uv_absorption_scaling_);
    GetParameter("WLSNPhotonsEndCapFoil", c.endcap_tpb_qe_);
    GetParameter("WLSNPhotonsSideFoil", c.side_tpb_qe_);
    GetParameter("WLSNPhotonsPMT", c.pmt_tpb_qe_);
    GetParameter("EndCapFoilTPBThickness", c.endcap_tpb_thickness_);
    GetParameter("SideFoilTPBThickness", c.side_tpb_thickness_);
    GetParameter("PMTTPBThickness", c.pmt_tpb_thickness_);
    GetParameter("TPBAbsorptionTau", c.tpb_abs_tau_);
    GetParameter("TPBAbsorptionNorm", c.tpb_abs_norm_);
    GetParameter("TPBAbsorptionScale", c.tpb_abs_scale_);
    GetParameter("MieGG", c.mie_gg_);
    GetParameter("MieRatio", c.mie_ratio_);
    GetParameter("Normalization", c.normalization_);
    GetParameter("PhotonSampling", c.photon_sampling_factor_);
    GetParameter("RindexGamma", c.refractive_index_gamma_UV_);
    GetParameter("BirksConstant", c.birks_constant_);
    GetParameter("MeanExcitationEnergy", c.mean_excitation_energy_);
    GetParameter("TPBWLSTimeConstant", c.tpb_wls_time_constant_);
    GetParameter("RandomSeed", s.random_seed_);

    // We read in these parameters in different units than we use them internally
    // So we have to convert them to the I3Units system here
    s.time_cut_ *= I3Units::ns;
    c.uv_absorption_a_ *= (1.0/I3Units::nanometer);
    c.uv_absorption_b_ *= I3Units::nanometer;
    c.uv_absorption_d_ *= I3Units::cm;
    c.tpb_abs_tau_ *= (1.0/I3Units::nanometer);
    c.tpb_abs_norm_ *= I3Units::nanometer;
    c.birks_constant_ *= (I3Units::cm/I3Units::MeV);
    c.mean_excitation_energy_ *= I3Units::eV;
    c.tpb_wls_time_constant_ *= I3Units::ns;

    // Now check if any of the parameters were changed from their default values
    if(!(s == defaultSimulationSettings_)) {
        log_info("Simulation settings changed from default:");
        if(s.save_all_energy_losses_tree_ != defaultSimulationSettings_.save_all_energy_losses_tree_)
            log_info("  SaveAllEnergyLossesTree : %s", s.save_all_energy_losses_tree_ ? "True" : "False");
        if(s.veto_sd_save_energy_losses_vector_ != defaultSimulationSettings_.veto_sd_save_energy_losses_vector_)
            log_info("  VetoSDSaveEnergyLossesVector : %s", s.veto_sd_save_energy_losses_vector_ ? "True" : "False");
        if(s.veto_sd_save_energy_losses_tree_ != defaultSimulationSettings_.veto_sd_save_energy_losses_tree_)
            log_info("  VetoSDSaveEnergyLossesTree : %s", s.veto_sd_save_energy_losses_tree_ ? "True" : "False");
        if(s.veto_sd_prune_tree_ != defaultSimulationSettings_.veto_sd_prune_tree_)
            log_info("  VetoSDPruneTree : %s", s.veto_sd_prune_tree_ ? "True" : "False");
        if(s.interior_sd_save_energy_losses_vector_ != defaultSimulationSettings_.interior_sd_save_energy_losses_vector_)
            log_info("  InteriorSDSaveEnergyLossesVector : %s", s.interior_sd_save_energy_losses_vector_ ? "True" : "False");
        if(s.interior_sd_save_energy_losses_tree_ != defaultSimulationSettings_.interior_sd_save_energy_losses_tree_)
            log_info("  InteriorSDSaveEnergyLossesTree : %s", s.interior_sd_save_energy_losses_tree_ ? "True" : "False");
        if(s.interior_sd_prune_tree_ != defaultSimulationSettings_.interior_sd_prune_tree_)
            log_info("  InteriorSDPruneTree : %s", s.interior_sd_prune_tree_ ? "True" : "False");
        if(s.kill_neutrinos_ != defaultSimulationSettings_.kill_neutrinos_)
            log_info("  KillNeutrinos : %s", s.kill_neutrinos_ ? "True" : "False");
        if(s.kill_photons_ != defaultSimulationSettings_.kill_photons_)
            log_info("  KillPhotons : %s", s.kill_photons_ ? "True" : "False");
        if(s.kill_scintillation_ != defaultSimulationSettings_.kill_scintillation_)
            log_info("  KillScintillation : %s", s.kill_scintillation_ ? "True" : "False");
        if(s.kill_cherenkov_ != defaultSimulationSettings_.kill_cherenkov_)
            log_info("  KillCherenkov : %s", s.kill_cherenkov_ ? "True" : "False");
        if(s.do_time_cut_ != defaultSimulationSettings_.do_time_cut_)
            log_info("  TimeCut : %s", s.do_time_cut_ ? "True" : "False");
        if(s.time_cut_ != defaultSimulationSettings_.time_cut_)
            log_info("  TimeCutValue : %f ns", s.time_cut_/I3Units::ns);
        if(s.detailed_photon_tracking_ != defaultSimulationSettings_.detailed_photon_tracking_)
            log_info("  DetailedPhotonTracking : %s", s.detailed_photon_tracking_ ? "True" : "False");
        if(s.track_particles_ != defaultSimulationSettings_.track_particles_)
            log_info("  TrackParticles : %s", s.track_particles_ ? "True" : "False");
        if(s.track_energy_losses_ != defaultSimulationSettings_.track_energy_losses_)
            log_info("  TrackEnergyLosses : %s", s.track_energy_losses_ ? "True" : "False");
        if(s.simulate_nuclear_recoils_ != defaultSimulationSettings_.simulate_nuclear_recoils_)
            log_info("  SimulateNuclearRecoils : %s", s.simulate_nuclear_recoils_ ? "True" : "False");
        if(s.g4_range_cut_ != defaultSimulationSettings_.g4_range_cut_)
            log_info("  G4RangeCut : %f mm", s.g4_range_cut_/I3Units::mm);
        if(s.g4_edep_min_ != defaultSimulationSettings_.g4_edep_min_)
            log_info("  G4EDepMin : %f keV", s.g4_edep_min_/I3Units::keV);
        if(s.g4_e_tracking_min_ != defaultSimulationSettings_.g4_e_tracking_min_)
            log_info("  G4ETrackingMin : %f keV", s.g4_e_tracking_min_/I3Units::keV);
        if(s.record_hits_ != defaultSimulationSettings_.record_hits_)
            log_info("  RecordHits : %s", s.record_hits_ ? "True" : "False");
        if(s.reset_time_for_radioactivation_ != defaultSimulationSettings_.reset_time_for_radioactivation_)
            log_info("  ResetTimeForRadioactivation : %s", s.reset_time_for_radioactivation_ ? "True" : "False");
        if(s.source_rod_in_ != defaultSimulationSettings_.source_rod_in_)
            log_info("  SourceRodIn : %s", s.source_rod_in_ ? "True" : "False");
        if(s.source_rod_location_ != defaultSimulationSettings_.source_rod_location_)
            log_info("  SourceRodLocation : %f cm", s.source_rod_location_/I3Units::cm);
        if(s.cobalt_source_run_ != defaultSimulationSettings_.cobalt_source_run_)
            log_info("  CobaltSourceRun : %s", s.cobalt_source_run_ ? "True" : "False");
        if(s.sodium_source_run_ != defaultSimulationSettings_.sodium_source_run_)
            log_info("  SodiumSourceRun : %s", s.sodium_source_run_ ? "True" : "False");
        if(s.training_source_ != defaultSimulationSettings_.training_source_)
            log_info("  TrainingSource : %s", s.training_source_ ? "True" : "False");
        if(s.decay_x_ != defaultSimulationSettings_.decay_x_)
            log_info("  DecayX : %f cm", s.decay_x_/I3Units::cm);
        if(s.decay_y_ != defaultSimulationSettings_.decay_y_)
            log_info("  DecayY : %f cm", s.decay_y_/I3Units::cm);
        if(s.decay_z_ != defaultSimulationSettings_.decay_z_)
            log_info("  DecayZ : %f cm", s.decay_z_/I3Units::cm);
        if(s.random_seed_ != defaultSimulationSettings_.random_seed_)
            log_info("  RandomSeed : %lu", s.random_seed_);
    }
    if(!(c == defaultDetectorConfig_)) {
        log_info("Detector configuration changed from default:");
        if(c.rayleigh_scattering_length_ != defaultDetectorConfig_.rayleigh_scattering_length_)
            log_info("  Rayleigh128Length : %f cm", c.rayleigh_scattering_length_/I3Units::cm);
        if(c.enable_uv_absorption_ != defaultDetectorConfig_.enable_uv_absorption_)
            log_info("  EnableUVAbsorption : %s", c.enable_uv_absorption_ ? "True" : "False");
        if(c.uv_absorption_a_ != defaultDetectorConfig_.uv_absorption_a_)
            log_info("  UVAbsA : %f 1/nm", c.uv_absorption_a_ * I3Units::nanometer);
        if(c.uv_absorption_b_ != defaultDetectorConfig_.uv_absorption_b_)
            log_info("  UVAbsB : %f nm", c.uv_absorption_b_ / I3Units::nanometer);
        if(c.uv_absorption_d_ != defaultDetectorConfig_.uv_absorption_d_)
            log_info("  UVAbsD : %f cm", c.uv_absorption_d_ / I3Units::cm);
        if(c.uv_absorption_scaling_ != defaultDetectorConfig_.uv_absorption_scaling_)
            log_info("  UVAbsScaling : %f", c.uv_absorption_scaling_);
        if(c.endcap_tpb_qe_ != defaultDetectorConfig_.endcap_tpb_qe_)
            log_info("  WLSNPhotonsEndCapFoil : %f", c.endcap_tpb_qe_);
        if(c.side_tpb_qe_ != defaultDetectorConfig_.side_tpb_qe_)
            log_info("  WLSNPhotonsSideFoil : %f", c.side_tpb_qe_);
        if(c.pmt_tpb_qe_ != defaultDetectorConfig_.pmt_tpb_qe_)
            log_info("  WLSNPhotonsPMT : %f", c.pmt_tpb_qe_);
        if(c.endcap_tpb_thickness_ != defaultDetectorConfig_.endcap_tpb_thickness_)
            log_info("  EndCapFoilTPBThickness : %f um", c.endcap_tpb_thickness_/I3Units::micrometer);
        if(c.side_tpb_thickness_ != defaultDetectorConfig_.side_tpb_thickness_)
            log_info("  SideFoilTPBThickness : %f um", c.side_tpb_thickness_/I3Units::micrometer);
        if(c.pmt_tpb_thickness_ != defaultDetectorConfig_.pmt_tpb_thickness_)
            log_info("  PMTTPBThickness : %f um", c.pmt_tpb_thickness_/I3Units::micrometer);
        if(c.tpb_abs_tau_ != defaultDetectorConfig_.tpb_abs_tau_)
            log_info("  TPBAbsorptionTau : %f", c.tpb_abs_tau_);
        if(c.tpb_abs_norm_ != defaultDetectorConfig_.tpb_abs_norm_)
            log_info("  TPBAbsorptionNorm : %e", c.tpb_abs_norm_);
        if(c.tpb_abs_scale_ != defaultDetectorConfig_.tpb_abs_scale_)
            log_info("  TPBAbsorptionScale : %f", c.tpb_abs_scale_);
        if(c.mie_gg_ != defaultDetectorConfig_.mie_gg_)
            log_info("  MieGG : %f", c.mie_gg_);
        if(c.mie_ratio_ != defaultDetectorConfig_.mie_ratio_)
            log_info("  MieRatio : %f", c.mie_ratio_);
        if(c.normalization_ != defaultDetectorConfig_.normalization_)
            log_info("  Normalization : %f", c.normalization_);
        if(c.photon_sampling_factor_ != defaultDetectorConfig_.photon_sampling_factor_)
            log_info("  PhotonSampling : %f", c.photon_sampling_factor_);
        if(c.refractive_index_gamma_UV_ != defaultDetectorConfig_.refractive_index_gamma_UV_)
            log_info("  RindexGamma : %f", c.refractive_index_gamma_UV_);
        if(c.refractive_index_a0_ != defaultDetectorConfig_.refractive_index_a0_)
            log_info("  RefractiveIndexA0 : %f", c.refractive_index_a0_);
        if(c.refractive_index_aUV_ != defaultDetectorConfig_.refractive_index_aUV_)
            log_info("  RefractiveIndexAUV : %e", c.refractive_index_aUV_);
        if(c.refractive_index_wavelength_UV_ != defaultDetectorConfig_.refractive_index_wavelength_UV_)
            log_info("  RefractiveIndexWavelengthUV : %f nm", c.refractive_index_wavelength_UV_/I3Units::nanometer);
        if(c.mie_scattering_length_200nm_ != defaultDetectorConfig_.mie_scattering_length_200nm_)
            log_info("  MieScatteringLength200nm : %f cm", c.mie_scattering_length_200nm_/I3Units::cm);
        if(c.mie_scattering_cutoff_ != defaultDetectorConfig_.mie_scattering_cutoff_)
            log_info("  MieScatteringCutoff : %f nm", c.mie_scattering_cutoff_/I3Units::nanometer);
        if(c.rayleigh_scattering_length_128nm_ != defaultDetectorConfig_.rayleigh_scattering_length_128nm_)
            log_info("  RayleighScatteringLength128nm : %f cm", c.rayleigh_scattering_length_128nm_/I3Units::cm);
        if(c.birks_constant_ != defaultDetectorConfig_.birks_constant_)
            log_info("  BirksConstant : %f cm/MeV", c.birks_constant_/(I3Units::cm/I3Units::MeV));
        if(c.mean_excitation_energy_ != defaultDetectorConfig_.mean_excitation_energy_)
            log_info("  MeanExcitationEnergy : %f eV", c.mean_excitation_energy_/I3Units::eV);
        if(c.tpb_wls_time_constant_ != defaultDetectorConfig_.tpb_wls_time_constant_)
            log_info("  TPBWLSTimeConstant : %f ns", c.tpb_wls_time_constant_/I3Units::ns);
    }
}

CCM200Response::~CCM200Response() {}

void CCM200Response::Initialize() {

    if(g4Interface_ == nullptr) {
        g4Interface_ = G4Interface::GetInstance();
    }

    g4Interface_->SetNumberOfThreads(n_threads_);

    // let's let's construct the detector
    g4Interface_->InstallDetector(simulationSettings_, detectorConfig_);
}

I3FrameObjectPtr CCM200Response::GetSimulationConfiguration() {
    DetectorResponseConfigPtr config = boost::make_shared<DetectorResponseConfig>();
    *config = detectorConfig_;
    return config;
}

void CCM200Response::SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries, I3MCTreePtr veto_tree, I3MCTreePtr inner_tree, I3VectorI3ParticlePtr veto_vector, I3VectorI3ParticlePtr inner_vector) {
    g4Interface_->SimulateEvent(primary, tree, mcpeseries, veto_tree, inner_tree, veto_vector, inner_vector);
}

void CCM200Response::SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> veto_trees, std::vector<I3MCTreePtr> inner_trees, std::vector<I3VectorI3ParticlePtr> veto_vectors, std::vector<I3VectorI3ParticlePtr> inner_vectors) {
    g4Interface_->SimulateEvents(primaries, trees, mcpeseries, veto_trees, inner_trees, veto_vectors, inner_vectors);
}

void CCM200Response::DestroyInterface() {

    g4Interface_->DestroyInstance();
}

typedef I3SingleServiceFactory<CCM200Response,CCMDetectorResponse> CCM200ResponseFactory;

I3_SERVICE_FACTORY(CCM200ResponseFactory);


