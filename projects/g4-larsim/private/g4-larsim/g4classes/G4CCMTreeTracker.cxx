
#include "g4-larsim/g4classes/G4CCMTreeTracker.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include <map>
#include <tuple>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include <G4ios.hh>
#include <G4Step.hh>
#include <G4Event.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>
#include <G4SDManager.hh>
#include <G4VTouchable.hh>
#include <G4ParticleTypes.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TouchableHistory.hh>
#include <G4ParticleDefinition.hh>
#include <G4EventManager.hh>
#include <G4TrackingManager.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::unordered_set<int> G4CCMTreeTracker::energyLossPDGCodes = {0, 2000000001, 2000000002, 2000000003, 2000000004, 2000000005, 2000000006, 2000000007, 2000000008, 2000000009, 2000000010, 2000000011, 2000000012, 2000000013, 2000000014, 2000000015, 2000000016, 2000000017, 2000000018};

const std::unordered_set<std::string> G4CCMTreeTracker::knownProcessNames = {"NoProcess", "Transportation"};

const std::unordered_map<std::string, int> G4CCMTreeTracker::energyLossToI3ParticlePDGCode = {{"phot", 2000000001}, {"compt", 2000000002}, {"conv", 2000000003},
                                                                                          {"Rayl", 2000000004}, {"msc", 2000000005}, {"eIoni", 2000000006},
                                                                                          {"eBrem", 2000000007}, {"ePairProd", 2000000008}, {"CoulombScat", 2000000009},
                                                                                          {"annihil", 2000000010}, {"Cerenkov", 2000000011}, {"Radioactivation", 2000000012},
                                                                                          {"Scintillation", 2000000013}, {"OpWLS", 2000000014}, {"ionIoni" , 2000000015},
                                                                                          {"hIoni", 2000000016}, {"neutronInelastic", 2000000017}, {"hadElastic", 2000000018}, {"Unknown", 0}};

const std::unordered_map<std::string, PhotonSummary::PhotonSource> G4CCMTreeTracker::processNameToPhotonSource = {{"Unknown", PhotonSummary::PhotonSource::Unknown},
                                                                                                              {"Scintillation", PhotonSummary::PhotonSource::Scintillation},
                                                                                                              {"Cerenkov", PhotonSummary::PhotonSource::Cherenkov},
                                                                                                              {"OpWLS", PhotonSummary::PhotonSource::OpWLS}};

const std::unordered_map<PhotonSummary::PhotonSource, std::string> G4CCMTreeTracker::photonSourceToProcessName = {{PhotonSummary::PhotonSource::Unknown, "Unknown"},
                                                                                                              {PhotonSummary::PhotonSource::Scintillation, "Scintillation"},
                                                                                                              {PhotonSummary::PhotonSource::Cherenkov, "Cerenkov"},
                                                                                                              {PhotonSummary::PhotonSource::OpWLS, "OpWLS"}};

const std::unordered_map<WLSLocation::WLSLoc, std::string> G4CCMTreeTracker::wlsLocationToProcessName = {{WLSLocation::WLSLoc::Unknown, "Unknown"},
                                                                                                             {WLSLocation::WLSLoc::PMT, "PMT"},
                                                                                                             {WLSLocation::WLSLoc::FoilTop, "FoilTop"},
                                                                                                             {WLSLocation::WLSLoc::FoilBottom, "FoilBottom"},
                                                                                                             {WLSLocation::WLSLoc::FoilSides, "FoilSides"}};

G4CCMTreeTracker::G4CCMTreeTracker(G4String name) : G4VSensitiveDetector(name) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMTreeTracker::Initialize(G4HCofThisEvent* hitsCE) {
    Reset();
    event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    primary_ = readout_->GetPrimary(event_id);
    mcTree = readout_->GetEDepMCTree(event_id);

    DaughterParticleMap.insert({1, primary_.GetID()});
    if(not I3MCTreeUtils::Has(*mcTree, primary_.GetID())) {
        I3MCTreeUtils::AddPrimary(*mcTree, primary_);
    }
}

void G4CCMTreeTracker::EndOfEvent(G4HCofThisEvent*) {
    for(std::pair<int const, std::tuple<double, std::vector<I3Particle>> > & sub_loss : sub_threshold_losses) {
        int track_id = sub_loss.first;
        double total_energy = std::get<0>(sub_loss.second);
        std::vector<I3Particle> & losses = std::get<1>(sub_loss.second);
        if(losses.size() == 0 or total_energy == 0)
            continue;

        I3Particle daughter = losses.back();
        I3Position position(0, 0, 0);
        I3Position direction(0, 0, 0);
        double time = 0.0;

        std::map<I3Particle::ParticleType, double> type_vote;

        for(I3Particle const & p : losses) {
            double const & e = p.GetEnergy();
            position += p.GetPos() * e;
            direction += p.GetDir() * e;
            time += p.GetTime() * e;
            std::map<I3Particle::ParticleType, double>::iterator it = type_vote.find(p.GetType());
            if(it == type_vote.end()) {
                it = std::get<0>(type_vote.insert(std::make_pair(p.GetType(), 0.0)));
            }
            it->second += e;
        }

        using pair_type = std::pair<I3Particle::ParticleType const, double>;

        I3Particle::ParticleType max_type = std::max_element(
            std::begin(type_vote),
            std::end(type_vote),
            [] (pair_type const & p1, pair_type const & p2) {
                return p1.second < p2.second;
            }
        )->first;

        position /= total_energy;
        direction /= total_energy;
        time /= total_energy;

        daughter.SetEnergy(total_energy);
        daughter.SetPos(position);
        daughter.SetDir(I3Direction(direction));
        daughter.SetTime(time);
        daughter.SetType(max_type);

        I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id), daughter);
    }

    if(DetailedPhotonTracking_) {
        readout_->LogTrackingResult(event_id, summary_map);
    } else {
        readout_->LogTrackingResult(event_id, source_map);
    }
}

void G4CCMTreeTracker::AddNewPhoton(int parent_id, int track_id, double time, double distance, double wavelength, std::string creation_process_name) {
    // We want to keep track of the parental lineage of any photon
    // But this photon does not come from another optical photon

    std::pair<I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator, bool> it = summary_map->insert({track_id, {ParentInfo(), PhotonSummary()}});

    PhotonSummary & s = std::get<1>(it.first->second);
    s.distance_uv = distance;
    s.original_wavelength = wavelength;
    // s.distance_visible = 0.0;
    s.time = time;
    //s.n_photons_per_wls = std::vector<size_t>();
    //s.wls_loc = WLSLocationSeries();
    s.photon_source = processNameToPhotonSource.at(creation_process_name);
    s.current_process = s.photon_source;
}

void G4CCMTreeTracker::AddPhotonTrack(int parent_id, int track_id, double delta_time, double delta_distance, std::string creation_process_name) {
    // Copy the existing entry from this parent track
    PhotonSummary const & parent_summary = std::get<1>(summary_map->at(parent_id));

    I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator it;
    bool b;
    std::tie(it, b) = summary_map->insert({track_id, {ParentInfo(), parent_summary}});

    PhotonSummary & this_photon_summary = std::get<1>(it->second);

    // Update the time, current process, and distance travelled
    this_photon_summary.time += delta_time;
    this_photon_summary.current_process = processNameToPhotonSource.at(creation_process_name);

    if(this_photon_summary.n_photons_per_wls.size() > 0)
        this_photon_summary.distance_visible += delta_distance;
    else
        this_photon_summary.distance_uv += delta_distance;
}

void G4CCMTreeTracker::UpdatePhoton(int parent_id, int track_id, double delta_time, double delta_distance, std::string creation_process_name) {
    // Grab the existing entry for this track ID
    I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator photon_summary_it = summary_map->find(track_id);
    if(photon_summary_it == summary_map->end()) {
        log_fatal("UpdatePhoton(%d,%d,%f,%f,%s): Could not find track_id %d in summary_map", parent_id, track_id, delta_time, delta_distance, creation_process_name.c_str(), track_id);
    }
    PhotonSummary & this_photon_summary = std::get<1>(photon_summary_it->second);

    // Update the time, current process, and distance travelled
    this_photon_summary.time += delta_time;
    this_photon_summary.current_process = processNameToPhotonSource.at(creation_process_name);

    if(this_photon_summary.n_photons_per_wls.size() > 0)
        this_photon_summary.distance_visible += delta_distance;
    else
        this_photon_summary.distance_uv += delta_distance;
}

void G4CCMTreeTracker::AddWLSPhotonTrack(int parent_id, int track_id, double delta_time, double delta_distance, G4Step* aStep) {
    // At the end of the day we really only care about photons that are hitting the PMTs
    // So we want to keep updating the photon's information until it hits a PMT

    // Copy the existing entry from this parent track
    I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator parent_it = summary_map->find(parent_id);
    if(parent_it == summary_map->end()) {
        log_fatal("AddWLSPhotonTrack(%d,%d,%f,%f,%p): Could not find parent_id %d in summary_map", parent_id, track_id, delta_time, delta_distance, static_cast<void*>(aStep), parent_id);
    }
    PhotonSummary const & parent_summary = std::get<1>(parent_it->second);

    bool b;
    I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator it;
    std::tie(it, b) = summary_map->insert({track_id, {ParentInfo(), parent_summary}});

    PhotonSummary & this_photon_summary = std::get<1>(it->second);
    this_photon_summary.time += delta_time;
    this_photon_summary.current_process = PhotonSummary::PhotonSource::OpWLS;

    // let's save where this wls occured
    std::string wls_loc_string = static_cast<std::string>(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());

    WLSLocation::WLSLoc wls_loc = WLSLocation::WLSLoc::Unknown;

    if(wls_loc_string.find("Coating") != std::string::npos) {
        wls_loc = WLSLocation::WLSLoc::PMT;
    } else if(wls_loc_string.find("FoilTop") != std::string::npos) {
        wls_loc = WLSLocation::WLSLoc::FoilTop;
    } else if(wls_loc_string.find("FoilBottom") != std::string::npos) {
        wls_loc = WLSLocation::WLSLoc::FoilBottom;
    } else if(wls_loc_string.find("FoilSides") != std::string::npos) {
        wls_loc = WLSLocation::WLSLoc::FoilSides;
    }

    // This is a new wavelength shift for this track so we need to keep track of the WLS locations
    size_t n_wls = this_photon_summary.n_photons_per_wls.size();

    n_wls += 1;
    if(n_wls == 1)
        this_photon_summary.distance_uv += delta_distance;
    else
        this_photon_summary.distance_visible += delta_distance;

    this_photon_summary.wls_loc.push_back(wls_loc);
    this_photon_summary.n_photons_per_wls.push_back(std::get<0>(parent_it->second).n_children);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CCMTreeTracker::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    // note -- this chunk of code resets global time to 0 in the case of radioactive decays
    // very important for retaining time structure of scintillation photons!!!
    G4Track* track = aStep->GetTrack();

    // Check if the particle has decayed
    // Check if it's a primary particle (ParentID == 0)
    if(track->GetParentID() == 0 and track->GetTrackStatus() == fStopAndKill) {
        G4String processName;
        const G4VProcess* currentProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep();
        if(currentProcess) {
            processName = static_cast<std::string>(currentProcess->GetProcessName());
        }
        if(processName == "Radioactivation") {
            // Get the list of secondaries
            const G4TrackVector* secondaries = aStep->GetSecondary();
            // Modify the start time of each secondary particle
            for(size_t i = 0; i < secondaries->size(); ++i) {
                G4Track* secondary = const_cast<G4Track*>(secondaries->at(i));
                secondary->SetGlobalTime(0.);
            }
        }
    }

    // Perform time cut if enabled
    if(TimeCut_) {
        G4double time = aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
        if(time > time_cut_value_) {
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return false;
        }
    }

    // Check the type of particle
    G4ParticleDefinition * particle_definition = aStep->GetTrack()->GetDefinition();

    // get pdg code
    G4int pdg = particle_definition->GetPDGEncoding();

    // Handle neutrinos
    if(KillNeutrinos_) {
        unsigned int abs_pdg = std::abs(int(pdg));
        bool pdg_even = abs_pdg % 2 == 0;
        bool is_nu = pdg_even and (abs_pdg >= 12 and abs_pdg <= 16);
        if(is_nu) {
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return false;
        }
    }

    // Handle optical photons
    if(particle_definition == G4OpticalPhoton::OpticalPhotonDefinition()) {
        if(KillPhotons_) {
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return false;
        }

        // check if we want to kill cerenkov photons
        if(KillCherenkov_ or KillScintillation_) {
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creation_process_name = "Unknown";
            if(creationProcess) {
                creation_process_name = static_cast<std::string>(creationProcess->GetProcessName());
            }
            if(KillCherenkov_ and (creation_process_name == "Cerenkov")) {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return false;
            } else if(KillScintillation_ and creation_process_name == "Scintillation") {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return false;
            }
        }

        if(DetailedPhotonTracking_) {
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creation_process_name = "Unknown";
            if(creationProcess) {
                creation_process_name = static_cast<std::string>(creationProcess->GetProcessName());
            }

            bool already_seen_photon = summary_map->find(aStep->GetTrack()->GetTrackID()) != summary_map->end();
            bool parent_seen = summary_map->find(aStep->GetTrack()->GetParentID()) != summary_map->end();
            bool brand_new_photon = not (already_seen_photon or parent_seen or creation_process_name == "OpWLS");
            bool photon_from_wls = (not already_seen_photon) and creation_process_name == "OpWLS";
            bool photon_from_photon = (not already_seen_photon) and parent_seen and creation_process_name != "OpWLS";

            size_t n_photon_secondaries = 0;
            std::vector<const G4Track*> const * secondaries = aStep->GetSecondaryInCurrentStep();
            if(secondaries != nullptr) {
                size_t n_secondaries = secondaries->size();
                for(size_t i = 0; i < n_secondaries; ++i) {
                    G4Track const * secondary = secondaries->at(i);
                    if(secondary->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
                        n_photon_secondaries += 1;
                    }
                }
            }
            bool has_photon_secondaries = n_photon_secondaries > 0;
            bool stopping = aStep->GetTrack()->GetTrackStatus() == fStopAndKill;

            // Child processing
            bool child_case_1_brand_new_photon = brand_new_photon;
            bool child_case_2_existing_photon = already_seen_photon;
            bool child_case_3_new_photon_from_wls = photon_from_wls;
            bool child_case_4_new_photon_from_photon = photon_from_photon;

            G4int parent_id = aStep->GetTrack()->GetParentID();
            G4int track_id = aStep->GetTrack()->GetTrackID();

            double g4_delta_distance = (aStep->GetStepLength() / mm) * I3Units::mm;
            double pre_step_global_time = aStep->GetPreStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
            double g4_delta_time_step = aStep->GetDeltaTime() / nanosecond * I3Units::nanosecond;

            double wavelength = hc / (aStep->GetTrack()->GetTotalEnergy() / electronvolt); // units of nanometer
            double original_wavelength = wavelength * I3Units::nanometer;

            if(child_case_1_brand_new_photon) {
                AddNewPhoton(parent_id, track_id, g4_delta_time_step, g4_delta_distance, original_wavelength, creation_process_name);
            } else if(child_case_2_existing_photon) {
                UpdatePhoton(parent_id, track_id, g4_delta_time_step, g4_delta_distance, creation_process_name);
            } else if(child_case_3_new_photon_from_wls) {
                AddWLSPhotonTrack(parent_id, track_id, g4_delta_time_step, g4_delta_distance, aStep);
                // Check if we still need the parent
                I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator it = summary_map->find(parent_id);
                if(it == summary_map->end()) {
                    log_fatal("ProcessHits child_case_3_new_photon_from_wls: Could not find parent_id %d in summary_map", parent_id);
                }
                ParentInfo & info = std::get<0>(it->second);
                info.n_children_remaining -= 1;
                if(info.n_children_remaining == 0) {
                    // this parent is done
                    summary_map->erase(parent_id);
                }
            } else if(child_case_4_new_photon_from_photon) {
                AddPhotonTrack(parent_id, track_id, g4_delta_time_step, g4_delta_distance, creation_process_name);
            }

            bool parent_case_1_photon_is_wls = has_photon_secondaries and stopping;
            bool parent_case_2_photon_not_stopping = has_photon_secondaries and (not stopping);
            if(parent_case_1_photon_is_wls) {
                I3Map<int, std::tuple<ParentInfo, PhotonSummary>>::iterator it = summary_map->find(track_id);

                ParentInfo & info = std::get<0>(it->second);
                info.n_children = n_photon_secondaries;
                info.n_children_remaining = n_photon_secondaries;
            }

            assert(not parent_case_2_photon_not_stopping);
        } else {
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creation_process_name = "Unknown";
            if(creationProcess) {
                creation_process_name = static_cast<std::string>(creationProcess->GetProcessName());
            }

            bool already_seen_photon = source_map->find(aStep->GetTrack()->GetTrackID()) != source_map->end();
            bool parent_seen = source_map->find(aStep->GetTrack()->GetParentID()) != source_map->end();
            bool brand_new_photon = not (already_seen_photon or parent_seen or creation_process_name == "OpWLS");
            bool photon_from_wls = (not already_seen_photon) and creation_process_name == "OpWLS";
            bool photon_from_photon = (not already_seen_photon) and parent_seen and creation_process_name != "OpWLS";

            size_t n_photon_secondaries = 0;
            std::vector<const G4Track*> const * secondaries = aStep->GetSecondaryInCurrentStep();
            if(secondaries != nullptr) {
                size_t n_secondaries = secondaries->size();
                for(size_t i = 0; i < n_secondaries; ++i) {
                    G4Track const * secondary = secondaries->at(i);
                    if(secondary->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
                        n_photon_secondaries += 1;
                    }
                }
            }
            bool has_photon_secondaries = n_photon_secondaries > 0;
            bool stopping = aStep->GetTrack()->GetTrackStatus() == fStopAndKill;

            // Child processing
            bool child_case_1_brand_new_photon = brand_new_photon;
            bool child_case_2_existing_photon = already_seen_photon;
            bool child_case_3_new_photon_from_wls = photon_from_wls;
            bool child_case_4_new_photon_from_photon = photon_from_photon;

            G4int parent_id = aStep->GetTrack()->GetParentID();
            G4int track_id = aStep->GetTrack()->GetTrackID();

            double g4_delta_distance = (aStep->GetStepLength() / mm) * I3Units::mm;
            double pre_step_global_time = aStep->GetPreStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
            double g4_delta_time_step = aStep->GetDeltaTime() / nanosecond * I3Units::nanosecond;

            double wavelength = hc / (aStep->GetTrack()->GetTotalEnergy() / electronvolt); // units of nanometer
            double original_wavelength = wavelength * I3Units::nanometer;

            if(child_case_1_brand_new_photon) {
                source_map->insert({track_id, {ParentInfo(), processNameToPhotonSource.at(creation_process_name)}});
                //photon_source_map->insert({track_id, processNameToPhotonSource.at(creation_process_name)});
                assert(creation_process_name != "OpWLS");
            } else if(child_case_2_existing_photon) {
            } else if(child_case_3_new_photon_from_wls) {
                I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>::iterator it = source_map->find(parent_id);
                if(it == source_map->end()) {
                    log_fatal("ProcessHits child_case_3_new_photon_from_wls: Could not find parent_id %d in source_map", parent_id);
                }
                source_map->insert({track_id, {ParentInfo(), std::get<1>(it->second)}});

                ParentInfo & info = std::get<0>(it->second);
                assert(std::get<1>(it->second) != PhotonSummary::PhotonSource::OpWLS);
                // Check if we still need the parent
                info.n_children_remaining -= 1;
                if(info.n_children_remaining == 0) {
                    // this parent is done
                    source_map->erase(parent_id);
                }
            } else if(child_case_4_new_photon_from_photon) {
                I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>::iterator it = source_map->find(parent_id);
                if(it == source_map->end()) {
                    log_fatal("ProcessHits child_case_34new_photon_from_photon: Could not find parent_id %d in source_map", parent_id);
                }
                source_map->insert({track_id, {ParentInfo(), std::get<1>(it->second)}});
                assert(std::get<1>(source_map->at(parent_id)) != PhotonSummary::PhotonSource::OpWLS);
                source_map->erase(parent_id);
            }

            bool parent_case_1_photon_is_wls = has_photon_secondaries and stopping;
            bool parent_case_2_photon_not_stopping = has_photon_secondaries and (not stopping);
            if(parent_case_1_photon_is_wls) {
                I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>::iterator it = source_map->find(track_id);
                if(it == source_map->end()) {
                    log_fatal("ProcessHit parent_case_1_photon_is_wls: Could not find track_id %d in source_map", track_id);
                }
                ParentInfo & info = std::get<0>(it->second);
                info.n_children = n_photon_secondaries;
                info.n_children_remaining = n_photon_secondaries;
            }

            assert(not parent_case_2_photon_not_stopping);
        }

        /*
        if(DetailedPhotonTracking_) {
            // ok if we survived all checks, add an entry to our photon summary
            // we need parent id, track id, distance travelled, and wavelength
            G4int parent_id = aStep->GetTrack()->GetParentID();
            G4int track_id = aStep->GetTrack()->GetTrackID();

            double g4_delta_distance = (aStep->GetStepLength() / mm) * I3Units::mm;
            double pre_step_global_time = aStep->GetPreStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
            double g4_delta_time_step = aStep->GetDeltaTime() / nanosecond * I3Units::nanosecond;

            double wavelength = hc / (aStep->GetTrack()->GetTotalEnergy() / electronvolt); // units of nanometer
            double original_wavelength = wavelength * I3Units::nanometer;

            G4VProcess const * process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
            std::string processName = (process) ? process->GetProcessName() : "Unknown";
            if(processName == "OpWLS") {
                ParentInfo info;
                std::vector<const G4Track*> * secondaries = aStep->GetSecondaryInCurrentStep();
                if(secondaries != nullptr) {
                    size_t n_secondaries = secondaries->size();
                    if(n_secondaries > 0) {
                        // this is a WLS photon that produced some secondaries
                        info.n_children = n_secondaries;
                        info.n_children_remaining = n_secondaries;
                    }
                }
            }

            // let's also check if this photon got wls
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creation_process_name = "Unknown";
            if(creationProcess) {
                creation_process_name = static_cast<std::string>(creationProcess->GetProcessName());
            }

            // ok now let's save
            // check to see if the track id is in our map
            std::map<int, PhotonSummary>::iterator it = photon_summary->find(track_id);

            bool already_seen_track = it != photon_summary->end();

            if(already_seen_track) {
                // child case 2
                UpdatePhoton(parent_id, track_id, g4_delta_time_step, g4_delta_distance, creation_process_name);
            } else {
                // check if this parent id is in our map
                std::map<int, PhotonSummary>::iterator parent_it = photon_summary->find(parent_id);
                bool already_seen_parent = parent_it != photon_summary->end();
                bool is_wavelength_shift = creation_process_name == "OpWLS";
                if(already_seen_parent) {
                    // Need to add a new photon track and update any siblings
                    if(is_wavelength_shift)
                        // child case 3
                        AddWLSPhotonTrack(parent_id, track_id, g4_delta_time_step, g4_delta_distance, aStep);
                    else
                        // child case 5
                        AddPhotonTrack(parent_id, track_id, g4_delta_time_step, g4_delta_distance, creation_process_name);
                } else {
                    // child case 1
                    // Just need to add a brand new photon track
                    AddNewPhoton(parent_id, track_id, g4_delta_time_step, g4_delta_distance, original_wavelength, creation_process_name);
                }

                assert(num_children[parent_id] > 0);
                assert((*num_siblings[track_id]) > 0);
                num_children[parent_id] -= 1;
                (*num_siblings[track_id]) -= 1;
                if(num_children[parent_id] == 0 and (*num_siblings[parent_id]) == 0) {
                    photon_summary->erase(parent_id);
                    num_children.erase(parent_id);
                    num_siblings.erase(parent_id);
                }
            }
        }
        */
        return false;
    }

    if(not (TrackParticles_ or DetailedPhotonTracking_ or TrackEnergyLosses_)) {
        return false;
    }

    // position
    G4ThreeVector g4Position = aStep->GetPostStepPoint()->GetPosition();
    I3Position position(g4Position.x() / mm * I3Units::mm, g4Position.y() / mm * I3Units::mm, g4Position.z() / mm * I3Units::mm);

    // direction
    G4ThreeVector g4Direction = aStep->GetPostStepPoint()->GetMomentumDirection();
    I3Position direction(g4Direction.x(), g4Direction.y(), g4Direction.z());

    // time
    double time = primary_.GetTime() + aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;

    // length
    double length = aStep->GetStepLength() / mm * I3Units::mm;

    // let's also grab parent id
    // if parent id == 0, that's our primary injected particle
    G4int parent_id = aStep->GetTrack()->GetParentID();
    G4int track_id = aStep->GetTrack()->GetTrackID();


    G4VProcess const * process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    std::string processName = (process) ? process->GetProcessName() : "Unknown";

    // now save to our MCTree!
    if(parent_id > 0) {
        // ok so we've created a new particle
        // if this is the first time we're seeing this particle -- add particle + energy loss
        // if we've already added this daughter particle -- only add energy loss

        if(DaughterParticleMap.find(track_id) == DaughterParticleMap.end()) {
            I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(pdg);
            I3Particle daughter(daughter_type);

            std::map<int, int>::iterator parent_it = parent_map.find(parent_id);
            if(parent_it != parent_map.end())
                parent_id = parent_it->second;

            // we have not added the daugher...let's do it now
            double energy = aStep->GetTrack()->GetVertexKineticEnergy() / electronvolt * I3Units::eV;
            if(energy < G4ETrackingMin_) {
                parent_map.insert({track_id, parent_id});
            } else {
                daughter.SetEnergy(energy);
                daughter.SetPos(position);
                daughter.SetDir(I3Direction(direction));
                daughter.SetTime(time);
                daughter.SetLength(length);
                I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(parent_id), daughter);
            }

            // update map
            DaughterParticleMap.insert({track_id, daughter.GetID()});
        } else {
            I3MCTree::iterator it = mcTree->find(DaughterParticleMap.at(track_id));
            if(it != mcTree->end()) {
                // All particle properties are defined at creation,
                // except length
                I3Particle & daughter = *it;
                daughter.SetLength(daughter.GetLength() + length);
            }
        }
    }

    if(TrackEnergyLosses_) {
        double edep = aStep->GetTotalEnergyDeposit() / electronvolt * I3Units::eV;
        if(edep == 0.0)
            return false;

        I3Particle::ParticleType daughter_type;
        // Check if an energy loss has occurred
        if(energyLossToI3ParticlePDGCode.find(processName) == energyLossToI3ParticlePDGCode.end()) {
            if(knownProcessNames.find(processName) == knownProcessNames.end()) {
                G4cout << "oops! no conversion for " << processName << std::endl;
            }
            daughter_type = I3Particle::unknown;
        } else {
            daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(processName));
        }

        // now add energy loss
        I3Particle daughter(daughter_type);

        daughter.SetEnergy(edep);
        daughter.SetPos(position);
        daughter.SetDir(I3Direction(direction));
        daughter.SetTime(time);
        daughter.SetLength(length);

        std::map<int, int>::iterator parent_it = parent_map.find(track_id);
        if(parent_it != parent_map.end())
            track_id = parent_it->second;

        if(DaughterParticleMap.find(track_id) == DaughterParticleMap.end()) {
            G4cout << "oops! trying to save energy deposition type but DaughterParticleMap does not have track id!" << std::endl;
        } else {
            if(edep < G4EDepMin_) {
                std::map<int, std::tuple<double, std::vector<I3Particle>>>::iterator sub_loss_it = sub_threshold_losses.find(track_id);
                if(sub_loss_it == sub_threshold_losses.end()) {
                    sub_loss_it = sub_threshold_losses.insert(sub_loss_it, {track_id, {0.0, std::vector<I3Particle>()}});
                }
                std::tuple<double, std::vector<I3Particle>> & sub_losses = sub_loss_it->second;
                double & total_energy = std::get<0>(sub_losses);
                std::vector<I3Particle> & losses = std::get<1>(sub_losses);
                std::map<I3Particle::ParticleType, double> type_vote;
                total_energy += edep;
                if(total_energy > G4EDepMin_) {
                    type_vote.insert(std::make_pair(daughter.GetType(), edep));
                    position = position * edep;
                    direction *= edep;
                    time *= edep;
                    for(I3Particle const & p : losses) {
                        double const & e = p.GetEnergy();
                        position += p.GetPos() * e;
                        direction += p.GetDir() * e;
                        time += p.GetTime() * e;
                        length += p.GetLength();
                        std::map<I3Particle::ParticleType, double>::iterator it = type_vote.find(p.GetType());
                        if(it == type_vote.end()) {
                            it = std::get<0>(type_vote.insert(std::make_pair(p.GetType(), 0.0)));
                        }
                        it->second += e;
                    }
                    position /= total_energy;
                    direction /= total_energy;
                    time /= total_energy;

                    using pair_type = std::pair<I3Particle::ParticleType const, double>;

                    I3Particle::ParticleType max_type = std::max_element(
                        std::begin(type_vote),
                        std::end(type_vote),
                        [] (pair_type const & p1, pair_type const & p2) {
                            return p1.second < p2.second;
                        }
                    )->first;

                    // ok we have a valid energy deposition -- let's add it to the tree
                    daughter.SetEnergy(total_energy);
                    daughter.SetPos(position);
                    daughter.SetDir(I3Direction(direction));
                    daughter.SetTime(time);
                    daughter.SetLength(length);
                    daughter.SetType(max_type);
                    I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id), daughter);
                    total_energy = 0.0;
                    losses.clear();
                } else {
                    losses.push_back(daughter);
                }
            } else {
                I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id), daughter);
            }
        }
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
