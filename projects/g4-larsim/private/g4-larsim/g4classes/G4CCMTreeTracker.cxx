
#include "g4-larsim/g4classes/G4CCMTreeTracker.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"

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

const std::unordered_set<int> G4CCMTreeTracker::energyLossPDGCodes = {0, 2000000001, 2000000002, 2000000003, 2000000004, 2000000005, 2000000006, 2000000007, 2000000008, 2000000009, 2000000010, 2000000011, 2000000012, 2000000013, 2000000014, 2000000015, 2000000016};

const std::unordered_map<std::string, int> G4CCMTreeTracker::energyLossToI3ParticlePDGCode = {{"phot", 2000000001}, {"compt", 2000000002}, {"conv", 2000000003},
                                                                                          {"Rayl", 2000000004}, {"msc", 2000000005}, {"eIoni", 2000000006},
                                                                                          {"eBrem", 2000000007}, {"ePairProd", 2000000008}, {"CoulombScat", 2000000009},
                                                                                          {"annihil", 2000000010}, {"Cerenkov", 2000000011}, {"Radioactivation", 2000000012},
                                                                                          {"Scintillation", 2000000013}, {"OpWLS", 2000000014}, {"ionIoni" , 2000000015},
                                                                                          {"hIoni", 2000000016}, {"Unknown", 0}};

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

    DaughterParticleMap[1] = primary_.GetID();
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

        for(I3Particle const & p : losses) {
            double const & e = p.GetEnergy();
            position += p.GetPos() * e;
            direction += p.GetDir() * e;
            time += p.GetTime() * e;
        }

        position /= total_energy;
        direction /= total_energy;
        time /= total_energy;

        daughter.SetEnergy(total_energy);
        daughter.SetPos(position);
        daughter.SetDir(I3Direction(direction));
        daughter.SetTime(time);

        I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id), daughter);
    }

    readout_->LogTrackingResult(event_id, photon_summary, optical_photon_map, DetailedPhotonTracking_);
}

void G4CCMTreeTracker::AddNewPhoton(int parent_id, int track_id, double time, double distance, double wavelength, std::string creation_process_name) {
    // We want to keep track of the parental lineage of any photon
    // But this photon does not come from another optical photon
    (*wls_daughter_parent_map)[track_id] = parent_id;

    std::vector<size_t> n_photons_per_wls = {};
    PhotonSummary this_photon_summary = PhotonSummary(distance, // distance travelled before wls
                                                      wavelength,
                                                      0.0, // distance travelled after wls
                                                      time,
                                                      n_photons_per_wls, WLSLocationSeries(),
                                                      processNameToPhotonSource.at(creation_process_name),
                                                      processNameToPhotonSource.at(creation_process_name));
    photon_summary->push_back(this_photon_summary);
    optical_photon_map->insert(std::make_pair(track_id, photon_summary->size() - 1));
}

void G4CCMTreeTracker::AddPhotonTrack(int parent_id, int track_id, size_t parent_index, double delta_time, double delta_distance, std::string creation_process_name) {
    // We want to keep track of the parental lineage of any photon
    (*wls_daughter_parent_map)[track_id] = parent_id;

    // Copy the existing entry from this parent track
    PhotonSummary this_photon_summary = photon_summary->at(parent_index);

    // Update the time, current process, and distance travelled
    this_photon_summary.time += delta_time;
    this_photon_summary.current_process = processNameToPhotonSource.at(creation_process_name);

    if(this_photon_summary.n_photons_per_wls.size() > 0)
        this_photon_summary.distance_visible += delta_distance;
    else
        this_photon_summary.distance_uv += delta_distance;

    // Now add the new track to the vector of photon summaries
    photon_summary->push_back(this_photon_summary);
    size_t new_pos = photon_summary->size() - 1;
    (*optical_photon_map)[track_id] = new_pos;
}

void G4CCMTreeTracker::UpdatePhoton(int parent_id, int track_id, size_t photon_index, double delta_time, double delta_distance, std::string creation_process_name) {
    // Technically we should never end up here since each step has a unique track ID

    // Grab the existing entry for this track ID
    PhotonSummary & this_photon_summary = photon_summary->at(photon_index);

    // Update the time, current process, and distance travelled
    this_photon_summary.time += delta_time;
    this_photon_summary.current_process = processNameToPhotonSource.at(creation_process_name);

    if(this_photon_summary.n_photons_per_wls.size() > 0)
        this_photon_summary.distance_visible += delta_distance;
    else
        this_photon_summary.distance_uv += delta_distance;
}

void G4CCMTreeTracker::AddWLSPhotonTrack(int parent_id, int track_id, size_t parent_index, double delta_time, double delta_distance, G4Step* aStep) {
    // At the end of the day we really only care about photons that are hitting the PMTs
    // So we want to keep updating the photon's information until it hits a PMT

    // We want to keep track of the parental lineage of any photon
    (*wls_daughter_parent_map)[track_id] = parent_id;

    // Copy the existing entry from this parent track
    PhotonSummary this_photon_summary = photon_summary->at(parent_index);
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

    // let's update the parent id and track id map
    // this keeps track of parent id and wls daughter track ids
    std::map<int, std::set<int>>::iterator wls_it = wls_parent_daughter_map->find(parent_id);
    if(wls_it != wls_parent_daughter_map->end()) {
        // ok this parent id is in our map! let's update daughter track ids
        wls_it->second.insert(track_id);
    } else {
        // this key is NOT in our map!! let's add a value
        (*wls_parent_daughter_map)[parent_id] = std::set<int> {track_id};
    }

    // now let's update n_photons_per_wls for all daughers tracks of this parent
    std::set<int> & sibling_track_ids = (*wls_parent_daughter_map)[parent_id];
    size_t n_siblings = sibling_track_ids.size();

    // This is a new wavelength shift for this track so we need to keep track of the WLS locations
    size_t n_wls = this_photon_summary.n_photons_per_wls.size();

    n_wls += 1;
    if(n_wls == 1)
        this_photon_summary.distance_uv += delta_distance;
    else
        this_photon_summary.distance_visible += delta_distance;

    if(this_photon_summary.n_photons_per_wls.size() >= n_wls) {
        this_photon_summary.n_photons_per_wls.at(n_wls-1) = n_siblings;
        this_photon_summary.wls_loc.at(n_wls-1) = wls_loc;
    } else {
        for(size_t n = this_photon_summary.n_photons_per_wls.size() + 1; n < n_wls - 1; ++n) {
            this_photon_summary.n_photons_per_wls.push_back(0);
            this_photon_summary.wls_loc.push_back(WLSLocation::WLSLoc::Unknown);
        }
        this_photon_summary.n_photons_per_wls.push_back(n_siblings);
        this_photon_summary.wls_loc.push_back(wls_loc);
    }

    // now update the photon_summary, delete from map, and add new entry to the map
    photon_summary->push_back(this_photon_summary);
    size_t new_pos = photon_summary->size() - 1;
    (*optical_photon_map)[track_id] = new_pos;

    std::deque<int> fifo(sibling_track_ids.begin(), sibling_track_ids.end());

    // loop over daughter tracks
    while(not fifo.empty()) {
        int sibling_track_id = fifo.front();
        fifo.pop_front();

        PhotonSummary & sibling_photon_summary = photon_summary->at((*optical_photon_map)[sibling_track_id]);

        int parent_id = (*wls_daughter_parent_map)[sibling_track_id];
        size_t n_siblings = (*wls_parent_daughter_map)[parent_id].size();
        size_t current_sibling_n_wls = sibling_photon_summary.wls_loc.size();

        if(current_sibling_n_wls >= n_wls) {
            sibling_photon_summary.n_photons_per_wls.at(n_wls - 1) = n_siblings;
            sibling_photon_summary.wls_loc.at(n_wls - 1) = wls_loc;
        } else {
            for(size_t n = current_sibling_n_wls + 1; n < n_wls - 1; ++n) {
                sibling_photon_summary.n_photons_per_wls.push_back(0);
                sibling_photon_summary.wls_loc.push_back(WLSLocation::WLSLoc::Unknown);
            }
            sibling_photon_summary.n_photons_per_wls.push_back(n_siblings);
            sibling_photon_summary.wls_loc.push_back(wls_loc);
        }

        std::set<int> const & sibling_track_ids = (*wls_parent_daughter_map)[parent_id];

        fifo.insert(fifo.end(), sibling_track_ids.cbegin(), sibling_track_ids.cend());
    }
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
            if(KillCherenkov_ and creation_process_name == "Cerenkov") {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return false;
            } else if(KillScintillation_ and creation_process_name == "Scintillation") {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return false;
            }
        }

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

            // let's also check if this photon got wls
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creation_process_name = "Unknown";
            if(creationProcess) {
                creation_process_name = static_cast<std::string>(creationProcess->GetProcessName());
            }

            // ok now let's save
            // check to see if the track id is in our map
            std::map<int, size_t>::iterator it = optical_photon_map->find(track_id);

            bool already_seen_track = it != optical_photon_map->end();

            if(already_seen_track) {
                UpdatePhoton(parent_id, track_id, it->second, g4_delta_time_step, g4_delta_distance, creation_process_name);
            } else {
                // check if this parent id is in our map
                std::map<int, size_t>::iterator parent_it = optical_photon_map->find(parent_id);
                bool already_seen_parent = parent_it != optical_photon_map->end();
                bool is_wavelength_shift = creation_process_name == "OpWLS";
                if(already_seen_parent) {
                    // Need to add a new photon track and update any siblings
                    if(is_wavelength_shift)
                        AddWLSPhotonTrack(parent_id, track_id, parent_it->second, g4_delta_time_step, g4_delta_distance, aStep);
                    else
                        AddPhotonTrack(parent_id, track_id, parent_it->second, g4_delta_time_step, g4_delta_distance, creation_process_name);
                } else {
                    // Just need to add a brand new photon track
                    AddNewPhoton(parent_id, track_id, g4_delta_time_step, g4_delta_distance, original_wavelength, creation_process_name);
                }
            }
        }
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
                parent_map[track_id] = parent_id;
            } else {
                daughter.SetEnergy(energy);
                daughter.SetPos(position);
                daughter.SetDir(I3Direction(direction));
                daughter.SetTime(time);
                I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(parent_id) , daughter);
            }

            // update map
            DaughterParticleMap[track_id] = daughter.GetID();
        }
    }

    if(TrackEnergyLosses_) {
        // Check if an energy loss has occurred
        if(energyLossToI3ParticlePDGCode.find(processName) != energyLossToI3ParticlePDGCode.end()) {
            // now add energy loss
            I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(processName));
            I3Particle daughter(daughter_type);

            double edep = aStep->GetTotalEnergyDeposit() / electronvolt * I3Units::eV;
            daughter.SetEnergy(edep);
            daughter.SetPos(position);
            daughter.SetDir(I3Direction(direction));
            daughter.SetTime(time);

            std::map<int, int>::iterator parent_it = parent_map.find(track_id);
            if(parent_it != parent_map.end())
                track_id = parent_it->second;

            if(DaughterParticleMap.find(track_id) == DaughterParticleMap.end()) {
                G4cout << "oops! trying to save energy deposition type but DaughterParticleMap does not have track id!" << std::endl;
            } else {
                if(edep < G4EDepMin_) {
                    if(sub_threshold_losses.find(track_id) == sub_threshold_losses.end()) {
                        sub_threshold_losses[track_id] = {0.0, std::vector<I3Particle>()};
                    }
                    std::tuple<double, std::vector<I3Particle>> & sub_losses = sub_threshold_losses[track_id];
                    double & total_energy = std::get<0>(sub_losses);
                    std::vector<I3Particle> & losses = std::get<1>(sub_losses);
                    total_energy += edep;
                    if(total_energy > G4EDepMin_) {
                        position = position * edep;
                        direction *= edep;
                        time *= edep;
                        for(I3Particle const & p : losses) {
                            double const & e = p.GetEnergy();
                            position += p.GetPos() * e;
                            direction += p.GetDir() * e;
                            time += p.GetTime() * e;
                        }
                        position /= total_energy;
                        direction /= total_energy;
                        time /= total_energy;

                        // ok we have a valid energy deposition -- let's add it to the tree
                        daughter.SetEnergy(total_energy);
                        daughter.SetPos(position);
                        daughter.SetDir(I3Direction(direction));
                        daughter.SetTime(time);
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
        } else {
            G4cout << "oops! no conversion for " << processName << std::endl;
        }
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
