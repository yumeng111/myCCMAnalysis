
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
                                                                                                              {"Cerenkov", PhotonSummary::PhotonSource::Cerenkov},
                                                                                                              {"OpWLS", PhotonSummary::PhotonSource::OpWLS}};

const std::unordered_map<PhotonSummary::PhotonSource, std::string> G4CCMTreeTracker::photonSourceToProcessName = {{PhotonSummary::PhotonSource::Unknown, "Unknown"},
                                                                                                              {PhotonSummary::PhotonSource::Scintillation, "Scintillation"},
                                                                                                              {PhotonSummary::PhotonSource::Cerenkov, "Cerenkov"},
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
    event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    primary_ = readout_->GetPrimary(event_id);
    mcTree = readout_->GetEDepMCTree(event_id);

    DaughterParticleMap[1] = primary_.GetID();
    if(not I3MCTreeUtils::Has(*mcTree, primary_.GetID())) {
        I3MCTreeUtils::AddPrimary(*mcTree, primary_);
    }
}

void G4CCMTreeTracker::EndOfEvent(G4HCofThisEvent*) {
    readout_->LogTrackingResult(event_id, photon_summary, optical_photon_map, DetailedPhotonTracking_);
    Reset();
}

void G4CCMTreeTracker::AddEntryToPhotonSummary(int parent_id, int track_id, double g4_delta_distance, double original_wavelength,
                                               double g4_time, std::string creationProcessName) {
    // map does not have key -- let's add our PhotonSummary then update map
    size_t n_wls = 0;
    std::vector<size_t> n_photons_per_wls = {0};
    PhotonSummary this_photon_summary = PhotonSummary(g4_delta_distance, // distance travelled before wls
                                                      original_wavelength,
                                                      0.0, // distance travelled after wls
                                                      g4_time,
                                                      n_wls, n_photons_per_wls, WLSLocationSeries(),
                                                      processNameToPhotonSource.at(creationProcessName),
                                                      processNameToPhotonSource.at(creationProcessName),
                                                      processNameToPhotonSource.at(creationProcessName));
    photon_summary->push_back(this_photon_summary);
    optical_photon_map->insert(std::make_pair(track_id, photon_summary->size() - 1));
}

void G4CCMTreeTracker::UpdatePhotonSummary(int parent_id, int track_id, double g4_time, std::string creationProcessName,
                                           std::map<int, size_t>::iterator it, bool new_process, G4Step* aStep, double g4_delta_distance) {
    // map contains this photon
    // so we need to grab PhotonSummary and update it -- then update key
    PhotonSummary this_photon_summary = photon_summary->at(it->second);
    this_photon_summary.g4_time += g4_time;
    this_photon_summary.current_process = processNameToPhotonSource.at(creationProcessName);

    if (new_process){
        this_photon_summary.temp_parent = photon_summary->at(optical_photon_map->find(parent_id)->second).current_process;
    }

    if(creationProcessName == "OpWLS" and new_process) {
        this_photon_summary.n_wls += 1;

        // let's save where this wls occured
        std::string wls_loc_string = static_cast<std::string>(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());

        WLSLocation::WLSLoc wls_loc = WLSLocation::WLSLoc::Unknown;

        if (wls_loc_string.find("Coating") != std::string::npos) {
            wls_loc = WLSLocation::WLSLoc::PMT;
        } else if (wls_loc_string.find("FoilTop") != std::string::npos) {
            wls_loc = WLSLocation::WLSLoc::FoilTop;
        } else if (wls_loc_string.find("FoilBottom") != std::string::npos) {
            wls_loc = WLSLocation::WLSLoc::FoilBottom;
        } else if (wls_loc_string.find("FoilSides") != std::string::npos) {
            wls_loc = WLSLocation::WLSLoc::FoilSides;
        }


        size_t saving_idx = this_photon_summary.n_wls - 1;
        if (this_photon_summary.wls_loc.size() > saving_idx){
             this_photon_summary.wls_loc.at(saving_idx) = wls_loc;
        } else {
            this_photon_summary.wls_loc.push_back(wls_loc);
        }
    }

    // now let's add our distance travelled based on whether or not we've wls
    if (this_photon_summary.n_wls == 1 and new_process and creationProcessName == "OpWLS"){
        // special case where wls triggered step, that means the distance travelled is before wls
        this_photon_summary.g4_distance_uv += g4_delta_distance;
    } else if (this_photon_summary.n_wls > 0) {
        // we're already wls at least once
        this_photon_summary.g4_distance_visible += g4_delta_distance;
    } else {
        // not yet wls!
        this_photon_summary.g4_distance_uv += g4_delta_distance;
    }

    // now update the photon_summary, delete from map, and add new entry to the map
    photon_summary->push_back(this_photon_summary);
    size_t new_pos = photon_summary->size() - 1;
    (*optical_photon_map)[track_id] = new_pos;

    // let's also update n photons per wls
    if (creationProcessName == "OpWLS" and new_process){

        // let's update the parent id and track id map
        // this keeps track of parent id and wls daughter track ids
        std::map<int, std::vector<int>>::iterator wls_it = wls_parent_daughter_map->find(parent_id);
        if (wls_it != wls_parent_daughter_map->end()) {
            // ok this parent id is in our map! let's update daughter track ids
            wls_it->second.push_back(track_id);
        } else {
            // this key is NOT in our map!! let's add a value
            (*wls_parent_daughter_map)[parent_id] = std::vector<int> {track_id};
        }

        // now let's update n_photons_per_wls for all daughers tracks of this parent
        std::vector<int> tracks_per_parent = (*wls_parent_daughter_map)[parent_id];

        // loop over daughter tracks
        for (size_t t = 0; t < tracks_per_parent.size(); t++){
            // grab n_wls this track has undergone
            size_t this_n_wls =  photon_summary->at((*optical_photon_map)[tracks_per_parent[t]]).n_wls;

            // looping over from current to original wls
            int prev_track;
            for (size_t n = this_n_wls; n > 0; n--){
                // special logic for first iteration
                if (n == this_n_wls) {
                    size_t this_n_photons_per_wls = (*wls_parent_daughter_map)[parent_id].size();

                    // now update our photon summary
                    size_t s = photon_summary->at((*optical_photon_map)[tracks_per_parent[t]]).n_photons_per_wls.size();
                    if (n <= s){
                        photon_summary->at((*optical_photon_map)[tracks_per_parent[t]]).n_photons_per_wls.at(n-1) = this_n_photons_per_wls;
                    } else {
                        for (size_t nn = s; nn < n; nn++){
                            photon_summary->at((*optical_photon_map)[tracks_per_parent[t]]).n_photons_per_wls.push_back(this_n_photons_per_wls);
                        }
                    }

                    prev_track = parent_id;

                } else {

                    // ok our "track id"  is prev_track, let's find the corresponding parent id
                    int prev_parent;
                    for (auto it = wls_parent_daughter_map->begin(); it != wls_parent_daughter_map->end(); ++it) {
                        std::vector<int>& vec = it->second;
                        if (std::find(vec.begin(), vec.end(), prev_track) != vec.end()) {
                            prev_parent = it->first;
                        }
                    }

                    // ok now let's update our photon summary
                    size_t this_n_photons_per_wls = (*wls_parent_daughter_map)[prev_parent].size();
                    photon_summary->at((*optical_photon_map)[tracks_per_parent[t]]).n_photons_per_wls.at(n-1) = this_n_photons_per_wls;

                    prev_track = prev_parent;
                }
            }
        }
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
        if (currentProcess) {
            processName = static_cast<std::string>(currentProcess->GetProcessName());
        }
        if(processName == "Radioactivation") {
            // Get the list of secondaries
            const G4TrackVector* secondaries = aStep->GetSecondary();
            // Modify the start time of each secondary particle
            for (size_t i = 0; i < secondaries->size(); ++i) {
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
            std::string creationProcessName = "Unknown";
            if(creationProcess) {
                creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
            }
            if(KillCherenkov_ and creationProcessName == "Cerenkov") {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return false;
            } else if (KillScintillation_ and creationProcessName == "Scintillation") {
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
            std::string creationProcessName = "Unknown";
            if (creationProcess) {
                creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
            }

            // ok now let's save
            // check to see if the track id is in our map
            bool new_process = false;
            std::map<int, size_t>::iterator it = optical_photon_map->find(track_id);

            if (it != optical_photon_map->end()) {
                // ok so this photon is in our map, let's just update
                UpdatePhotonSummary(parent_id, track_id, g4_delta_time_step, creationProcessName, it, new_process, aStep, g4_delta_distance);
            } else {
                // check if this parent id is in our map
                std::map<int, size_t>::iterator parent_it = optical_photon_map->find(parent_id);

                bool new_process = true;
                if (parent_it != optical_photon_map->end()){
                    // this is a new process! let's update our map
                    UpdatePhotonSummary(parent_id, track_id, g4_delta_time_step, creationProcessName, parent_it, new_process, aStep, g4_delta_distance);

                } else {
                    // need to add a new photon to our map
                    AddEntryToPhotonSummary(parent_id, track_id, g4_delta_distance, original_wavelength, primary_.GetTime() + pre_step_global_time + g4_delta_time_step, creationProcessName);
                }
            }
        }
        return false;
    }

    if(not (TrackParticles_ or DetailedPhotonTracking_ or TrackEnergyLosses_))
        return false;

    // position
    G4ThreeVector g4Position = aStep->GetPostStepPoint()->GetPosition();
    I3Position position(g4Position.x() / mm * I3Units::mm, g4Position.y() / mm * I3Units::mm, g4Position.z() / mm * I3Units::mm);

    // direction
    G4ThreeVector g4Direction = aStep->GetPostStepPoint()->GetMomentumDirection();
    I3Direction direction(g4Direction.x(), g4Direction.y(), g4Direction.z());

    // time
    double time = primary_.GetTime() + aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;

    // let's also grab parent id
    // if parent id == 0, that's our primary injected particle
    G4int parent_id = aStep->GetTrack()->GetParentID();
    G4int track_id = aStep->GetTrack()->GetTrackID();

    G4VProcess const * process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    std::string processName = (process) ? process->GetProcessName() : "Unknown";

    // now save to our MCTree!
    if(parent_id == 0) {
        // let's create and fill our I3Particle
        // since parent id = 0, we need to add daughter energy loss (aka processName)

        // process name -- use for parent id == 0!
        if(energyLossToI3ParticlePDGCode.find(processName) != energyLossToI3ParticlePDGCode.end()) {
            if(TrackEnergyLosses_) {
                I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(processName));
                I3Particle daughter(daughter_type);

                double edep = aStep->GetTotalEnergyDeposit() / electronvolt * I3Units::eV;
                daughter.SetEnergy(edep);
                daughter.SetPos(position);
                daughter.SetDir(direction);
                daughter.SetTime(time);

                if (DaughterParticleMap.find(1) == DaughterParticleMap.end()){
                    G4cout << "oops! no primary particle in the map!" << std::endl;
                } else {
                    I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(1), daughter); // append energy deposition to primary particle
                }
            }
        } else {
            G4cout << "oops! no conversion for " << processName << std::endl;
        }
    } else if (parent_id > 0) {
        // ok so we've created a new particle
        // if this is the first time we're seeing this particle -- add particle + energy loss
        // if we've already added this daughter particle -- only add energy loss

        if(DaughterParticleMap.find(track_id) == DaughterParticleMap.end()) {
            // we have not added the daugher...let's do it now
            I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(pdg);
            I3Particle daughter(daughter_type);
            daughter.SetEnergy(aStep->GetTrack()->GetVertexKineticEnergy() / electronvolt * I3Units::eV);
            daughter.SetPos(position);
            daughter.SetDir(direction);

            I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(parent_id) , daughter);

            // update map
            DaughterParticleMap[track_id] = daughter.GetID();
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
                daughter.SetDir(direction);
                if(DaughterParticleMap.find(track_id) == DaughterParticleMap.end()) {
                    G4cout << "oops! trying to save energy deposition type but DaughterParticleMap does not have track id!" << std::endl;
                } else {
                    I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id), daughter);
                }
            } else {
                G4cout << "oops! no conversion for " << processName << std::endl;
            }
        }
    }

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
