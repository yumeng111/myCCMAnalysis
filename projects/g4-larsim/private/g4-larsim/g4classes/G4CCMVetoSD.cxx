
#include "g4-larsim/g4classes/G4CCMVetoSD.h"
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


G4CCMVetoSD::G4CCMVetoSD(G4String name) : G4VSensitiveDetector(name) {
    collectionName.insert("vetoScintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMVetoSD::Initialize(G4HCofThisEvent* hitsCE) {
    fScintCollection = new G4CCMScintHitsCollection(SensitiveDetectorName, collectionName[0]);

    if(fHitsCID < 0) {
        fHitsCID = G4SDManager::GetSDMpointer()->GetCollectionID(fScintCollection);
    }
    hitsCE->AddHitsCollection(fHitsCID, fScintCollection);

    event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    primary_ = readout_->GetPrimary(event_id);
    mcTree = readout_->GetMCTree(event_id);

    if (mcTree == nullptr) {
        mcTree = I3MCTreePtr(new I3MCTree());
    }

    DaughterParticleMap[1] = primary_.GetID();
    if(not I3MCTreeUtils::Has(*mcTree, primary_.GetID())) {
        I3MCTreeUtils::AddPrimary(*mcTree, primary_);
    }
}

void G4CCMVetoSD::EndOfEvent(G4HCofThisEvent*) {
    readout_->AddEntry(G4Threading::G4GetThreadId(), event_id, mcTree, photon_summary, optical_photon_map, FullPhotonTracking_);
    Reset();
}

void G4CCMVetoSD::AddEntryToPhotonSummary(int parent_id, int track_id, double g4_uv_distance, double g4_vis_distance,
                                           double calculated_uv_distance, double calculated_vis_distance,
                                           double g4_time, double calculated_time, std::string creationProcessName){
    // map does not have key -- let's add our PhotonSummary then update map
    size_t n_wls = 0;
    std::vector<size_t> n_photons_per_wls = {0};
    PhotonSummary this_photon_summary = PhotonSummary(g4_uv_distance, //g4_vis_distance,
                                                      //calculated_uv_distance, calculated_vis_distance,
                                                      g4_time, //calculated_time,
                                                      n_wls, n_photons_per_wls, WLSLocationSeries(),
                                                      processNameToPhotonSource.at(creationProcessName),
                                                      processNameToPhotonSource.at(creationProcessName),
                                                      processNameToPhotonSource.at(creationProcessName));
    photon_summary->push_back(this_photon_summary);
    optical_photon_map->insert(std::make_pair(track_id, photon_summary->size() - 1));
}

void G4CCMVetoSD::UpdatePhotonSummary(int parent_id, int track_id, double g4_uv_distance, double g4_vis_distance,
                                       double calculated_uv_distance, double calculated_vis_distance,
                                       double g4_time, double calculated_time, std::string creationProcessName,
                                       std::map<int, size_t>::iterator it, bool new_process, G4Step* aStep) {
    // map contains this photon
    // so we need to grab PhotonSummary and update it -- then update key
    PhotonSummary this_photon_summary = photon_summary->at(it->second);
    this_photon_summary.g4_distance_uv += g4_uv_distance;
    //this_photon_summary.g4_distance_visible += g4_vis_distance;
    //this_photon_summary.calculated_distance_uv += calculated_uv_distance;
    //this_photon_summary.calculated_distance_visible += calculated_vis_distance;
    this_photon_summary.g4_time += g4_time;
    //this_photon_summary.calculated_time += calculated_time;
    this_photon_summary.current_process = processNameToPhotonSource.at(creationProcessName);

    if (new_process){
        this_photon_summary.temp_parent = photon_summary->at(optical_photon_map->find(parent_id)->second).current_process;
    }

    if (creationProcessName == "OpWLS" and new_process){
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

    //if  (creationProcessName == "OpWLS"){
    //    std::cout << "optical photon parent id = " << parent_id << ", track id = " << track_id
    //        << ", creation process = " << creationProcessName
    //        //<< ", photon source = " << photonSourceToProcessName.at(photon_summary->at((*optical_photon_map)[track_id]).photon_source)
    //        //<< ", temp parent = " << photonSourceToProcessName.at(photon_summary->at((*optical_photon_map)[track_id]).temp_parent)
    //        << ", new_process = " << new_process
    //        << ", n_wls = " << photon_summary->at((*optical_photon_map)[track_id]).n_wls
    //        << ", photons produced per wls = " << photon_summary->at((*optical_photon_map)[track_id]).n_photons_per_wls
    //        << ", pre step vol = " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName()
    //        << ", photon summary wls loc = " ;
    //        WLSLocationSeries this_wls_loc = photon_summary->at((*optical_photon_map)[track_id]).wls_loc;
    //        for (size_t i = 0; i < this_wls_loc.size(); i++){
    //            std::cout << wlsLocationToProcessName.at(this_wls_loc.at(i).wls_loc) << ", ";
    //        }
    //        std::cout << std::endl;
    //}

}

double G4CCMVetoSD::InterpolateRindex(double wavelength){
    // this function takes a wavelength
    // and interpolates to find closed rindex

    auto it = std::min_element(rindex_wavelength.begin(), rindex_wavelength.end(), [wavelength](double a, double b) {
                                   return std::abs(a - wavelength) < std::abs(b - wavelength); });

    size_t closest_idx = std::distance(rindex_wavelength.begin(), it);

    size_t upper_idx;
    size_t lower_idx;

    if (rindex_wavelength.at(closest_idx) >= wavelength){
        upper_idx = closest_idx;
        if (upper_idx > 0){
            lower_idx = upper_idx - 1;
        } else {
            lower_idx = upper_idx;
        }

    } else {
        lower_idx = closest_idx;
        if (lower_idx < (rindex_wavelength.size() - 1)){
            upper_idx = lower_idx + 1;
        } else {
            upper_idx = lower_idx;
        }
    }

    double wavelength_above = rindex_wavelength.at(upper_idx);
    double rindex_above = rindex.at(upper_idx);

    double wavelength_below = rindex_wavelength.at(lower_idx);
    double rindex_below = rindex.at(lower_idx);


    // now interpolate!
    double interpolated_rindex = rindex_below + (wavelength - wavelength_below) * ((rindex_above - rindex_below) / (wavelength_above - wavelength_below));

    return interpolated_rindex;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CCMVetoSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    // note -- this chunk of code resets global time to 0 in the case of radioactive decays
    // very important for retaining time structure of scintillation photons!!!
    G4Track* track = aStep->GetTrack();

    // get pdg code
    G4int pdg = particle_definition->GetPDGEncoding();

    if(particle_definition == G4OpticalPhoton::OpticalPhotonDefinition()) {
        return false;
    }

    if(track->GetTrackStatus() == fStopAndKill) {
        return false;
    }

    G4int parent_id = aStep->GetTrack()->GetParentID();
    G4int track_id = aStep->GetTrack()->GetTrackID();

    // now we want to grab energy deposited, location, direction, time, and process type to save to MCTree

    bool source_from_tree = SaveEnergyLossesTree_ or SaveEnergyLosses_ and tree_tracker->GetTrackEnergyLosses();

    if(source_from_tree) {
        auto tree_tracker->DaughterParticleMap::const_iterator it = DaughterParticleMap.find(track_id);
        if(it != DaughterParticleMap.cend()) {
            std::vector<I3Particle> daughters = I3MCTreeUtils::GetDaughters(*(tree_tracker->mcTree), *it);
            I3Particle energy_loss = std::max_element(daughters.begin(), duaghters.end(), [](I3Particle const & p0, I3Particle const & p1) -> bool {return p0.GetID() < p1.GetID();})*;
            if(SaveEnergyLossesTree_) {
                energy_loss_ids.push_back(energy_loss.GetID());
            }
            if(SaveEnergyLosses_) {
                energy_losses->push_back(energy_loss);
            }
        }
        return false;
    }

    // now let's check energy deposited
    G4double edep = aStep->GetTotalEnergyDeposit() / electronvolt * I3Units::eV;
    G4double ekin = aStep->GetTrack()->GetKineticEnergy() / electronvolt * I3Units::eV;

    // position
    G4ThreeVector photonPosition = aStep->GetPostStepPoint()->GetPosition();
    I3Position position(photonPosition.x() / mm * I3Units::mm, photonPosition.y() / mm * I3Units::mm, photonPosition.z() / mm * I3Units::mm);

    // direction
    G4ThreeVector photonDirection = aStep->GetPostStepPoint()->GetMomentumDirection();
    I3Direction direction(photonDirection.x(), photonDirection.y(), photonDirection.z());

    // time
    G4double photonTime = aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;

    // creation process -- use for parent id > 0
    const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
    std::string creationProcessName = "Unknown";
    if (creationProcess) {
        creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
    }

    // process name -- use for parent id == 0!
    std::string processName = "Unknown";
    const G4VProcess* currentProcess = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    if (currentProcess) {
        processName = static_cast<std::string>(currentProcess->GetProcessName());
    }

    // let's also grab parent id
    // if parent id == 0, that's our primary injected particle
    G4int parent_id = aStep->GetTrack()->GetParentID();
    G4int track_id = aStep->GetTrack()->GetTrackID();

    // get name and pdg code
    G4ParticleDefinition* fParticleDefinition = aStep->GetTrack()->GetDefinition();
    G4int pdg = fParticleDefinition->GetPDGEncoding();
    G4String particleName = fParticleDefinition->GetParticleName();

    // kill neutrinos
    if (fParticleDefinition == G4NeutrinoE::NeutrinoE()){
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        return false;
    }

    //std::cout << "creation process name = " << creationProcessName << ", processName = " << processName << ", parent id = " << parent_id
    //    << ", track id = " << track_id << ", name = " << particleName << ", edep = "  << edep << ", e kin = " << ekin << ", and time = " << photonTime << std::endl;

    // now save to our MCTree!
    if (parent_id == 0) {
        //std::cout << "energy deposition name = " << processName << ", parent id = " << parent_id << ", track id = " << track_id
        //    << ", particle name = " << particleName << ", edep = "  << edep << std::endl;
        // let's create and fill our I3Particle
        // since parent id = 0, we need to add daughter energy loss (aka processName)
        if (energyLossToI3ParticlePDGCode.find(processName) != energyLossToI3ParticlePDGCode.end()){
            I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(processName));
            I3Particle daughter(daughter_type);
            daughter.SetEnergy(edep);
            daughter.SetPos(position);
            daughter.SetDir(direction);

            if (DaughterParticleMap.find(1) == DaughterParticleMap.end()){
                std::cout << "oops! no primary particle in the map!" << std::endl;
            } else {
                I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(1), daughter); // append energy deposition to primary particle
            }
        } else {
            std::cout << "oops! no conversion for " << processName << std::endl;
        }
    } else if (parent_id > 0) {
        //std::cout << "energy deposition name = " << creationProcessName << ", parent id = " << parent_id << ", track id = " << track_id
        //    << ", particle name = " << particleName << ", edep = "  << edep << std::endl;
        // ok so we've created a new particle
        // if this is the first time we're seeing this particle -- add particle + energy loss
        // if we've already added this daughter particle -- only add energy loss

        if (DaughterParticleMap.find(track_id) == DaughterParticleMap.end()) {
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
        if (energyLossToI3ParticlePDGCode.find(creationProcessName) != energyLossToI3ParticlePDGCode.end()){
            // now add energy loss
            I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(creationProcessName));
            I3Particle daughter(daughter_type);
            daughter.SetEnergy(edep);
            daughter.SetPos(position);
            daughter.SetDir(direction);
            if (DaughterParticleMap.find(track_id) == DaughterParticleMap.end()){
                std::cout << "oops! trying to save energy deposition type but DaughterParticleMap does not have track id!" << std::endl;
            } else {
                I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id) , daughter);
            }
        } else {
            std::cout << "oops! no conversion for " << creationProcessName << std::endl;
        }

    }

    // now back to scint hit things
    G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
    auto theTouchable = (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
    G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();

    G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

    // Get the average position of the hit
    G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
    pos /= 2.;

    auto scintHit = new G4CCMScintHit(thePrePV);

    scintHit->SetEdep(edep);
    scintHit->SetPos(pos);

    fScintCollection->insert(scintHit);

    return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
