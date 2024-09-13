//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file optical/LXe/src/LXeScintSD.cc
/// \brief Implementation of the LXeScintSD class
//
//
#include "g4-larsim/g4classes/G4CCMScintSD.h"
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
const std::unordered_map<std::string, int> G4CCMScintSD::energyLossToI3ParticlePDGCode = {{"phot", 2000000001}, {"compt", 2000000002}, {"conv", 2000000003},
                                                                                          {"Rayl", 2000000004}, {"msc", 2000000005}, {"eIoni", 2000000006},
                                                                                          {"eBrem", 2000000007}, {"ePairProd", 2000000008}, {"CoulombScat", 2000000009},
                                                                                          {"annihil", 2000000010}, {"Cerenkov", 2000000011}, {"Radioactivation", 2000000012},
                                                                                          {"Scintillation", 2000000013}, {"OpWLS", 2000000014}};

const std::unordered_map<std::string, PhotonSummary::PhotonSource> G4CCMScintSD::processNameToPhotonSource = {{"Unknown", PhotonSummary::PhotonSource::Unknown},
                                                                                                              {"Scintillation", PhotonSummary::PhotonSource::Scintillation},
                                                                                                              {"Cerenkov", PhotonSummary::PhotonSource::Cerenkov}};

G4CCMScintSD::G4CCMScintSD(G4String name) : G4VSensitiveDetector(name) {
    collectionName.insert("scintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMScintSD::Initialize(G4HCofThisEvent* hitsCE) {
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

void G4CCMScintSD::EndOfEvent(G4HCofThisEvent*) {
    readout_->AddEntry(G4Threading::G4GetThreadId(), event_id, mcTree, photon_summary, optical_photon_map);
    Reset();
}

void G4CCMScintSD::AddEntryToPhotonSummary(int parent_id, int track_id, double g4_uv_distance, double g4_vis_distance,
                                           double calculated_uv_distance, double calculated_vis_distance,
                                           double g4_time, double calculated_time, std::string creationProcessName){
    // map does not have key -- let's add our PhotonSummary then update map
    size_t n_wls = 0;
    if (creationProcessName == "OpWLS"){
        n_wls = 1;
    }
    PhotonSummary this_photon_summary = PhotonSummary(g4_uv_distance, g4_vis_distance,
                                                      calculated_uv_distance, calculated_vis_distance,
                                                      g4_time, calculated_time, n_wls, processNameToPhotonSource.at(creationProcessName));
    photon_summary->push_back(this_photon_summary);
    //std::cout << "adding entry to optical photon map for photon at track id = " << track_id << ", distance uv = " << this_photon_summary.distance_uv
    //          << ", distance vis = " << this_photon_summary.distance_visible <<  ", and n_wls = " << this_photon_summary.n_wls << std::endl;
    //std::cout << "" << std::endl;
    optical_photon_map->insert(std::make_pair(track_id, photon_summary->size() - 1));
}

void G4CCMScintSD::UpdatePhotonSummary(int parent_id, int track_id, double g4_uv_distance, double g4_vis_distance,
                                       double calculated_uv_distance, double calculated_vis_distance,
                                       double g4_time, double calculated_time, std::string creationProcessName,
                                       std::map<int, size_t>::iterator it, bool new_process) {
    // map contains this photon
    // so we need to grab PhotonSummary and update it -- then update key
    PhotonSummary this_photon_summary = photon_summary->at(it->second);
    this_photon_summary.g4_distance_uv += g4_uv_distance;
    this_photon_summary.g4_distance_visible += g4_vis_distance;
    this_photon_summary.calculated_distance_uv += calculated_uv_distance;
    this_photon_summary.calculated_distance_visible += calculated_vis_distance;
    this_photon_summary.g4_time += g4_time;
    this_photon_summary.calculated_time += calculated_time;

    if (creationProcessName == "OpWLS" and new_process){
        this_photon_summary.n_wls += 1;
    }

    // now update the photon_summary, delete from map, and add new entry to the map
    size_t pos = it->second;
    photon_summary->at(pos) = this_photon_summary;
    //std::cout << "updating optical photon map for photon at track id = " << track_id << ", distance uv = " << this_photon_summary.distance_uv
    //          << ", distance vis = " << this_photon_summary.distance_visible <<  ", and n_wls = " << this_photon_summary.n_wls << std::endl;
    //std::cout << "" << std::endl;
    //optical_photon_map->erase(it);
    optical_photon_map->insert(std::make_pair(track_id, pos));
}

double G4CCMScintSD::InterpolateRindex(double wavelength){
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

G4bool G4CCMScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
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
        //std::cout << "process name = " << processName << std::endl;
        if(processName == "Radioactivation") {
            // G4String parentName = track->GetDefinition()->GetParticleName();
            // std::cout << "for parent particle = " << parentName << ", secondaries = " << std::endl;

            // Get the list of secondaries
            const G4TrackVector* secondaries = aStep->GetSecondary();
            // Modify the start time of each secondary particle
            for (size_t i = 0; i < secondaries->size(); ++i) {
                G4Track* secondary = const_cast<G4Track*>(secondaries->at(i));
                //std::cout << secondary->GetDefinition()->GetParticleName() << std::endl;
                secondary->SetGlobalTime(0.);
            }
        }
    }

    // let's do a check on the time cut
    if (TimeCut_) {
        // check time ... doesnt matter what type of particle it is
        G4double time = aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
        if (time > 200.0){
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return false;
        }
    }

    // ok back to SD logic
    // our scint SD is tracking energy deposited in the fiducial argon
    // we don't care about optical photons for getting energy deposited in LAr
    // and if we don't care about PMTs, then we can kill any optical photon particle tracks
    // (this will make the simulation faster)
    if(aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        if (!PMTSDStatus_){
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
            return false;
        }

        // check if we want to kill cerenkov photons
        if (KillCherenkov_){
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creationProcessName = "Unknown";
            if (creationProcess) {
                creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
            }
            if (creationProcessName == "Cerenkov"){
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return false;
            }
        }

        // ok if we survived all checks, add an entry to our photon summary
        // we need parent id, track id, distance travelled, and wavelength
        G4int parent_id = aStep->GetTrack()->GetParentID();
        G4int track_id = aStep->GetTrack()->GetTrackID();

        double g4_delta_distance = (aStep->GetStepLength() / mm) * I3Units::mm;
        double pre_step_global_time = aStep->GetPreStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
        double g4_delta_time_step = aStep->GetDeltaTime() / nanosecond * I3Units::nanosecond;

        double wavelength = hc / (aStep->GetTrack()->GetTotalEnergy() / electronvolt); // units of nanometer
        double interpolated_rindx = InterpolateRindex(wavelength);

        // based on the wavelength, let's classify as uv or vis
        double g4_vis_distance = 0.0;
        double g4_uv_distance = 0.0;
        double calculated_vis_distance = 0.0;
        double calculated_uv_distance = 0.0;

        if (wavelength <= 325.0){
            g4_uv_distance = g4_delta_distance;
            calculated_uv_distance = (c_mm_per_nsec * g4_delta_time_step) / interpolated_rindx;
        } else {
            g4_vis_distance = g4_delta_distance;
            calculated_vis_distance = (c_mm_per_nsec * g4_delta_time_step) / interpolated_rindx;
        }

        // let's also calculate the travel time based on the g4 distance
        double calculated_delta_time_step = (g4_delta_distance * interpolated_rindx) / c_mm_per_nsec;

        // let's also check if this photon got wls
        const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
        std::string creationProcessName = "Unknown";
        if (creationProcess) {
            creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
        }

        //std::cout << "optical photon parent id = " << parent_id << ", track id = " << track_id << ", step lenght = " << aStep->GetStepLength()
        //    << ", delta local time = " << aStep->GetPostStepPoint()->GetLocalTime() - aStep->GetPreStepPoint()->GetLocalTime()
        //    << ", and calculated delta time = " << (uv_distance / I3Units::cm) / (c_cm_per_nsec / uv_index_of_refraction) + (vis_distance / I3Units::cm)/ (c_cm_per_nsec / vis_index_of_refraction)
        //    << std::endl;

        // ok now let's save
        // check to see if the parent id is in our map
        bool new_process = true;
        std::map<int, size_t>::iterator it = optical_photon_map->find(parent_id);
        if (it != optical_photon_map->end()) {
            // update our optical photon map
            UpdatePhotonSummary(parent_id, track_id, g4_uv_distance, g4_vis_distance,
                                calculated_uv_distance, calculated_vis_distance,
                                g4_delta_time_step, calculated_delta_time_step, creationProcessName, it, new_process);
        } else {
            // check if this track id is in our map
            std::map<int, size_t>::iterator track_it = optical_photon_map->find(track_id);
            bool new_process = false;
            if (track_it != optical_photon_map->end()){
                // ok so this photon is in our map, let's just update
                UpdatePhotonSummary(parent_id, track_id, g4_uv_distance, g4_vis_distance,
                                    calculated_uv_distance, calculated_vis_distance,
                                    g4_delta_time_step, calculated_delta_time_step, creationProcessName, track_it, new_process);
            } else {
                // need to add a new photon to our map
                AddEntryToPhotonSummary(parent_id, track_id, g4_uv_distance, g4_vis_distance, calculated_uv_distance, calculated_vis_distance,
                                        primary_.GetTime() + pre_step_global_time + g4_delta_time_step, primary_.GetTime() + pre_step_global_time + calculated_delta_time_step, creationProcessName);
            }
        }

        //std::cout << "new process = " << new_process << std::endl;


        return false;
    }

    // now we want to grab energy deposited, location, direction, time, and process type to save to MCTree

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
