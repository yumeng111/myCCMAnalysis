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
const std::unordered_map<std::string, PhotonSummary::CreationProcess> G4CCMScintSD::processNameToCreationProcess = {{"Unknown", PhotonSummary::CreationProcess::Unknown},
                                                                                                                    {"Scintillation", PhotonSummary::CreationProcess::Scintillation},
                                                                                                                    {"Cerenkov", PhotonSummary::CreationProcess::Cerenkov},
                                                                                                                    {"OpWLS", PhotonSummary::CreationProcess::OpWLS}};

const std::unordered_map<std::string, int> G4CCMScintSD::energyLossToI3ParticlePDGCode = {{"phot", 2000000001}, {"compt", 2000000002}, {"conv", 2000000003},
                                                                                          {"Rayl", 2000000004}, {"msc", 2000000005}, {"eIoni", 2000000006},
                                                                                          {"eBrem", 2000000007}, {"ePairProd", 2000000008}, {"CoulombScat", 2000000009},
                                                                                          {"annihil", 2000000010}, {"Cerenkov", 2000000011}, {"Radioactivation", 2000000012}};

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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CCMScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

    // note -- this chunk of code resets global time to 0 in the case of radioactive decays
    // very important for retaining time structure of scintillation photons!!!
    G4Track* track = aStep->GetTrack();

    // Check if the particle has decayed
    if (track->GetTrackStatus() == fStopAndKill) {
        // Check if it's a primary particle (ParentID == 0)
        if (track->GetParentID() == 0) {
            // Get the list of secondaries
            G4String parentName = track->GetDefinition()->GetParticleName();
            //std::cout << "for parent particle = " << parentName << ", secondaries = " << std::endl;
            const G4TrackVector* secondaries = aStep->GetSecondary();
            // Modify the start time of each secondary particle
            for (size_t i = 0; i < secondaries->size(); ++i) {
                G4Track* secondary = const_cast<G4Track*>(secondaries->at(i));
                //std::cout << secondary->GetDefinition()->GetParticleName() << std::endl;
                secondary->SetGlobalTime(0.);
            }
        }
    }
    // ok back to SD logic

    // our scint SD is tracking energy deposited in the fiducial argon 
    // this will be used for voxelization
    
    // we don't care about optical photons for getting energy deposited in LAr
    // and if we don't care about PMTs, then we can kill any optical photon particle tracks
    // (this will make the simulation faster)
    if(aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        if (!PMTSDStatus_){
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        } else {
            // ok add an entry to our photon summary
            // we need parent id, track id, time, position, and creation process
            G4int parent_id = aStep->GetTrack()->GetParentID();
            G4int track_id = aStep->GetTrack()->GetTrackID();
            double time = static_cast<double>(aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond) * I3Units::nanosecond;
            G4ThreeVector photonPosition = aStep->GetPostStepPoint()->GetPosition();
            I3Position position(photonPosition.x() / mm * I3Units::mm, photonPosition.y() / mm * I3Units::mm, photonPosition.z() / mm * I3Units::mm);
            const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
            std::string creationProcessName = "Unknown";
            if (creationProcess) {
                creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
            }
            // ok now ready to make a PhotonSummary and then save
            PhotonSummary this_photon_summary = PhotonSummary(parent_id, track_id, time, position, processNameToCreationProcess.at(creationProcessName));
            // and now save
            photon_summary->push_back(this_photon_summary);

        }
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
    
    // do not add entry to MCTree for no energy deposition
    if(edep == 0.){
        return false; 
    }
    
    //std::cout << "creation process name = " << creationProcessName << ", processName = " << processName << ", parent id = " << parent_id
    //    << ", track id = " << track_id << ", name = " << particleName << ", edep = "  << edep << ", e kin = " << ekin << ", and time = " << photonTime << std::endl; 

    // now save to our MCTree!
    if (parent_id == 0){
        std::cout << "energy deposition name = " << processName << ", parent id = " << parent_id << ", track id = " << track_id
            << ", particle name = " << particleName << ", edep = "  << edep << std::endl; 
        // let's create and fill our I3Particle
        // since parent id = 0, we need to add daughter energy loss (aka processName) 
        I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(processName));
        I3Particle daughter(daughter_type);
        daughter.SetEnergy(edep);
        daughter.SetPos(position);
        daughter.SetDir(direction);

        I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(1), daughter); // append energy deposition to primary particle
    }
    else if (parent_id > 0) {
        std::cout << "energy deposition name = " << creationProcessName << ", parent id = " << parent_id << ", track id = " << track_id
            << ", particle name = " << particleName << ", edep = "  << edep << std::endl; 
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

        // now add energy loss 
        I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(energyLossToI3ParticlePDGCode.at(creationProcessName));
        I3Particle daughter(daughter_type);
        daughter.SetEnergy(edep);
        daughter.SetPos(position);
        daughter.SetDir(direction);
        I3MCTreeUtils::AppendChild(*mcTree, DaughterParticleMap.at(track_id) , daughter);
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
