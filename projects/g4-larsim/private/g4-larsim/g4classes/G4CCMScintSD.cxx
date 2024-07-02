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

    // our scint SD is tracking energy deposited in the fiducial argon 
    // this will be used for voxelization
    
    // we don't care about optical photons for getting energy deposited in LAr
    // and if we don't care about PMTs, then we can kill any optical photon particle tracks
    // (this will make the simulation faster)
    if(aStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        if (!PMTSDStatus_){
            aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        }
        return false;
    }

    // now let's check energy deposited
    G4double edep = aStep->GetTotalEnergyDeposit() / electronvolt * I3Units::eV;
    G4double ekin = aStep->GetTrack()->GetKineticEnergy() / electronvolt * I3Units::eV;

    // now we want to grab energy deposited, location, direction, time, and process type to save to MCTree

    // position
    G4ThreeVector photonPosition = aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector prePos = aStep->GetPreStepPoint()->GetPosition();
    double delta_pos = std::sqrt(std::pow(photonPosition.x()/mm - prePos.x()/mm, 2) + std::pow(photonPosition.y()/mm - prePos.y()/mm, 2) + std::pow(photonPosition.z()/mm - prePos.z()/mm, 2));
    I3Position position(photonPosition.x() / mm * I3Units::mm, photonPosition.y() / mm * I3Units::mm, photonPosition.z() / mm * I3Units::mm);

    // direction
    G4ThreeVector photonDirection = aStep->GetPostStepPoint()->GetMomentumDirection();
    I3Direction direction(photonDirection.x(), photonDirection.y(), photonDirection.z());

    // time
    G4double photonTime = aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;

    // process type
    const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
    std::string creationProcessName = "Unknown";
    if (creationProcess) {
        creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
    }

    // let's also grab parent id
    // if parent id == 0, that's our primary injected particle
    G4int parent_id = aStep->GetTrack()->GetParentID();

    G4ParticleDefinition* fParticleDefinition = aStep->GetTrack()->GetDefinition();
    G4int pdg = fParticleDefinition->GetPDGEncoding();
    G4String particleName = fParticleDefinition->GetParticleName();
    
    //if (parent_id == 0){
    //    G4TrackVector* secTracks = G4EventManager::GetEventManager()->GetTrackingManager()->GimmeSecondaries();
    //    for (size_t i = 0; i < secTracks->size(); ++i) {
    //        G4Track* secondaryTrack = (*secTracks)[i];
    //        G4String secondaryName = secondaryTrack->GetDefinition()->GetParticleName();

    //        G4cout << "Secondary particle name: " << secondaryName << G4endl;
    //    }
    //}

    //if (creationProcessName == "Radioactivation"){
    //    // get name of parent particle 
    //    G4Track* parentTrack = aStep->GetTrack()->GetParentTrack();
    //    std::cout << particleName << " was created via " << creationProcessName << " and the parent particle = " << parentTrack->GetDefinition()->GetParticleName() << std::endl;
    //}

    // we can kill neutrinos 
    if (fParticleDefinition == G4NeutrinoE::NeutrinoE()){
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        return false;
    } 
    
    if(edep == 0.)
        return false;  // No edep so don't count as hit

    std::cout << "creation process name = " << creationProcessName << ", parent id = " << parent_id << ", track id = " << aStep->GetTrack()->GetTrackID() << ", name = " 
        << particleName << ", edep = "  << edep << ", and e kin = " << ekin << std::endl; 
   

    // now save to our MCTree!
    if (parent_id == 0){
        // let's create and fill our I3Particle
        // since parent id = 0, this will also be a primary
        I3Particle primary(primaryParticleType_);
        primary.SetEnergy(edep);
        primary.SetPos(position);
        primary.SetDir(direction);
 
        prev_parent_particle_id_ = primary.GetID(); 
        I3MCTreeUtils::AddPrimary(*mcTree, primary);
    }
    else if (parent_id > 0) {
        // this will be the daughter in our mcTree
        // let's grab the particle id of the last I3Particle in previous parent id
        I3ParticleID prev_parent_particle_id = prev_parent_particle_id_;
        
        // now let's get particle type for this daughter
        I3Particle::ParticleType daughter_type = static_cast<I3Particle::ParticleType>(pdg);
        I3Particle daughter(daughter_type);
        daughter.SetEnergy(edep);
        daughter.SetPos(position);
        daughter.SetDir(direction);

        I3MCTreeUtils::AppendChild(*mcTree, prev_parent_particle_id , daughter);
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
