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

#include "G4ios.hh"
#include "G4LogicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"
#include "G4VTouchable.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::unordered_map<std::string, CCMMCPE::PhotonSource> G4CCMScintSD::processNameToPhotonSource = {{"Unknown", CCMMCPE::PhotonSource::Unknown},
                                                                                                      {"Scintillation", CCMMCPE::PhotonSource::Scintillation},
                                                                                                      {"Cerenkov", CCMMCPE::PhotonSource::Cerenkov}};

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

    // our scint SD is tracking scintillation photons produced in the fiducial argon
    // this will be used for voxelization
    
    // let's make sure this is an optical photon
    if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;
    
    // if we don't care about PMTs, then we can kill any non-primary particle tracks
    if(aStep->GetTrack()->GetParentID() > 1 and !PMTSDStatus_) {
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);
        return false;
    }

    // regardless of saving PMT hits or not, we only want to save ParentID = 1 in the LAr
    if(aStep->GetTrack()->GetParentID() != 1) {
        return false;
    }

    // let's also do some scint hits stuff idk

    G4double edep = aStep->GetTotalEnergyDeposit();
    if(edep == 0.)
        return false;  // No edep so don't count as hit
    
    // let's grab everything from our step
    G4ThreeVector photonPosition = aStep->GetPostStepPoint()->GetPosition();
    I3Position position(photonPosition.x(), photonPosition.y(), photonPosition.z());

    G4ThreeVector photonDirection = aStep->GetPostStepPoint()->GetMomentumDirection();
    I3Direction direction(photonDirection.x(), photonDirection.y(), photonDirection.z());

    G4double photonTime = aStep->GetPostStepPoint()->GetGlobalTime();
    G4double photonEnergy = aStep->GetTrack()->GetTotalEnergy();
    G4double photonWavelength = h_Planck * c_light / photonEnergy;
    const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
    std::string creationProcessName = "Unknown";
    if (creationProcess) {
        creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
    }

    // now save to CCMMCPE!
    CCMMCPE this_mc_pe = CCMMCPE(photonTime, photonWavelength, position, direction, processNameToPhotonSource.at(creationProcessName));

    // now push back
    CCMMCPEList->push_back(this_mc_pe);

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
