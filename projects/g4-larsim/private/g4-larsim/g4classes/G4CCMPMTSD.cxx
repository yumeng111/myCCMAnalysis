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
/// \file optical/LXe/src/G4CCMPMTSD.cc
/// \brief Implementation of the G4CCMPMTSD class
//
//
#include "g4-larsim/g4classes/G4CCMPMTHit.h"
#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
const std::unordered_map<std::string, CCMMCPE::PhotonSource> G4CCMPMTSD::processNameToPhotonSource = {{"Unknown", CCMMCPE::PhotonSource::Unknown},
                                                                                                      {"Scintillation", CCMMCPE::PhotonSource::Scintillation},
                                                                                                      {"Cerenkov", CCMMCPE::PhotonSource::Cerenkov}};

G4CCMPMTSD::G4CCMPMTSD(G4String name) : G4VSensitiveDetector(name) {
    collectionName.insert("pmtHitCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMPMTSD::~G4CCMPMTSD() {
    delete fPMTPositionsX;
    delete fPMTPositionsY;
    delete fPMTPositionsZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMPMTSD::SetPmtPositions(const std::vector<G4ThreeVector>& positions) {
    for(size_t i = 0; i < positions.size(); ++i) {
        if(fPMTPositionsX)
            fPMTPositionsX->push_back(positions[i].x());
        if(fPMTPositionsY)
            fPMTPositionsY->push_back(positions[i].y());
        if(fPMTPositionsZ)
        fPMTPositionsZ->push_back(positions[i].z());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMPMTSD::Initialize(G4HCofThisEvent* hitsCE) {
    fPMTHitCollection = new G4CCMPMTHitsCollection(SensitiveDetectorName, collectionName[0]);

    if(fHitCID < 0) {
        fHitCID = G4SDManager::GetSDMpointer()->GetCollectionID(fPMTHitCollection);
    }
    hitsCE->AddHitsCollection(fHitCID, fPMTHitCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CCMPMTSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    
    // need to know if this is an optical photon
    if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

    // User replica number 1 since photocathode is a daughter volume to the pmt which was replicated
    G4int pmtNumber = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
    //G4VPhysicalVolume* physVol = aStep->GetPostStepPoint()->GetTouchable()->GetVolume(1);
    G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(1);
    std::string physVolName = static_cast<std::string>(physVol->GetName()); // the name is row_pmtNumber 

    // let's convert physVolName to a row and pmt number to make a CCMPMTKey
    std::stringstream string_stream(physVolName);
    std::string segment;
    std::getline(string_stream, segment, '_');
    int row = std::stoi(segment);
    std::getline(string_stream, segment);
    int pmt_number = std::stoi(segment);

    CCMPMTKey key = CCMPMTKey(row, pmt_number);
    // so now we need to save the info to a CCMMCPE
    // then we will save associated with correct PMT
    // we will be saving time, wavelength, positions, direction, and source
    
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


    // Find the correct hit collection
    size_t n = fPMTHitCollection->entries();
    G4CCMPMTHit* hit = nullptr;
    for(size_t i = 0; i < n; ++i) {
        if((*fPMTHitCollection)[i]->GetPMTNumber() == pmtNumber) {
            hit = (*fPMTHitCollection)[i];
            break;
        }
    }

    if(hit == nullptr) { // this pmt wasn't previously hit in this event
        hit = new G4CCMPMTHit();  // so create new hit
        hit->SetPMTNumber(pmtNumber);
        hit->SetPMTPhysVol(physVol);
        fPMTHitCollection->insert(hit);
        hit->SetPMTPos((*fPMTPositionsX)[pmtNumber], (*fPMTPositionsY)[pmtNumber], (*fPMTPositionsZ)[pmtNumber]);
    }

    hit->IncPhotonCount();  // increment hit for the selected pmt
    hit->SetDrawit(true);
    
    // let's add this CCMPMTKey to CCMMCPEMap if it does not exist already
    if (CCMMCPEMap->find(key) == CCMMCPEMap->end()) {
        (*CCMMCPEMap)[key] = CCMMCPESeries ();
    }

    // so we have the CCMPMTKey and we have the CCMMCPE, let's save!
    for (auto it = CCMMCPEMap->begin(); it != CCMMCPEMap->end(); ++it) {
        CCMPMTKey mc_pe_map_key = it->first;
        if (mc_pe_map_key == key){
            CCMMCPEMap->at(mc_pe_map_key).push_back(this_mc_pe);
        }
    }   

    // now kill photon after registering hit
    aStep->GetTrack()->SetTrackStatus(fStopAndKill); 

    return true;
}

