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
                                                                                                      {"Cerenkov", CCMMCPE::PhotonSource::Cherenkov},
                                                                                                      {"OpWLS", CCMMCPE::PhotonSource::OpWLS}};

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

    event_id = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    CCMMCPEMap = readout_->GetMCPESeries(event_id);
}

void G4CCMPMTSD::EndOfEvent(G4HCofThisEvent*) {
    int thread_id = G4Threading::G4GetThreadId();
    readout_->LogPMTResult(event_id, CCMMCPEMap, FullPhotonTracking_);
    Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4CCMPMTSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    // need to know if this is an optical photon
    if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

    G4VProcess const * process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
    std::string processName = (process) ? process->GetProcessName() : "Unknown";

    // We only want to save photons that are in the PMT
    // These processes can be registered with the SD even if they are not in the PMT
    if(processName == "Transportation" || processName == "OpBoundary" || processName == "SurfaceRefl")
        return false;

    // User replica number 1 since photocathode is a daughter volume to the pmt which was replicated
    G4int pmtNumber = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
    G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetTouchable()->GetVolume(0);
    std::string physVolName = static_cast<std::string>(physVol->GetName()); // the name is CoatedPMT_row_pmtNumber

    // let's convert physVolName to a row and pmt number to make a CCMPMTKey
    CCMPMTKey key;
    std::map<std::string, CCMPMTKey>::iterator vol_it = volumeToKey.find(physVolName);
    if (vol_it == volumeToKey.end()) {
        std::stringstream string_stream(physVolName);
        std::string segment;
        std::getline(string_stream, segment, '_');
        std::getline(string_stream, segment, '_');
        int row = std::stoi(segment);
        std::getline(string_stream, segment);
        int pmt_number = std::stoi(segment);

        key = CCMPMTKey(row, pmt_number);
        volumeToKey.insert(std::pair<std::string, CCMPMTKey>(physVolName, key));
    } else {
        key = vol_it->second;
    }
    // so now we need to save the info to a CCMMCPE
    // then we will save associated with correct PMT
    // we will be saving time, wavelength, positions, direction, and source

    // let's grab everything from our step
    //G4ThreeVector photonPosition = aStep->GetPostStepPoint()->GetPosition();
    //I3Position position(photonPosition.x() / mm * I3Units::mm, photonPosition.y() / mm * I3Units::mm, photonPosition.z() / mm * I3Units::mm);

    //G4ThreeVector photonDirection = aStep->GetPostStepPoint()->GetMomentumDirection();
    //I3Direction direction(photonDirection.x(), photonDirection.y(), photonDirection.z());

    G4double globalTime = aStep->GetPostStepPoint()->GetGlobalTime() / nanosecond * I3Units::nanosecond;
    //G4double localTime = aStep->GetPostStepPoint()->GetLocalTime() / nanosecond * I3Units::nanosecond;

    G4double photonEnergy = aStep->GetTrack()->GetTotalEnergy() / electronvolt;

    G4double photonWavelength = hc / photonEnergy * I3Units::nanometer;

    const G4VProcess* creationProcess = aStep->GetTrack()->GetCreatorProcess();
    std::string creationProcessName = "Unknown";
    if (creationProcess) {
        creationProcessName = static_cast<std::string>(creationProcess->GetProcessName());
    }
    G4int parent_id = aStep->GetTrack()->GetParentID();
    G4int track_id = aStep->GetTrack()->GetTrackID();

    // now save to CCMMCPE
    if (processNameToPhotonSource.find(creationProcessName) == processNameToPhotonSource.end()){
        std::cout << "oops!!! no photon process for " << creationProcessName << std::endl;
        return false;
    }

    // let's add this CCMPMTKey to CCMMCPEMap if it does not exist already
    if (CCMMCPEMap->find(key) == CCMMCPEMap->end()) {
        (*CCMMCPEMap)[key] = CCMMCPESeries();
    }

    // so we have the CCMPMTKey and we have the CCMMCPE, let's save!
    CCMMCPEMap->at(key).emplace_back(
        parent_id, track_id,
        std::vector<size_t>{}, // Number of photons per WLS, updated later
        WLSLocationSeries(), // WLS locations, updated later
        globalTime, photonWavelength,
        0.0, // Distance in UV, updated later
        0.0, // Original wavelength, updated later
        0.0, // Distance in visible, updated later
        processNameToPhotonSource.at(creationProcessName) // Source
    );

    // now kill photon after registering hit
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    return true;
}

