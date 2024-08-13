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
/// \file optical/LXe/src/LXeSteppingAction.cc
/// \brief Implementation of the LXeSteppingAction class
//
//
#include "G4CCMSteppingAction.h"
#include "G4CCMEventAction.h"
#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "G4CCMSteppingMessenger.h"
#include "G4CCMTrajectory.h"
#include "G4CCMUserTrackInformation.h"

#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include <G4ParticleTypes.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMSteppingAction::G4CCMSteppingAction(G4CCMEventAction* ea) : fEventAction(ea) {
    fSteppingMessenger = new G4CCMSteppingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMSteppingAction::~G4CCMSteppingAction() { delete fSteppingMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMSteppingAction::UserSteppingAction(const G4Step* theStep) {
    G4Track* theTrack = theStep->GetTrack();
    const G4ParticleDefinition* part = theTrack->GetDefinition();
    G4String particleName = part->GetParticleName();
    G4int pdg = part->GetPDGEncoding();

    //if(theTrack->GetCurrentStepNumber() == 1)
    //    fExpectedNextStatus = Undefined;

    auto trackInformation = static_cast<G4CCMUserTrackInformation*>(theTrack->GetUserInformation());

    G4StepPoint* thePrePoint    = theStep->GetPreStepPoint();
    G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
    std::string thePrePVName = static_cast<std::string>(thePrePV->GetName());

    G4StepPoint* thePostPoint    = theStep->GetPostStepPoint();
    G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
    std::string thePostPVName = static_cast<std::string>(thePostPV->GetName());

    std::string desired_str ("PMT");
    std::string desired_str2 ("TPB");

    if(theTrack->GetParentID() == 0) {
        // This is a primary track
        trackInformation->SetTrackStatusFlags(primaryLArDeposition);
    }

    if(nullptr == thePostPV) {  // out of world
        fExpectedNextStatus = Undefined;
        return;
    }
    
    std::cout << particleName << " in " << thePostPVName << std::endl;
    
    //kill neutrinos 
    if (theTrack->GetDefinition() == G4NeutrinoE::NeutrinoE()){
        theTrack->SetTrackStatus(fStopAndKill);
    } 

    // Optical photon only
    if(pdg == -22) {
        // grab photons on PMT and then kill! 
        if (thePrePVName.find(desired_str) != std::string::npos) {
            trackInformation->AddTrackStatusFlag(hitPMT);
            theTrack->SetTrackStatus(fStopAndKill);
            
        } 
        if (thePostPVName.find(desired_str) != std::string::npos) {
            trackInformation->AddTrackStatusFlag(hitPMT);
            theTrack->SetTrackStatus(fStopAndKill);
        } 
        if (thePrePVName.find(desired_str2) != std::string::npos) {
            trackInformation->AddTrackStatusFlag(hitPMT);
            theTrack->SetTrackStatus(fStopAndKill);
            
        } 
        if (thePostPVName.find(desired_str2) != std::string::npos) {
            trackInformation->AddTrackStatusFlag(hitPMT);
            theTrack->SetTrackStatus(fStopAndKill);
        } 
        //if (thePrePV->GetName() == "photocath" or thePostPV->GetName() == "photocath"){
        //  // draw anything in our pmts
        //  //trackInformation->AddTrackStatusFlag(hitPMT);
        //}
        if(thePrePVName == "expHall") {
            // Kill photons entering expHall from something other than Slab
            theTrack->SetTrackStatus(fStopAndKill);
        }
    }
}

