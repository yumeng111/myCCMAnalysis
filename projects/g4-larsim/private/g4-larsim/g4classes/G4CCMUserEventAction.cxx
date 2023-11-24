/*
   Event Action for CCM simulation.

   This code offers the ability to perform actions at the start/end of events. 
   Currently handles the output to ROOT, using a file defined in the BeginOfEventAction, 
   trigger information obtained in the EndOfEventAction, and defines a method for adding
   a hit (optical photon to PMT).
*/
#include <sstream>
#include <iostream>

#include "g4-larsim/g4classes/G4CCMUserEventAction.hh"
#include "g4-larsim/g4classes/G4CCMUserRunAction.hh"
#include "g4-larsim/g4classes/G4CCMUserPrimaryGenerator.hh"
#include "g4-larsim/g4classes/G4CCMUserDetectorConstruction.hh"

#include "CCMAnalysis/CCMDataStructures/MCTruth.h"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

G4CCMUserEventAction::G4CCMUserEventAction(G4CCMUserRunAction* runaction)
    : G4UserEventAction(), fRunAction(runaction),
    mctruth(new MCTruth()), rootIO(std::make_shared<CCMRootIO>()) {
    rootSet = false;
}

G4CCMUserEventAction::~G4CCMUserEventAction() {
    // clear data and root file structures
    delete mctruth;
    rootIO->Close();
}

// Within this function define things to be done at the start of an event (example: the debugging line currently commented out).
// The code currently does not access the detector or the primary generator, so information produced there is not accessible.
void G4CCMUserEventAction::BeginOfEventAction(const G4Event* event) {
    // get the event ID, reset the rootIO for that event, and determine if the root filename is already set.
    const G4int eventID = event->GetEventID();

    mctruth->Reset(eventID);
    rootSet = fRunAction->GetRootSet();

    // if rootSet is false, pull the root file name and use it to open a new root file and rootIO structure
    if (!rootSet){
        const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        G4String filename = detector->GetRootFile();
        std::ostringstream oss;
        oss << filename << "thread" << eventID << ".root";
        filename = oss.str();

        rootIO->SetParameter("SaveBranches","mcTruth");
        rootIO->SetOutFileName(filename);
        rootIO->SetupOutputFile();
        fRunAction->SetRootSet(true);
    }
}

// Same as BeginOfEventAction, but for end of event.
void G4CCMUserEventAction::EndOfEventAction(G4Event const * event) {
    const G4int eventID = event->GetEventID();

    // pull the initial conditions from the primary Generator for usage in ROOT output.
    const G4CCMUserPrimaryGenerator* generator = static_cast<const G4CCMUserPrimaryGenerator*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    G4double xpos = generator->GetXpos();
    G4double ypos = generator->GetYpos();
    G4double zpos = generator->GetZpos();
    G4double momx = generator->GetXmom();
    G4double momy = generator->GetYmom();
    G4double momz = generator->GetZmom();
    G4double partEneg = generator->GetPartEneg();
    G4String partName = generator->GetPartName();

    // test output string for matching to primary generator notes
    // output the initial conditions to the rootIO and write for the Trigger.
    mctruth->SetParticlePosition(xpos,ypos,zpos);
    mctruth->SetParticleMomentum(momx,momy,momz);
    mctruth->SetParticleID(eventID);
    mctruth->SetTotalEnergyDeposited(partEneg);
    mctruth->SetQuenchingFactor(1.0);

    rootIO->SetMCTruth(*mctruth);
    rootIO->WriteTrigger();
}

// Method to AddHitInformation to the root output (used in steppingAction).
void G4CCMUserEventAction::AddHit(G4int row, G4int col, G4bool coat, G4double eneg, G4double time, G4double angle, G4String creatorProcess) {
    // test output string for matching to AddHit output.
    mctruth->AddHitInformation(row,col,coat,eneg,time,angle,true);
}
