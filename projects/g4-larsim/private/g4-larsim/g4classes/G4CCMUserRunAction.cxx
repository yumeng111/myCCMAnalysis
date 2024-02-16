/*
Run action for CCM simulation

This code defines things to be done at the start or end of a run.
Currently it prints an event number every 1000 runs.
Also does a few variable translations and print statements
 */

#include "g4-larsim/g4classes/G4CCMUserRunAction.hh"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

G4CCMUserRunAction::G4CCMUserRunAction()
  : G4UserRunAction() {
  // set printing event number per each 1000 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

G4CCMUserRunAction::~G4CCMUserRunAction() {}

// Define things to do at the beginning of a run.
// Currently determines if the OM values have been randomized and prints them to a text file.
// Also pulls a detector variable that tells if the output root file has been set and to what
void G4CCMUserRunAction::BeginOfRunAction(const G4Run*) {
  // calls the detector and obtains the root file name.
  const G4CCMDetectorConstruction* detector = static_cast<const G4CCMDetectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // determines if the detector OM is randomized and prints the variable list if it is
  if (detector->IsRandom()) {
    G4String randoms = detector->GetRandoms();
    G4cout << randoms << G4endl;
  }

  // Do not save the random seeds for every event.
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

// Define things to do at the end of a run. currently nothing.
void G4CCMUserRunAction::EndOfRunAction(const G4Run*) {}
