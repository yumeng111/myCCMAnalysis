/*
Run action for CCM simulation

This code defines things to be done at the start or end of a run.
Currently it prints an event number every 1000 runs.
Also does a few variable translations and print statements
 */

#include "runAction.hh"
#include "primaryGenerator.hh"
#include "detectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

//constructor. 
runAction::runAction()
  : G4UserRunAction()
{
  //set printing event number per each 1000 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//Deconstructor
runAction::~runAction()
{}

//Define things to do at the beginning of a run.
//Currently determines if the OM values have been randomized and prints them to a text file.
//Also pulls a detector variable that tells if the output root file has been set and to what
void runAction::BeginOfRunAction(const G4Run*)
{
  //G4cout << "maincodenote: starting runAction" << G4endl;
  //calls the detector and obtains the root file name.
  const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  rootSet = detector->GetRootSet();

  //determines if the detector OM is randomized and prints the variable list if it is
  if (detector->IsRandom()) {
    G4String randoms = detector->GetRandoms();
    G4cout << randoms << G4endl;
  }

  //Do not save the random seeds for every event.
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//Define things to do at the end of a run. currently nothing.
void runAction::EndOfRunAction(const G4Run*)
{  //G4cout << "maincodenote: ending runAction" << G4endl;
}
