/*
Run action for CCM simulation

This code defines things to be done at the start or end of a run.
Currently it prints an event number every 1000 runs.
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

//Define things to do at the beginning of a run. currently does nothing.
void runAction::BeginOfRunAction(const G4Run*)
{
  //G4cout << "maincodenote: starting runAction" << G4endl;
  const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //detectorConstruction* detector = new detectorConstruction;

  rootSet = detector->GetRootSet();

  if (detector->IsRandom()) {
    G4String randoms = detector->GetRandoms();
    G4cout << randoms << G4endl;
  }//*/

  //Do not save the random seeds for every event.
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

//Define things to do at the end of a run. 
void runAction::EndOfRunAction(const G4Run*)
{  //G4cout << "maincodenote: ending runAction" << G4endl;
}
