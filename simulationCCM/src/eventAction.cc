/*
Event Action for CCM simulation.

This code offers the ability to perform actions at the start/end of events. 
Currently it does nothing, as most outputs are handled through the detector and stepping action.
*/
#include "eventAction.hh"
#include "runAction.hh"
#include "primaryGenerator.hh"
#include "detectorConstruction.hh"

#include "MCTruth.h"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

//Constructor
eventAction::eventAction(runAction* runaction)
  : G4UserEventAction(), fRunAction(runaction),
    mctruth(new MCTruth())
{}

//Deconstructor
eventAction::~eventAction()
{
  // clear data
  delete mctruth;
}

//Within this function define things to be done at the start of an event (example: the debugging line currently commented out).
//The code currently does not access the detector or the primary generator, so information produced there is not accessible.
void eventAction::BeginOfEventAction(const G4Event* event)
{
  const G4int eventID = event->GetEventID();

  mctruth->Reset(eventID);
  //G4cout << "maincodenote: starting eventAction" << G4endl;
}

//Same as BeginOfEventAction, but for end of event.
void eventAction::EndOfEventAction(const G4Event* event)
{
  auto rootIO = CCMRootIO::GetInstance();

  const G4int eventID = event->GetEventID();

  //pull the initial conditions from the primary Generator for usage in ROOT output.
  const primaryGenerator* generator = static_cast<const primaryGenerator*> (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

  G4double xpos = generator->GetXpos();
  G4double ypos = generator->GetYpos();
  G4double zpos = generator->GetZpos();
  G4double momx = generator->GetXmom();
  G4double momy = generator->GetYmom();
  G4double momz = generator->GetZmom();
  G4double partEneg = generator->GetPartEneg();
  G4String partName = generator->GetPartName();
  
  //test output string for matching to primary generator notes
  G4cout << "EventActionInitializationString:\t" << xpos << '\t' << ypos << '\t' << zpos << '\t' << momx << '\t' << momy << '\t' << momz << '\t' << partEneg << '\t' << partName << '\t' << G4endl;
  
  const detectorConstruction* detector = static_cast<const detectorConstruction*> (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4String filename = detector->GetRootFile();

  mctruth->SetParticlePosition(xpos,ypos,zpos);
  mctruth->SetParticleMomentum(momx,momy,momz);
  mctruth->SetParticleID(eventID);
  mctruth->SetTotalEnergyDeposited(partEneg);

  rootIO->SetParameter("SaveBranches","mctruth");
  rootIO->SetOutFileName(filename);
  rootIO->SetupOutputFile();
  
  rootIO->SetMCTruth(*mctruth);
  rootIO->WriteTrigger();

  rootIO->Close();
  //G4cout << "maincodenote: ending eventAction" << G4endl;
  //G4cout << "Event Ended" << G4endl;*/

  //delete rootIO;
}


void eventAction::AddHit( G4int row, G4int col, G4bool coat, G4double eneg, G4double time, G4double angle)
{
  G4cout << row << '\t' << col << '\t' << coat << '\t' << eneg << '\t' << time << '\t' << angle << G4endl;
  mctruth->AddHitInformation(row,col,coat,eneg,time,angle);
}
