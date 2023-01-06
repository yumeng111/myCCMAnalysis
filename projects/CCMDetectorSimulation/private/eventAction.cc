/*
Event Action for CCM simulation.

This code offers the ability to perform actions at the start/end of events. 
Currently handles the output to ROOT, using a file defined in the BeginOfEventAction, 
trigger information obtained in the EndOfEventAction, and defines a method for adding
a hit (optical photon to PMT).
*/
#include <sstream>
#include <iostream>

#include "CCMAnalysis/CCMDetectorSimulation/eventAction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/runAction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/primaryGenerator.hh"
#include "CCMAnalysis/CCMDetectorSimulation/detectorConstruction.hh"

#include "CCMAnalysis/CCMDataStructures/MCTruth.h"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

//Constructor
eventAction::eventAction(runAction* runaction)
  : G4UserEventAction(), fRunAction(runaction),
    mctruth(new MCTruth()), rootIO(std::make_shared<CCMRootIO>())
{
  rootSet = false;
}

//Deconstructor
eventAction::~eventAction()
{
// clear data and root file structures
  delete mctruth;
  rootIO->Close();
}

//Within this function define things to be done at the start of an event (example: the debugging line currently commented out).
//The code currently does not access the detector or the primary generator, so information produced there is not accessible.
void eventAction::BeginOfEventAction(const G4Event* event)
{
  //get the event ID, reset the rootIO for that event, and determine if the root filename is already set.
  const G4int eventID = event->GetEventID();

  //  G4cout << "Beginning event: " << eventID << G4endl;

  mctruth->Reset(eventID);
  rootSet = fRunAction->GetRootSet();

  //if rootSet is false, pull the root file name and use it to open a new root file and rootIO structure
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

//Same as BeginOfEventAction, but for end of event.
void eventAction::EndOfEventAction(const G4Event* event)
{
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
  //G4cout << "EventActionInitializationString:\t" << xpos << '\t' << ypos << '\t' << zpos << '\t' << momx << '\t' << momy << '\t' << momz << '\t' << partEneg << '\t' << partName << '\t' << G4endl;
  
  //output the initial conditions to the rootIO and write for the Trigger.
  mctruth->SetParticlePosition(xpos,ypos,zpos);
  mctruth->SetParticleMomentum(momx,momy,momz);
  mctruth->SetParticleID(eventID);
  mctruth->SetTotalEnergyDeposited(partEneg);
  mctruth->SetQuenchingFactor(1.0);

  rootIO->SetMCTruth(*mctruth);
  rootIO->WriteTrigger();
}

//Method to AddHitInformation to the root output (used in steppingAction).
void eventAction::AddHit( G4int row, G4int col, G4bool coat, G4double eneg, G4double time, G4double angle, G4String creatorProcess)
{
  //test output string for matching to AddHit output.
  //G4cout << row << '\t' << col << '\t' << coat << '\t' << eneg << '\t' << time << '\t' << angle << G4endl;
  mctruth->AddHitInformation(row,col,coat,eneg,time,angle,true);
}
