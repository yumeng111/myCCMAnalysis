/*
Event Action for CCM simulation.

This code offers the ability to perform actions at the start/end of events. 
Currently it does nothing, as most outputs are handled through the detector and stepping action.
*/
#include "eventAction.hh"
#include "runAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"

//Constructor
eventAction::eventAction(runAction* runaction)
  : G4UserEventAction(),
    fRunAction(runaction)
{}

//Deconstructor
eventAction::~eventAction()
{}

//Within this function define things to be done at the start of an event (example: the debugging line currently commented out).
//The code currently does not access the detector or the primary generator, so information produced there is not accessible.
void eventAction::BeginOfEventAction(const G4Event*)
{
  //G4cout << "maincodenote: starting eventAction" << G4endl;
}

//Same as BeginOfEventAction, but for end of event.
void eventAction::EndOfEventAction(const G4Event*)
{
  //G4cout << "maincodenote: ending eventAction" << G4endl;
  //G4cout << "Event Ended" << G4endl;
}
