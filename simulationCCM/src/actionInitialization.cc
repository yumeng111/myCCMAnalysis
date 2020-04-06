/*
Action initialization code for the CCM simulation.

This code is mostly a container that connects primary generator and run, event, and stepping action 
(called at the start of each run, start of each event, and start of each step respectively)
rather than a code which does anything itself. 

*/

#include "actionInitialization.hh"
#include "primaryGenerator.hh"
#include "runAction.hh"
#include "eventAction.hh"
#include "steppingAction.hh"
//#include "stackingAction.hh"
#include "detectorConstruction.hh"

//Constructor
actionInitialization::actionInitialization()
  : G4VUserActionInitialization()
{}

//Deconstructor
actionInitialization::~actionInitialization()
{}

//initialization of runAction
void actionInitialization::BuildForMaster() const
{
  runAction* runact = new runAction;
  SetUserAction(runact);
}

//Build the action for a run. All relatively self-explanatory, with two debugging lines.
void actionInitialization::Build() const
{
  //G4cout << "maincodenote: starting actionInitialization" << G4endl;

  SetUserAction(new primaryGenerator);

  runAction* runact = new runAction;
  SetUserAction(runact);
  
  eventAction* eventaction = new eventAction(runact);
  SetUserAction(eventaction);
  
  SetUserAction(new steppingAction(eventaction));

  //G4cout << "maincodenote: ending actionInitialization" << G4endl;
  
}
