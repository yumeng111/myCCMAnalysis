/*
  header file for event actions.

Does not contain any actions.
 */

#ifndef eventAction_h
#define eventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class runAction;

//Event action class

class eventAction : public G4UserEventAction
{
public:
  eventAction(runAction* runact);
  virtual ~eventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

private:
  runAction* fRunAction;
};

#endif
