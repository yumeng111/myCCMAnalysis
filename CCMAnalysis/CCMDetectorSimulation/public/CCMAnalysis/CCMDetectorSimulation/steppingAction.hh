/*
  header file for stepping action

does not contain any actions. All methods done in the stepping action made from the .cc code in the UserSteppingAction method
 */

#ifndef steppingAction_h
#define steppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class eventAction;
class G4LogicalVolume;

class steppingAction : public G4UserSteppingAction
{
public:
  steppingAction(eventAction* eventAction);
  virtual ~steppingAction();
  
  virtual void UserSteppingAction(const G4Step*);
  
private:
  eventAction* fEventAction;
};

#endif
