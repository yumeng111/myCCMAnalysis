/*
   header file for stepping action

   does not contain any actions. All methods done in the stepping action made from the .cc code in the UserSteppingAction method
*/

#ifndef g4_larsim_G4CCMUserSteppingAction_H
#define g4_larsim_G4CCMUserSteppingAction_H

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class eventAction;
class G4LogicalVolume;

class G4CCMUserSteppingAction : public G4UserSteppingAction {
public:
    G4CCMUserSteppingAction(eventAction* eventAction);
    virtual ~G4CCMUserSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

private:
    eventAction* fEventAction;
};

#endif // g4_larsim_G4CCMUserSteppingAction_H
