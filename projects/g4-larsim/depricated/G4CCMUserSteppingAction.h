#ifndef G4CCMUSERSTEPPINGACTION_H_INCLUDED
#define G4CCMUSERSTEPPINGACTION_H_INCLUDED

#include "G4UserSteppingAction.hh"

/**
 * Implementation of G4UserSteppingAction. This class kills gammas below threshold (set by G4IceTopTank).
 */
class G4CCMUserSteppingAction : public G4UserSteppingAction {

  public:
    G4CCMUserSteppingAction();
   ~G4CCMUserSteppingAction() {}

    void UserSteppingAction(const G4Step*);
};

#endif  // G4CCMUSERSTEPPINGACTION_H_INCLUDED
