#ifndef G4CCMUSERTRACKINGACTION_H_INCLUDED
#define G4CCMUSERTRACKINGACTION_H_INCLUDED

#include "G4UserTrackingAction.hh"

/**
 * Implementation of G4UserTrackingAction. This class kills gammas below threshold (set by G4CCMTank).
 */
class G4CCMUserTrackingAction : public G4UserTrackingAction {

  public:
    G4CCMUserTrackingAction();
   ~G4CCMUserTrackingAction() {}

    void PreUserTrackingAction(const G4Track*);
    void PostUserTrackingAction(const G4Track*);
};

#endif  // G4CCMUSERTRACKINGACTION_H_INCLUDED
