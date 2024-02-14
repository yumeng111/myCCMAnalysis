/*
   header file for event actions.

   Does not contain any actions.
*/

#ifndef g4_larsim_G4CCMUserEventAction_H
#define g4_larsim_G4CCMUserEventAction_H

#include <memory>

#include "globals.hh"

#include "G4UserEventAction.hh"

class TObjArray;
class TFile;
class TTree;
class TBranch;

class MCTruth;

class G4CCMUserRunAction;

// Event action class
class G4CCMUserEventAction : public G4UserEventAction {
public:
    G4CCMUserEventAction(G4CCMUserRunAction* runact);
    virtual ~G4CCMUserEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void AddHit( G4int, G4int, G4bool, G4double, G4double, G4double, G4String);

private:
    G4CCMUserRunAction * fRunAction;

    MCTruth * mctruth;
};

#endif
