/*
  header file for event actions.

Does not contain any actions.
 */

#ifndef eventAction_h
#define eventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "CCMRootIO.h"

class TObjArray;
class TFile;
class TTree;
class TBranch;

class MCTruth;

class runAction;

//Event action class

class eventAction : public G4UserEventAction
{
public:
  eventAction(runAction* runact);
  virtual ~eventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  void AddHit( G4int , G4int , G4bool , G4double , G4double , G4double );

private:
  runAction* fRunAction;

  MCTruth    * mctruth;
  //  CCMRootIO  * rootIO;

};

#endif
