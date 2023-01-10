/*
  header file for event actions.

Does not contain any actions.
 */

#ifndef eventAction_h
#define eventAction_h 1

#include <memory>

#include "globals.hh"

#include "CCMAnalysis/CCMIO/CCMRootIO.h"
#include "G4UserEventAction.hh"

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

  void AddHit( G4int , G4int , G4bool , G4double , G4double , G4double , G4String );

private:
  runAction* fRunAction;

  MCTruth    * mctruth;
  std::shared_ptr<CCMRootIO> rootIO;

  G4bool rootSet;

};

#endif
