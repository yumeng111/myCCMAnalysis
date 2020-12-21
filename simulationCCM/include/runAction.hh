/*
  header file for runAction

does nto contain any actions.
 */

#ifndef runAction_h
#define runAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class runAction : public G4UserRunAction
{
public:
  runAction();
  virtual ~runAction();
  
  virtual void BeginOfRunAction(const G4Run* run);
  virtual void   EndOfRunAction(const G4Run* run);

  void SetRootSet(G4bool newset) { rootSet = newset; }
  bool GetRootSet() { return rootSet; }

private:
  G4bool rootSet;
};

#endif
