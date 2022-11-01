/*
  Header file for action Initialization.

Does not contain any actions.
 */

#ifndef actionInitialization_h
#define actionInitialization_h 1

#include "G4VUserActionInitialization.hh"

//action initialization class.

class actionInitialization : public G4VUserActionInitialization{

public:
  actionInitialization();
  virtual ~actionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;
};

#endif
