/*
  header file for primary Generator

initiates the various methods used by the primary generator. 
Also defines a method to call the variable fParticleGun
 */

#ifndef primaryGenerator_h
#define primaryGenerator_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;

class primaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  primaryGenerator();
  virtual ~primaryGenerator();
  //constructor and deconstructor.
  
  virtual void GeneratePrimaries(G4Event*);
  //method for actually generating primary particles

  const G4ParticleGun* GetParticleGun() const {return fParticleGun; }
  //defines a method to access the fParticleGun variable from any part of the code is the primary generator is included.

private:
  G4ParticleGun* fParticleGun;
};

#endif
