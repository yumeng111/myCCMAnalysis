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

  void shootArgon(G4Event*);
  void shootDarkMatter(G4Event*);
  void shootSodium(G4Event*);
  void shootCosmic(G4Event*);

  const G4ParticleGun* GetParticleGun() const {return fParticleGun; }
  //defines a method to access the fParticleGun variable from any part of the code is the primary generator is included.

  G4double GetXpos() const { return xpos; }
  G4double GetYpos() const { return ypos; }
  G4double GetZpos() const { return zpos; }
  G4double GetXmom() const { return momx; }
  G4double GetYmom() const { return momy; }
  G4double GetZmom() const { return momz; }
  G4double GetPartEneg() const { return partEneg; }
  G4String GetPartName() const { return nameString; }

private:
  G4ParticleGun* fParticleGun;

  G4double xpos;
  G4double ypos;
  G4double zpos;
  G4double momx;
  G4double momy;
  G4double momz;
  G4double partEneg;
  G4String nameString;
};

#endif
