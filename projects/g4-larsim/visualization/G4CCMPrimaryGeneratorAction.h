//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file optical/LXe/include/LXePrimaryGeneratorAction.hh
/// \brief Definition of the LXePrimaryGeneratorAction class
//
//
#ifndef G4CCMPrimaryGeneratorAction_h
#define G4CCMPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

class G4ParticleGun;
class G4Event;

class G4CCMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  G4CCMPrimaryGeneratorAction();
  ~G4CCMPrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* anEvent) override;

  // method to access particle gps
  const G4GeneralParticleSource *GetGeneralParticleSource() const { return fGeneralParticleSource; }

 private:
  G4GeneralParticleSource *fGeneralParticleSource; // pointer a to G4 GeneralParticleSource class
  G4ParticleGun* fParticleGun = nullptr;
};

#endif

