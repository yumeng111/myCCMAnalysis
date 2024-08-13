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
/// \file optical/LXe/include/LXeRunAction.hh
/// \brief Definition of the LXeRunAction class
//
//
#include "G4UserRunAction.hh"
#include "G4CCMPrimaryGeneratorAction.h"

#ifndef G4CCMRunAction_h
#define G4CCMRunAction_h 1

class G4CCMRun;
//class G4CCMHistoManager;

class G4Run;

class G4CCMRunAction : public G4UserRunAction
{
 public:
  G4CCMRunAction(G4CCMPrimaryGeneratorAction*);
  ~G4CCMRunAction() override;

  G4Run* GenerateRun() override;
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  G4CCMRun* fRun = nullptr;
  G4CCMPrimaryGeneratorAction* fPrimary = nullptr;
  //G4CCMHistoManager* fHistoManager = nullptr;
};

#endif

