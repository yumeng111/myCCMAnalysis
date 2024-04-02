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
/// \file optical/LXe/src/LXeRunAction.cc
/// \brief Implementation of the LXeRunAction class
//
//
#include "g4-larsim/g4classes/G4CCMRunAction.h"

//#include "g4-larsim/g4classes/G4CCMHistoManager.h"
#include "g4-larsim/g4classes/G4CCMRun.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMRunAction::G4CCMRunAction()
{
  // Book predefined histograms
  //fHistoManager = new G4CCMHistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4CCMRunAction::~G4CCMRunAction() { delete fHistoManager; }
G4CCMRunAction::~G4CCMRunAction() { }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Run* G4CCMRunAction::GenerateRun()
{
  fRun = new G4CCMRun();
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMRunAction::BeginOfRunAction(const G4Run*)
{
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //if(analysisManager->IsActive())
  //{
  //  analysisManager->OpenFile();
  //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMRunAction::EndOfRunAction(const G4Run*)
{
  if(isMaster)
    fRun->EndOfRun();

  // save histograms
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //if(analysisManager->IsActive())
  //{
  //  analysisManager->Write();
  //  analysisManager->CloseFile();
  //}
}

