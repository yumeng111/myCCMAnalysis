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
/// \file optical/LXe/LXe.cc
/// \brief Main program of the optical/LXe example
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "../visualization/G4CCMVisActionInitialization.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"
#include "g4-larsim/g4classes/G4CCMPhysicsList.h"

#include "FTFP_BERT.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include <G4SystemOfUnits.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if(argc == 1)
  {
    ui = new G4UIExecutive(argc, argv);
  }

  auto runManager = G4RunManagerFactory::CreateRunManager();

  auto det = new G4CCMDetectorConstruction(8.2 * ns, 743.0 * ns, 94.0 * cm, 0.605, 0.605, 0.605, 0.002 * mm, 0.002 * mm, 0.002 * mm, 0.13457, 8.13914e-21, 1.0, 0.99, 0.8, true);
  runManager->SetUserInitialization(det);

  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  //physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  //auto opticalPhysics = new G4OpticalPhysics();
  //auto opticalParams  = G4OpticalParameters::Instance();

  //opticalParams->SetWLSTimeProfile("delta");

  //opticalParams->SetScintTrackSecondariesFirst(true);

  //opticalParams->SetCerenkovMaxPhotonsPerStep(100);
  //opticalParams->SetCerenkovMaxBetaChange(10.0);
  //opticalParams->SetCerenkovTrackSecondariesFirst(true);

  //physicsList->RegisterPhysics(opticalPhysics);
  //runManager->SetUserInitialization(physicsList);
  runManager->SetUserInitialization(new G4CCMPhysicsList(1));

  runManager->SetUserInitialization(new G4CCMVisActionInitialization(det));

  // initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(ui)
  {
    // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    if(ui->IsGUI())
    {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
  else
  {
    // batch mode
    G4String command  = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }

  // job termination
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
