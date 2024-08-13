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
/// \file optical/LXe/src/LXeEventMessenger.cc
/// \brief Implementation of the LXeEventMessenger class
//
//
#include "G4CCMEventMessenger.h"
#include "G4CCMEventAction.h"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMEventMessenger::G4CCMEventMessenger(G4CCMEventAction* event)
  : fG4CCMEvent(event)
{
  fVerboseCmd = new G4UIcmdWithAnInteger("/CCM/eventVerbose", this);
  fVerboseCmd->SetGuidance("Set the verbosity of event data.");
  fVerboseCmd->SetParameterName("verbose", true);
  fVerboseCmd->SetDefaultValue(1);

  fPmtThresholdCmd = new G4UIcmdWithAnInteger("/CCM/pmtThreshold", this);
  fPmtThresholdCmd->SetGuidance("Set the pmtThreshold (in # of photons)");

  fForceDrawPhotonsCmd = new G4UIcmdWithABool("/CCM/forceDrawPhotons", this);
  fForceDrawPhotonsCmd->SetGuidance("Force drawing of photons.");
  fForceDrawPhotonsCmd->SetGuidance(
    "(Higher priority than /CCM/forceDrawNoPhotons)");

  fForceDrawNoPhotonsCmd =
    new G4UIcmdWithABool("/CCM/forceDrawNoPhotons", this);
  fForceDrawNoPhotonsCmd->SetGuidance("Force no drawing of photons.");
  fForceDrawNoPhotonsCmd->SetGuidance(
    "(Lower priority than /CCM/forceDrawPhotons)");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMEventMessenger::~G4CCMEventMessenger()
{
  delete fVerboseCmd;
  delete fPmtThresholdCmd;
  delete fForceDrawPhotonsCmd;
  delete fForceDrawNoPhotonsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMEventMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fVerboseCmd)
  {
    fG4CCMEvent->SetEventVerbose(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if(command == fPmtThresholdCmd)
  {
    fG4CCMEvent->SetPMTThreshold(fPmtThresholdCmd->GetNewIntValue(newValue));
  }
  else if(command == fForceDrawPhotonsCmd)
  {
    fG4CCMEvent->SetForceDrawPhotons(
      fForceDrawPhotonsCmd->GetNewBoolValue(newValue));
  }
  else if(command == fForceDrawNoPhotonsCmd)
  {
    fG4CCMEvent->SetForceDrawNoPhotons(
      fForceDrawNoPhotonsCmd->GetNewBoolValue(newValue));
  }
}

