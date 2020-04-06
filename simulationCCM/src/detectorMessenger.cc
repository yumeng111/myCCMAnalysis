/*
  detector Messenger for CCM simulation

this code creates the commands that can be used to modify the detector from macro files.
it should be noted that the more complicated these commands are made, the higher the chance of memory leaks.
Be careful trying to make a detector that can be completely controlled from the ui.

*/
#include "detectorMessenger.hh"
#include "detectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Scintillation.hh"

#include "G4RunManager.hh"

//Constructor
detectorMessenger::detectorMessenger(detectorConstruction* detector)
 : fdetector(detector)
{
  //Setup a command directory for detector controls with guidance
  fDetectorDir = new G4UIdirectory("/ccm/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fVolumesDir = new G4UIdirectory("/ccm/detector/volumes/");
  fVolumesDir->SetGuidance("Enable/disable volumes");
 
  //Various commands for modifying detector geometry
  fPmtRadiusCmd = new G4UIcmdWithADoubleAndUnit
    ("/ccm/detector/pmtRadius",this);
  fPmtRadiusCmd->SetGuidance("Set the radius of the PMTs.");
  fPmtRadiusCmd->SetParameterName("radius",false);
  fPmtRadiusCmd->SetDefaultUnit("cm");
  fPmtRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fPmtRadiusCmd->SetToBeBroadcasted(false);

  fLaserCmd = new G4UIcmdWithADoubleAndUnit
    ("/ccm/detector/rodHeight",this);
  fLaserCmd->SetGuidance("Set the height of the laser rod from the middle.");
  fLaserCmd->SetParameterName("height",false);
  fLaserCmd->SetDefaultUnit("cm");
  fLaserCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fLaserCmd->SetToBeBroadcasted(false);

  fPMTsCmd = new G4UIcmdWithABool("/ccm/detector/pmts",this);
  fPMTsCmd->SetGuidance("Enable/Disable the PMTs.");
  fPMTsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fPMTsCmd->SetToBeBroadcasted(false);

  fSodiumCmd = new G4UIcmdWithABool("/ccm/detector/sodium",this);
  fSodiumCmd->SetGuidance("Enable/Disable the Sodium Source.");
  fSodiumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSodiumCmd->SetToBeBroadcasted(false);

  f200Cmd = new G4UIcmdWithABool("/ccm/detector/set200",this);
  f200Cmd->SetGuidance("Enable/Disable the 200 pmt detector.");
  f200Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  f200Cmd->SetToBeBroadcasted(false);

  fRodCmd = new G4UIcmdWithABool("/ccm/detector/setrod",this);
  fRodCmd->SetGuidance("Enable/Disable the Calibration Rod.");
  fRodCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fRodCmd->SetToBeBroadcasted(false);

  fTPBfoilCmd = new G4UIcmdWithABool("/ccm/detector/tpbfoil",this);
  fTPBfoilCmd->SetGuidance("Enable/Disable the TPBfoil. For LBOC design, Enable/Disable the triple cylinders.");
  fTPBfoilCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fTPBfoilCmd->SetToBeBroadcasted(false);

  fReflectorCmd = new G4UIcmdWithABool("/ccm/detector/reflector",this);
  fReflectorCmd->SetGuidance("Enable/Disable the Reflector.");
  fReflectorCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fReflectorCmd->SetToBeBroadcasted(false);

  fDefaultsCmd = new G4UIcommand("/ccm/detector/defaults",this);
  fDefaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  fDefaultsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDefaultsCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Deconstructor, deletes all commands from the ui so they don't impact other geant4 simulations.
detectorMessenger::~detectorMessenger()
{
  delete fPmtRadiusCmd;
  delete fLaserCmd;
  delete fDefaultsCmd;
  delete fPMTsCmd;
  delete fSodiumCmd;
  delete fRodCmd;
  delete f200Cmd;
  delete fTPBfoilCmd;
  delete fReflectorCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Set the actions to be taken for each command.
void detectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fPmtRadiusCmd){
    fdetector->SetPMTRadius(fPmtRadiusCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fLaserCmd){
    fdetector->SetRodHeight(fLaserCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fPMTsCmd){
    fdetector->SetPMTsOn(fPMTsCmd->GetNewBoolValue(newValue));
  }
  else if (command == fSodiumCmd){
    fdetector->SetSodiumOn(fSodiumCmd->GetNewBoolValue(newValue));
  }
  else if (command == f200Cmd){
    fdetector->SetCCM200(f200Cmd->GetNewBoolValue(newValue));
  }
  else if (command == fRodCmd){
    fdetector->SetRodin(fRodCmd->GetNewBoolValue(newValue));
  }
  else if (command == fTPBfoilCmd){
    fdetector->SetTPBfoilOn(fTPBfoilCmd->GetNewBoolValue(newValue));
  }
  else if (command == fReflectorCmd){
    fdetector->SetReflectorOn(fReflectorCmd->GetNewBoolValue(newValue));
  }
  else if (command == fDefaultsCmd){
    fdetector->SetDefaults();

    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}
