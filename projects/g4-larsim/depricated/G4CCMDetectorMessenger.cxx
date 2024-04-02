#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"


#include "G4RunManager.hh"
#include "G4Scintillation.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMDetectorMessenger::G4CCMDetectorMessenger(G4CCMDetectorConstruction* detector)
  : fG4CCMDetector(detector)
{
  // Setup a command directory for detector controls with guidance
  fDetectorDir = new G4UIdirectory("/CCM/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fVolumesDir = new G4UIdirectory("/CCM/detector/volumes/");
  fVolumesDir->SetGuidance("Enable/disable volumes");

  fLArCmd = new G4UIcmdWithABool("/CCM/detector/volumes/lar", this);
  fLArCmd->SetGuidance("Enable/Disable the main detector volume.");
  fLArCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fLArCmd->SetToBeBroadcasted(false);

  fDefaultsCmd = new G4UIcommand("/CCM/detector/defaults", this);
  fDefaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  fDefaultsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDefaultsCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMDetectorMessenger::~G4CCMDetectorMessenger()
{
  delete fLArCmd;
  delete fDefaultsCmd;
  delete fDetectorDir;
  delete fVolumesDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fDefaultsCmd) {
    fG4CCMDetector->SetDefaults();
    G4RunManager::GetRunManager()->ReinitializeGeometry();
  }
}


