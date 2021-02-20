/*
  header file for detector Messenger

initiates variables for all commands defined for the detector through the UI
Also initiates the methods needed for using the commands and the format of each command.
 */
#ifndef detectorMessenger_h
#define detectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class detectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class detectorMessenger: public G4UImessenger
{
  public:

    detectorMessenger(detectorConstruction*);
    virtual ~detectorMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:

    detectorConstruction*        fdetector;
    G4UIdirectory*               fDetectorDir;
    G4UIdirectory*               fVolumesDir;
    G4UIcmdWithADoubleAndUnit*   fPmtRadiusCmd;
    G4UIcmdWithADoubleAndUnit*   fLaserCmd;
    G4UIcmdWithADoubleAndUnit*   fRodCmd;
    G4UIcommand*                 fDefaultsCmd;
    G4UIcommand*                 fRandomsCmd;
    G4UIcommand*                 fCorrelateCmd;
    G4UIcommand*                 fCleanCmd;
    G4UIcmdWithABool*            fPMTsCmd;
    G4UIcmdWithABool*            fSodiumCmd;
    G4UIcmdWithABool*            fAr39Cmd;
    G4UIcmdWithABool*            fDMCmd;
    G4UIcmdWithABool*            f200Cmd;
    G4UIcmdWithABool*            fTPBfoilCmd;
    G4UIcmdWithABool*            fReflectorCmd;
    G4UIcmdWithAString*          fRootCmd;
    G4UIcmdWithAnInteger*        fOrCmd;

};

#endif
