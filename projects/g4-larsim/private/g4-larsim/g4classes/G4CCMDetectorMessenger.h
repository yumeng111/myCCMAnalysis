/*
  header file for detector Messenger

initiates variables for all commands defined for the detector through the UI
Also initiates the methods needed for using the commands and the format of each command.
 */
#ifndef G4CCMDetectorMessenger_h
#define G4CCMDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4CCMDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcommand;
class G4UIcmdWithABool;
class G4UIcmdWithAString;

class G4CCMDetectorMessenger: public G4UImessenger {
  public:

    G4CCMDetectorMessenger(G4CCMDetectorConstruction*);
    virtual ~G4CCMDetectorMessenger();
 
    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:

    G4CCMDetectorConstruction*        fdetector;
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
    G4UIcmdWithABool*            fCosmicCmd;
    G4UIcmdWithABool*            fInvBetaCmd;
    G4UIcmdWithABool*            fDMCmd;
    G4UIcmdWithABool*            fALPCmd;
    G4UIcmdWithABool*            f200Cmd;
    G4UIcmdWithABool*            fTPBfoilCmd;
    G4UIcmdWithABool*            fReflectorCmd;
    G4UIcmdWithAString*          fRootCmd;
    G4UIcmdWithAnInteger*        fOrCmd;

};

#endif
