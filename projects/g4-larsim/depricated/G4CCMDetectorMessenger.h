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

class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;
class G4UIcommand;
class G4UIdirectory;

class G4CCMDetectorMessenger: public G4UImessenger {
  public:

    G4CCMDetectorMessenger(G4CCMDetectorConstruction*);
    ~G4CCMDetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    G4CCMDetectorConstruction* fG4CCMDetector = nullptr;
    G4UIdirectory* fDetectorDir = nullptr;
    G4UIdirectory* fVolumesDir = nullptr;
    G4UIcmdWithADoubleAndUnit* fPmtRadiusCmd = nullptr;
    G4UIcmdWithABool* fLArCmd = nullptr;
    G4UIcommand* fDefaultsCmd = nullptr;


};

#endif
