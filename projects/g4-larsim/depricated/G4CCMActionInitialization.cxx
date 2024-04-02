#include "G4CCMActionInitialization.h"
#include "G4CCMDetectorConstruction.h"
#include "G4CCMEventAction.h"
#include "G4CCMPrimaryGeneratorAction.h"
#include "G4CCMRunAction.h"
#include "G4CCMStackingAction.h"
#include "G4CCMSteppingAction.h"
#include "G4CCMTrackingAction.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMActionInitialization::G4CCMActionInitialization(
    G4CCMDetectorConstruction const * det)
    : G4VUserActionInitialization()
    , fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMActionInitialization::BuildForMaster() const
{
    SetUserAction(new G4CCMRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMActionInitialization::Build() const
{
    SetUserAction(new G4CCMPrimaryGeneratorAction());

    auto eventAction = new G4CCMEventAction(fDetector);
    SetUserAction(eventAction);
    SetUserAction(new G4CCMStackingAction(eventAction));

    SetUserAction(new G4CCMRunAction());
    SetUserAction(new G4CCMTrackingAction());
    SetUserAction(new G4CCMSteppingAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



