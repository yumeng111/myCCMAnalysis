#include <G4CCMActionInitialization.h>
#include <G4CCMDetectorConstruction.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMActionInitialization::G4CCMActionInitialization(
    G4CCMDetectorConstruction const * det)
    : G4VUserActionInitialization()
    , fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMActionInitialization::BuildForMaster() const
{
    SetUserAction(new G4CCMPrimaryGeneratorAction());
    //SetUserAction(new G4CCMRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMActionInitialization::Build() const
{
    SetUserAction(new G4CCMPrimaryGeneratorAction());
    //SetUserAction(new G4CCMRunAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



