#include "G4CCMEventAction.h"

#include "G4CCMDetectorConstruction.h"
#include "G4CCMHistoManager.h"
#include "G4CCMPMTHit.h"
#include "G4CCMRun.h"
#include "G4CCMScintHit.h"
#include "G4CCMTrajectory.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::LXeEventAction(const LXeDetectorConstruction* det)
      : fDetector(det)
{
      fEventMessenger = new LXeEventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeEventAction::~LXeEventAction() { delete fEventMessenger; }

void LXeEventAction::BeginOfEventAction(const G4Event*)
{
    fHitCount                = 0;
    fPhotonCount_Scint       = 0;
    fPhotonCount_Ceren       = 0;
    fAbsorptionCount         = 0;
    fBoundaryAbsorptionCount = 0;
    fTotE                    = 0.0;

    fConvPosSet = false;
    fEdepMax    = 0.0;

    fPMTsAboveThreshold = 0;

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    if(fScintCollID < 0)
        fScintCollID = SDman->GetCollectionID("scintCollection");
    if(fPMTCollID < 0)
        fPMTCollID = SDman->GetCollectionID("pmtHitCollection");
}




