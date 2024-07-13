#ifndef G4CCMMainVolume_h
#define G4CCMMainVolume_h 1

#include "dataclasses/I3Map.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Orientation.h"
#include "dataclasses/geometry/CCMGeometry.h"

#include "g4-larsim/g4classes/J4PMTSolidMaker.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

#include "icetray/CCMPMTKey.h"

#include <map>

#include <G4PVPlacement.hh>

class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;

class G4CCMMainVolume : public G4PVPlacement
{
  public:
    G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                G4CCMDetectorConstruction* c, G4bool SodiumSourceOn, G4double SodiumSourceLocation);

    G4LogicalVolume* GetLogPMTCoated() { return fPMTCoated_log; }
    G4LogicalVolume* GetLogPMTUncoated() { return fPMTUncoated_log; }
    G4LogicalVolume* GetLogScint() { return fFiducialAr_log; }
    G4LogicalVolume* GetLogSodiumPellet() { return fSodiumSourcePellet_log; }
    G4LogicalVolume* GetLogSourceRod() { return fSourceRod_log; }
    G4LogicalVolume* GetLogTPBCoating() { return fTPBCoating_log; }
    G4LogicalVolume* GetLogTPBFoil() { return fTPBFoil_log; }
    std::vector<G4ThreeVector> GetPMTPositions() { return fPMTPositions; }

  private:
    void VisAttributes();
    void SurfaceProperties();
    G4CCMDetectorConstruction* fConstructor = nullptr;

    // Basic Volumes
    G4Tubs* fCryoVessel = nullptr;
    G4Tubs* fVacuum = nullptr;
    G4Tubs* fInnerJacket = nullptr;
    G4Tubs* fArgonOuter = nullptr;
    G4Tubs* fInnerFrame = nullptr;
    G4Tubs* fTPBFoil = nullptr;
    G4Tubs* fFiducialAr = nullptr;
    G4VSolid* fPMTCoated = nullptr; 
    G4VSolid* fPMTUncoated = nullptr; 
    G4VSolid* fTPBCoating = nullptr; 
    G4VSolid* fPhotocathCoated = nullptr; 
    G4VSolid* fPhotocathUncoated = nullptr; 
    G4Tubs* fSourceRod = nullptr;
    G4Tubs* fSodiumSourcePellet = nullptr;

    // Logical volumes
    G4LogicalVolume* fCryoVessel_log = nullptr;
    G4LogicalVolume* fVacuum_log = nullptr;
    G4LogicalVolume* fInnerJacket_log = nullptr;
    G4LogicalVolume* fArgonOuter_log = nullptr;
    G4LogicalVolume* fInnerFrame_log = nullptr;
    G4LogicalVolume* fTPBFoil_log = nullptr;
    G4LogicalVolume* fFiducialAr_log = nullptr;
    G4LogicalVolume* fPMTCoated_log = nullptr;
    G4LogicalVolume* fPMTUncoated_log = nullptr;
    G4LogicalVolume* fTPBCoating_log = nullptr; 
    G4LogicalVolume* fPhotocathCoated_log = nullptr;
    G4LogicalVolume* fPhotocathUncoated_log = nullptr;
    G4LogicalVolume* fSourceRod_log = nullptr;
    G4LogicalVolume* fSodiumSourcePellet_log = nullptr;

    // Sensitive Detectors positions
    std::vector<G4ThreeVector> fPMTPositions;

};

#endif

