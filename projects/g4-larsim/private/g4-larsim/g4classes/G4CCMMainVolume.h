#ifndef G4CCMMainVolume_h
#define G4CCMMainVolume_h 1

#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"
#include "g4-larsim/g4classes/J4PMTSolidMaker.h"

#include <map>

#include "dataclasses/I3Position.h"
#include "dataclasses/I3Orientation.h"
#include "dataclasses/I3Map.h"
#include "icetray/CCMPMTKey.h"
#include "dataclasses/geometry/CCMGeometry.h"

#include "G4PVPlacement.hh"

class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;

class G4CCMMainVolume : public G4PVPlacement
{
  public:
    G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                G4CCMDetectorConstruction* c);

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
    G4VSolid* fPMT = nullptr; 

    // Logical volumes
    G4LogicalVolume* fCryoVessel_log = nullptr;
    G4LogicalVolume* fVacuum_log = nullptr;
    G4LogicalVolume* fInnerJacket_log = nullptr;
    G4LogicalVolume* fArgonOuter_log = nullptr;
    G4LogicalVolume* fInnerFrame_log = nullptr;
    G4LogicalVolume* fTPBFoil_log = nullptr;
    G4LogicalVolume* fFiducialAr_log = nullptr;
    G4LogicalVolume* fPMT_log = nullptr;

    // Sensitive Detectors positions
    std::vector<G4ThreeVector> fPMTPositions;
};

#endif

