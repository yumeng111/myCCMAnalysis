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
                G4CCMDetectorConstruction* c, G4bool SourceRodIn, G4double SourceRodLocation, G4bool CobaltSourceRun, G4bool SodiumSourceRun);

    G4LogicalVolume* GetLogPMTCoatedWall() { return fPMTCoatedWall_log; }
    G4LogicalVolume* GetLogPMTCoatedCaps() { return fPMTCoatedCaps_log; }
    G4LogicalVolume* GetLogPMTUncoatedWall() { return fPMTUncoatedWall_log; }
    G4LogicalVolume* GetLogPMTUncoatedCaps() { return fPMTUncoatedCaps_log; }

    G4LogicalVolume* GetLogScint() { return fFiducialAr_log; }
    G4LogicalVolume* GetLogSourcePellet() { return fSourcePellet_log; }
    G4LogicalVolume* GetLogSourceRod() { return fSourceRod_log; }

    G4LogicalVolume* GetLogTPBCoatingWall() { return fTPBCoatingWall_log; }
    G4LogicalVolume* GetLogTPBCoatingCaps() { return fTPBCoatingCaps_log; }

    G4LogicalVolume* GetLogTPBFoilSides() { return fTPBFoilSides_log; }
    G4LogicalVolume* GetLogTPBFoilTop() { return fTPBFoilTop_log; }
    G4LogicalVolume* GetLogTPBFoilBottom() { return fTPBFoilBottom_log; }
    G4LogicalVolume* GetLogReflectorFoil() { return fReflectorFoil_log; }
    G4LogicalVolume* GetShinyC406R0() { return fShinyC406R0_log; }
    G4LogicalVolume* GetShinyTop() { return fShinyTop_log; }
    G4LogicalVolume* GetShinyBottom() { return fShinyBottom_log; }

    std::vector<G4ThreeVector> GetPMTPositions() { return fPMTPositions; }

  private:
    void VisAttributes(G4bool SourceRodIn);
    void SurfaceProperties();
    G4CCMDetectorConstruction* fConstructor = nullptr;

    // Basic Volumes
    G4Tubs* fCryoVessel = nullptr;
    G4Tubs* fVacuum = nullptr;
    G4Tubs* fInnerJacket = nullptr;
    G4Tubs* fArgonOuter = nullptr;
    G4Tubs* fInnerFrame = nullptr;
    G4Tubs* fReflectorFoil = nullptr;
    G4Tubs* fTPBFoilSides = nullptr;
    G4Tubs* fTPBFoilTopBottom = nullptr;
    G4Tubs* fFiducialAr = nullptr;

    G4VSolid* fPMTCoatedWall = nullptr;
    G4VSolid* fPMTCoatedCaps = nullptr;
    G4VSolid* fPMTUncoatedWall = nullptr;
    G4VSolid* fPMTUncoatedCaps = nullptr;

    G4VSolid* fBridleWall = nullptr;
    G4VSolid* fBridleCaps = nullptr;
    G4VSolid* fFrillWall = nullptr;
    G4VSolid* fFrillCaps = nullptr;

    G4VSolid* fTPBCoatingWall = nullptr; 
    G4VSolid* fTPBCoatingCaps = nullptr; 

    G4VSolid* fPhotocathCoated = nullptr; 
    G4VSolid* fPhotocathUncoated = nullptr; 
    G4Tubs* fSourceRod = nullptr;
    G4Tubs* fSourcePellet = nullptr;

    G4Tubs* fShinyC406R0 = nullptr;
    G4Tubs* fShinyTop = nullptr;
    G4Tubs* fShinyBottom = nullptr;

    // Logical volumes
    G4LogicalVolume* fCryoVessel_log = nullptr;
    G4LogicalVolume* fVacuum_log = nullptr;
    G4LogicalVolume* fInnerJacket_log = nullptr;
    G4LogicalVolume* fArgonOuter_log = nullptr;
    G4LogicalVolume* fInnerFrame_log = nullptr;
    G4LogicalVolume* fReflectorFoil_log = nullptr;
    G4LogicalVolume* fTPBFoilSides_log = nullptr;
    G4LogicalVolume* fTPBFoilTop_log = nullptr;
    G4LogicalVolume* fTPBFoilBottom_log = nullptr;
    G4LogicalVolume* fFiducialAr_log = nullptr;

    G4LogicalVolume* fPMTCoatedWall_log = nullptr;
    G4LogicalVolume* fPMTCoatedCaps_log = nullptr;
    G4LogicalVolume* fPMTUncoatedWall_log = nullptr;
    G4LogicalVolume* fPMTUncoatedCaps_log = nullptr;
    
    G4LogicalVolume* fBridleWall_log = nullptr;
    G4LogicalVolume* fBridleCaps_log = nullptr;
    G4LogicalVolume* fFrillWall_log = nullptr;
    G4LogicalVolume* fFrillCaps_log = nullptr;

    G4LogicalVolume* fTPBCoatingWall_log = nullptr; 
    G4LogicalVolume* fTPBCoatingCaps_log = nullptr; 

    G4LogicalVolume* fPhotocathCoated_log = nullptr;
    G4LogicalVolume* fPhotocathUncoated_log = nullptr;
    G4LogicalVolume* fSourceRod_log = nullptr;
    G4LogicalVolume* fSourcePellet_log = nullptr;
    
    G4LogicalVolume* fShinyC406R0_log = nullptr;
    G4LogicalVolume* fShinyTop_log = nullptr;
    G4LogicalVolume* fShinyBottom_log = nullptr;

    // Physical volumes (necessary for borders)
    G4VPhysicalVolume* fFiducialAr_phys = nullptr;
    G4VPhysicalVolume* fTPBFoilSides_phys = nullptr;
    G4VPhysicalVolume* fTPBFoilTop_phys = nullptr;
    G4VPhysicalVolume* fTPBFoilBottom_phys = nullptr;
    G4VPhysicalVolume* fTPBPMT_phys = nullptr;
    G4VPhysicalVolume* fReflectorFoil_phys = nullptr;

    // Sensitive Detectors positions
    std::vector<G4ThreeVector> fPMTPositions;

};

#endif

