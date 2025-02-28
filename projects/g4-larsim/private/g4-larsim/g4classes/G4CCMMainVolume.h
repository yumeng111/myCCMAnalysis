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
                G4CCMDetectorConstruction* c, G4bool SourceRodIn, G4double SourceRodLocation, G4bool CobaltSourceRun, G4bool SodiumSourceRun,
                G4bool TrainingSource, G4double DecayX, G4double DecayY, G4double DecayZ,
                G4double EndCapFoilTPBThickness, G4double SideFoilTPBThickness, G4double PMTTPBThickness);

    G4LogicalVolume* GetLogPMTCoatedWall() { return fPMTCoatedWall_log; }
    G4LogicalVolume* GetLogPMTCoatedCaps() { return fPMTCoatedCaps_log; }
    G4LogicalVolume* GetLogPMTUncoatedWall() { return fPMTUncoatedWall_log; }
    G4LogicalVolume* GetLogPMTUncoatedCaps() { return fPMTUncoatedCaps_log; }

    G4LogicalVolume* GetLogScint() { return fFiducialAr_log; }
    G4LogicalVolume* GetLogSourceRod() { return fSourceRod_log; }
    G4LogicalVolume* GetLogSourcePellet() { return fSourcePellet_log; }
    G4LogicalVolume* GetLogSourcePelletHousing() { return fSourcePelletHousing_log; }

    G4LogicalVolume* GetLogTPBCoatingWall() { return fTPBCoatingWall_log; }
    G4LogicalVolume* GetLogTPBCoatingCaps() { return fTPBCoatingCaps_log; }

    G4LogicalVolume* GetLogTPBFoilSides() { return fTPBFoilSides_log; }
    G4LogicalVolume* GetLogTPBFoilTop() { return fTPBFoilTop_log; }
    G4LogicalVolume* GetLogTPBFoilBottom() { return fTPBFoilBottom_log; }
    G4LogicalVolume* GetLogReflectorFoil() { return fReflectorFoil_log; }
    G4LogicalVolume* GetShinyC406R0() { return fShinyC406R0_log; }
    G4LogicalVolume* GetShinyTop() { return fShinyTop_log; }
    G4LogicalVolume* GetShinyBottom() { return fShinyBottom_log; }
    
    G4LogicalVolume* GetLogBridleCoatedWall() {return fBridleCoatedWall_log;}
    G4LogicalVolume* GetLogBridleUncoatedWall() {return fBridleUncoatedWall_log;}
    G4LogicalVolume* GetLogBridleCoatedCaps() {return fBridleCoatedCaps_log;}
    G4LogicalVolume* GetLogBridleUncoatedCaps() {return fBridleUncoatedCaps_log;}
    G4LogicalVolume* GetLogFrillCoatedWall() {return fFrillCoatedWall_log;}
    G4LogicalVolume* GetLogFrillUncoatedWall() {return fFrillUncoatedWall_log;}
    G4LogicalVolume* GetLogFrillCoatedCaps() {return fFrillCoatedCaps_log;}
    G4LogicalVolume* GetLogFrillUncoatedCaps() {return fFrillUncoatedCaps_log;}

    G4LogicalVolume* GetLogInnerFrame() {return fInnerFrame_log;}
    G4LogicalVolume* GetLogArgonOuter() {return fArgonOuter_log;}

    G4LogicalVolume* GetLogInnerJacket() {return fInnerJacket_log;}
    G4LogicalVolume* GetLogVacuum() {return fVacuum_log;}
    G4LogicalVolume* GetLogCryoVessel() {return fCryoVessel_log;}

    std::vector<G4LogicalVolume*> GetPMTLogicalVolumes() {
        return {fPMTCoatedWall_log, fPMTCoatedCaps_log, fPMTUncoatedWall_log, fPMTUncoatedCaps_log};
    }

    std::vector<G4LogicalVolume*> GetScintLogicalVolumes() {
        return {fFiducialAr_log, fArgonOuter_log};
    }

    std::vector<G4LogicalVolume*> GetSourceLogicalVolumes() {
        return {fSourceRod_log, fSourcePellet_log, fSourcePelletHousing_log};
    }

    std::vector<G4LogicalVolume*> GetTPBLogicalVolumes() {
        return {fTPBCoatingWall_log, fTPBCoatingCaps_log, fTPBFoilSides_log, fTPBFoilTop_log, fTPBFoilBottom_log};
    }

    std::vector<G4LogicalVolume*> GetReflectorLogicalVolumes() {
        return {fReflectorFoil_log, fShinyC406R0_log, fShinyTop_log, fShinyBottom_log};
    }

    std::vector<G4LogicalVolume*> GetBridleLogicalVolumes() {
        return {fBridleCoatedWall_log, fBridleUncoatedWall_log, fBridleCoatedCaps_log, fBridleUncoatedCaps_log};
    }

    std::vector<G4LogicalVolume*> GetFrillLogicalVolumes() {
        return {fFrillCoatedWall_log, fFrillUncoatedWall_log, fFrillCoatedCaps_log, fFrillUncoatedCaps_log};
    }

    std::vector<G4LogicalVolume*> GetFrameLogicalVolumes() {
        return {fInnerFrame_log};
    }

    std::vector<G4LogicalVolume*> GetCryoLogicalVolumes() {
        return {fInnerJacket_log, fVacuum_log, fCryoVessel_log};
    }

    std::vector<G4LogicalVolume*> GetAllLogicalVolumes() {
        std::vector<G4LogicalVolume*> allLogicalVolumes;
        std::vector<std::vector<G4LogicalVolume*>> logicalVolumes = {GetPMTLogicalVolumes(), GetSourceLogicalVolumes(), GetTPBLogicalVolumes(), GetReflectorLogicalVolumes(), GetBridleLogicalVolumes(), GetFrillLogicalVolumes(), GetFrameLogicalVolumes(), GetCryoLogicalVolumes()};
        for(auto logicalVolume : logicalVolumes) {
            allLogicalVolumes.insert(allLogicalVolumes.end(), logicalVolume.begin(), logicalVolume.end());
        }
        return allLogicalVolumes;
    }

    std::vector<G4LogicalVolume*> GetInteriorLArLogicalVolumes() {
        return {fFiducialAr_log};
    }

    std::vector<G4LogicalVolume*> GetVetoLArLogicalVolumes() {
        return {fOuterLAr_log};
    }

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
    G4Tubs* fTPBFoilTop = nullptr;
    G4Tubs* fTPBFoilBottom = nullptr;
    G4Tubs* fFiducialAr = nullptr;

    G4VSolid* fPMTCoatedWall = nullptr;
    G4VSolid* fPMTCoatedCaps = nullptr;
    G4VSolid* fPMTUncoatedWall = nullptr;
    G4VSolid* fPMTUncoatedCaps = nullptr;

    G4VSolid* fBridleCoatedWall = nullptr;
    G4VSolid* fBridleCoatedCaps = nullptr;
    G4VSolid* fFrillCoatedWall = nullptr;
    G4VSolid* fFrillCoatedCaps = nullptr;
    G4VSolid* fBridleUncoatedWall = nullptr;
    G4VSolid* fBridleUncoatedCaps = nullptr;
    G4VSolid* fFrillUncoatedWall = nullptr;
    G4VSolid* fFrillUncoatedCaps = nullptr;

    G4VSolid* fTPBCoatingWall = nullptr;
    G4VSolid* fTPBCoatingCaps = nullptr;

    G4VSolid* fPhotocathCoated = nullptr;
    G4VSolid* fPhotocathUncoated = nullptr;
    G4Tubs* fSourceRod = nullptr;
    G4Tubs* fSourcePellet = nullptr;
    G4Tubs* fSourcePelletHousing = nullptr;

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

    G4LogicalVolume* fBridleCoatedWall_log = nullptr;
    G4LogicalVolume* fBridleCoatedCaps_log = nullptr;
    G4LogicalVolume* fFrillCoatedWall_log = nullptr;
    G4LogicalVolume* fFrillCoatedCaps_log = nullptr;
    G4LogicalVolume* fBridleUncoatedWall_log = nullptr;
    G4LogicalVolume* fBridleUncoatedCaps_log = nullptr;
    G4LogicalVolume* fFrillUncoatedWall_log = nullptr;
    G4LogicalVolume* fFrillUncoatedCaps_log = nullptr;

    G4LogicalVolume* fTPBCoatingWall_log = nullptr;
    G4LogicalVolume* fTPBCoatingCaps_log = nullptr;

    G4LogicalVolume* fPhotocathCoated_log = nullptr;
    G4LogicalVolume* fPhotocathUncoated_log = nullptr;
    G4LogicalVolume* fSourceRod_log = nullptr;
    G4LogicalVolume* fSourcePellet_log = nullptr;
    G4LogicalVolume* fSourcePelletHousing_log = nullptr;

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
    G4VPhysicalVolume* fPMTUncoatedWall_phys = nullptr;
    G4VPhysicalVolume* fPMTUncoatedCaps_phys = nullptr;
    G4VPhysicalVolume* fPMTCoatedWall_phys = nullptr;
    G4VPhysicalVolume* fPMTCoatedCaps_phys = nullptr;

    // Sensitive Detectors positions
    std::vector<G4ThreeVector> fPMTPositions;

};

#endif

