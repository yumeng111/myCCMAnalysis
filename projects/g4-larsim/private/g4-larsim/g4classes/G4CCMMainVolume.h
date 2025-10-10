#ifndef G4CCMMainVolume_h
#define G4CCMMainVolume_h 1

#include "g4-larsim/g4classes/J4PMTSolidMaker.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

#include "simclasses/CCMSimulationSettings.h"
#include "simclasses/DetectorResponseConfig.h"

#include <map>

#include <G4PVPlacement.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>

class G4Box;
class G4LogicalVolume;
class G4Sphere;
class G4Tubs;

class G4CCMMainVolume : public G4PVPlacement {
    std::vector<std::shared_ptr<G4RotationMatrix>> placement_rotations;
  public:
    G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                G4LogicalVolume* pMotherLogical, G4bool pMany, G4int pCopyNo,
                G4CCMDetectorConstruction* c,
                CCMSimulationSettings const & settings,
                DetectorResponseConfig const & config);

    G4LogicalVolume* GetLogPMTCoatedWall() { return fPMTCoatedWall_log; }
    G4LogicalVolume* GetLogPMTCoatedCaps() { return fPMTCoatedCaps_log; }
    G4LogicalVolume* GetLogPMTUncoatedWall() { return fPMTUncoatedWall_log; }
    G4LogicalVolume* GetLogPMTUncoatedCaps() { return fPMTUncoatedCaps_log; }

    G4LogicalVolume* GetLogScint() { return fFiducialLAr_log; }
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

    G4LogicalVolume* GetLogInnerFrame() {return fInnerFrame_log;}
    G4LogicalVolume* GetLogArgonOuter() {return fOuterLAr_log;}

    G4LogicalVolume* GetLogInnerJacket() {return fInnerJacket_log;}
    G4LogicalVolume* GetLogVacuum() {return fVacuum_log;}
    G4LogicalVolume* GetLogCryoVessel() {return fCryoVessel_log;}

    std::vector<G4LogicalVolume*> GetPMTLogicalVolumes() {
        std::vector<G4LogicalVolume*> res = {fPMTCoatedWall_log, fPMTCoatedCaps_log, fPMTUncoatedWall_log, fPMTUncoatedCaps_log};
        for(std::pair<std::tuple<bool, G4VSolid *> const, std::shared_ptr<G4LogicalVolume>> const & p : fPMTLogicalVolumes) {
            res.push_back(p.second.get());
        }
        return res;
    }

    std::vector<G4LogicalVolume*> GetScintLogicalVolumes() {
        return {fFiducialLAr_log, fOuterLAr_log};
    }

    std::vector<G4LogicalVolume*> GetSourceLogicalVolumes() {
        return {fSourceRod_log, fSourcePellet_log, fSourcePelletHousing_log};
    }

    std::vector<G4LogicalVolume*> GetTPBLogicalVolumes() {
        std::vector<G4LogicalVolume*> res = {fTPBCoatingWall_log, fTPBCoatingCaps_log, fTPBFoilSides_log, fTPBFoilTop_log, fTPBFoilBottom_log};
        for(auto const & p : fTPBLogicalVolumes) {
            res.push_back(p.second.get());
        }
        return res;
    }

    std::vector<G4LogicalVolume*> GetPMTVacuumLogicalVolumes() {
        std::vector<G4LogicalVolume*> res = {fPMTVacuumWall_log, fPMTVacuumCaps_log};
        for(auto const & p : fVacuumLogicalVolumes) {
            res.push_back(p.second.get());
        }
        return res;
    }

    std::vector<G4LogicalVolume*> GetReflectorLogicalVolumes() {
        return {fReflectorFoil_log, fShinyC406R0_log, fShinyTop_log, fShinyBottom_log};
    }

    std::vector<G4LogicalVolume*> GetBridleLogicalVolumes() {
        std::vector<G4LogicalVolume * > logicalVolumes;
        for(auto logicalVolume : fBridleLogicalVolumes) {
            logicalVolumes.push_back(logicalVolume.second.get());
        }
        return logicalVolumes;
    }

    std::vector<G4LogicalVolume*> GetFrillLogicalVolumes() {
        std::vector<G4LogicalVolume * > logicalVolumes;
        for(auto logicalVolume : fFrillLogicalVolumes) {
            logicalVolumes.push_back(logicalVolume.second.get());
        }
        return logicalVolumes;
    }

    std::vector<G4LogicalVolume*> GetFrameLogicalVolumes() {
        return {fInnerFrame_log};
    }

    std::vector<G4LogicalVolume*> GetCryoLogicalVolumes() {
        return {fInnerJacket_log, fVacuum_log, fCryoVessel_log};
    }

    std::vector<G4LogicalVolume*> GetAllLogicalVolumes() {
        std::vector<G4LogicalVolume*> allLogicalVolumes;
        allLogicalVolumes.push_back(fMother_log);
        std::vector<std::vector<G4LogicalVolume*>> logicalVolumes = {GetPMTLogicalVolumes(), GetScintLogicalVolumes(), GetSourceLogicalVolumes(), GetTPBLogicalVolumes(), GetPMTVacuumLogicalVolumes(), GetReflectorLogicalVolumes(), GetBridleLogicalVolumes(), GetFrillLogicalVolumes(), GetFrameLogicalVolumes(), GetCryoLogicalVolumes()};
        for(auto logicalVolume : logicalVolumes) {
            allLogicalVolumes.insert(allLogicalVolumes.end(), logicalVolume.begin(), logicalVolume.end());
        }
        return allLogicalVolumes;
    }

    std::vector<G4LogicalVolume*> GetInteriorLArLogicalVolumes() {
        return {fFiducialLAr_log};
    }

    std::vector<G4LogicalVolume*> GetVetoLArLogicalVolumes() {
        return {fOuterLAr_log};
    }

    std::vector<G4ThreeVector> GetPMTPositions() { return fPMTPositions; }

  private:

    J4PMTSolidMaker pmt_solid_maker;

    void VisAttributes(G4bool SourceRodIn);
    void SurfaceProperties();
    G4CCMDetectorConstruction* fConstructor = nullptr;

    DetectorResponseConfig detectorConfig_;
    CCMSimulationSettings simulationSettings_;

    // Basic Volumes
    G4Tubs* fCryoVessel = nullptr;
    G4Tubs* fVacuum = nullptr;
    G4Tubs* fInnerJacket = nullptr;
    G4Tubs* fOuterLAr = nullptr;
    G4Tubs* fInnerFrame = nullptr;
    G4Tubs* fReflectorFoil = nullptr;
    G4Tubs* fTPBFoilSides = nullptr;
    G4Tubs* fTPBFoilTop = nullptr;
    G4Tubs* fTPBFoilBottom = nullptr;
    G4Tubs* fFiducialLAr = nullptr;

    G4VSolid* fPMTVacuumWall = nullptr;
    G4VSolid* fPMTVacuumCaps = nullptr;

    G4VSolid* fPMTCoatedWall = nullptr;
    G4VSolid* fPMTCoatedCaps = nullptr;
    G4VSolid* fPMTUncoatedWall = nullptr;
    G4VSolid* fPMTUncoatedCaps = nullptr;

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

    std::map<G4VSolid *, std::shared_ptr<G4LogicalVolume>> fBridleLogicalVolumes;
    std::map<G4VSolid *, std::shared_ptr<G4LogicalVolume>> fFrillLogicalVolumes;
    std::map<G4VSolid *, std::shared_ptr<G4LogicalVolume>> fTPBLogicalVolumes;
    std::map<G4VSolid *, std::shared_ptr<G4LogicalVolume>> fVacuumLogicalVolumes;
    std::map<std::tuple<bool, G4VSolid *>, std::shared_ptr<G4LogicalVolume>> fPMTLogicalVolumes;
    std::map<std::tuple<bool, G4VSolid *>, std::shared_ptr<G4LogicalSkinSurface>> fLogicalSkinSurfaces;
    std::vector<std::shared_ptr<G4LogicalBorderSurface>> fLogicalBorderSurfaces;
    std::vector<std::shared_ptr<G4PVPlacement>> fPlacements;

    // Logical volumes
    G4LogicalVolume* fMother_log = nullptr;
    G4LogicalVolume* fCryoVessel_log = nullptr;
    G4LogicalVolume* fVacuum_log = nullptr;
    G4LogicalVolume* fInnerJacket_log = nullptr;
    G4LogicalVolume* fOuterLAr_log = nullptr;
    G4LogicalVolume* fInnerFrame_log = nullptr;
    G4LogicalVolume* fReflectorFoil_log = nullptr;
    G4LogicalVolume* fTPBFoilSides_log = nullptr;
    G4LogicalVolume* fTPBFoilTop_log = nullptr;
    G4LogicalVolume* fTPBFoilBottom_log = nullptr;
    G4LogicalVolume* fFiducialLAr_log = nullptr;

    G4LogicalVolume* fPMTVacuumWall_log = nullptr;
    G4LogicalVolume* fPMTVacuumCaps_log = nullptr;

    G4LogicalVolume* fPMTCoatedWall_log = nullptr;
    G4LogicalVolume* fPMTCoatedCaps_log = nullptr;
    G4LogicalVolume* fPMTUncoatedWall_log = nullptr;
    G4LogicalVolume* fPMTUncoatedCaps_log = nullptr;

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
    G4VPhysicalVolume* fFiducialLAr_phys = nullptr;
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

