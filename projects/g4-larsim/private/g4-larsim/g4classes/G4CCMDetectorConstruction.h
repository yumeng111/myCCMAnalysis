/*
  header file for the detector construction

initiates the methods used in G4CCMDetectorConstruction.cc and defines a number of global variables
Methods include: geometry modifications, pmt placement, material definition, and boolean swtichers
Global variables include most of the important volumes and materials.
*/
#ifndef G4CCMDetectorConstruction_H
#define G4CCMDetectorConstruction_H 1

#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "g4-larsim/g4classes/G4CCMReadout.h"
#include "g4-larsim/g4classes/G4CCMTreeTracker.h"
#include "g4-larsim/g4classes/G4CCMEDepSD.h"
#include "simclasses/CCMSimulationSettings.h"
#include "simclasses/DetectorResponseConfig.h"

#include <G4Cache.hh>
#include <G4Material.hh>
#include <G4RunManager.hh>
#include <G4VisAttributes.hh>
#include <G4RotationMatrix.hh>
#include <CLHEP/Units/SystemOfUnits.h>
#include <G4VUserDetectorConstruction.hh>
#include <G4SystemOfUnits.hh>

class G4CCMMainVolume;
class G4CCMPMTSD;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;

class G4CCMDetectorConstruction : public G4VUserDetectorConstruction {
  public:

    // constructor and destructor
    G4CCMDetectorConstruction(CCMSimulationSettings const & settings, DetectorResponseConfig const & config);
    ~G4CCMDetectorConstruction() override;

    void SetReadout(G4CCMReadout * readout);

    // build the detector
    G4VPhysicalVolume* Construct() override;

    // Add SD
    void ConstructSDandField() override;

    //Functions to modify the geometry
    void SetDefaults();//Method to set the default values of all geometry modifications

    //Construct main volume
    void SetMainVolumeOn(G4bool b);
    G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

    double HarmonicOscillatorRefractiveIndex(double a0, double aUV, double gammaUV, double wlUV, double wl);
    double HarmonicOscillatorRefractiveIndexDerivative(double aUV, double lambda, double lambdaUV, double gammaUV);

  private:

    G4CCMReadout * readout_ = nullptr;

    void DefineMaterials();
    G4CCMDetectorMessenger* fDetectorMessenger = nullptr;

    G4Box* fExperimentalHall_box = nullptr;
    G4LogicalVolume* fExperimentalHall_log = nullptr;
    G4VPhysicalVolume* fExperimentalHall_phys = nullptr;

    // Materials & Elements
    G4Element* fH = nullptr;
    G4Element* fC = nullptr;
    G4Element* fN = nullptr;
    G4Element* fO = nullptr;
    G4Element* elFe = nullptr;
    G4Element* elCr = nullptr;
    G4Element* elNi = nullptr;
    G4Element* elC = nullptr;
    G4Element* elMn = nullptr;
    G4Element* elF = nullptr;
    G4Material* fGlass = nullptr;
    G4Material* fPlastic = nullptr;
    G4Material* fBlackPlastic = nullptr;
    G4Material* fLAr = nullptr;
    G4Material* fAlum = nullptr;
    G4Material* fSteel = nullptr;
    G4Material* fVacuum = nullptr;
    G4Material* fPTFE = nullptr;
    G4Material* fTPBFoilSides = nullptr;
    G4Material* fTPBFoilTopBottom = nullptr;
    G4Material* fTPBPMT = nullptr;

    // Geometry
    G4bool fMainVolumeOn = true;

    G4CCMMainVolume* fMainVolume = nullptr;

    G4MaterialPropertiesTable* fGlass_mt = nullptr;
    G4MaterialPropertiesTable* fPlastic_mt = nullptr;
    G4MaterialPropertiesTable* fBlackPlastic_mt = nullptr;
    G4MaterialPropertiesTable* fLAr_mt = nullptr;
    G4MaterialPropertiesTable* fAlum_mt = nullptr;
    G4MaterialPropertiesTable* fSteel_mt = nullptr;
    G4MaterialPropertiesTable* fVacuum_mt = nullptr;
    G4MaterialPropertiesTable* fPTFE_mt = nullptr;
    G4MaterialPropertiesTable* fTPBFoilSides_mt = nullptr;
    G4MaterialPropertiesTable* fTPBFoilTopBottom_mt = nullptr;
    G4MaterialPropertiesTable* fTPBPMT_mt = nullptr;

    // Sensitive Detector
    G4Cache<G4CCMPMTSD*> fPMT_SD;
    G4Cache<G4CCMTreeTracker*> fTreeTracker_SD;
    G4Cache<G4CCMEDepSD*> fVeto_SD;
    G4Cache<G4CCMEDepSD*> fInterior_SD;

    // Settings
    CCMSimulationSettings simulationSettings_;
    DetectorResponseConfig detectorConfig_;
};

#endif
