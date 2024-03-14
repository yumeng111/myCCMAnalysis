/*
  header file for the detector construction

initiates the methods used in G4CCMDetectorConstruction.cc and defines a number of global variables
Methods include: geometry modifications, pmt placement, material definition, and boolean swtichers
Global variables include most of the important volumes and materials.
*/
#ifndef G4CCMDetectorConstruction_H
#define G4CCMDetectorConstruction_H 1

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


#include "G4Material.hh"
#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include <CLHEP/Units/SystemOfUnits.h>

class G4CCMDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    // constructor and destructor
    G4CCMDetectorConstruction();
    ~G4CCMDetectorConstruction() override;

    // build the detector
    G4VPhysicalVolume* Construct() override;

    // Add SD
    void ConstructSD();

    //Functions to modify the geometry
    void SetDefaults();//Method to set the default values of all geometry modifications

    
    //Construct main volume
    void SetMainVolumeOn(G4bool b);
    G4bool GetMainVolumeOn() const { return fMainVolumeOn; }

  private:
    
    void DefineMaterials();

    G4CCMDetectorMessenger* fDetectorMessenger = nullptr;

    G4Box* fExperimentalHall_box = nullptr;
    G4LogicalVolume* fExperimentalHall_log = nullptr;
    G4VPhysicalVolume* fExperimentalHall_phys = nullptr;

    // Materials & Elements
    G4Element* fH;
    G4Element* fC;
    G4Element* fN;
    G4Element* fO;
    G4Element* elFe;
    G4Element* elCr;
    G4Element* elNi;
    G4Element* elC;
    G4Element* elMn;
    G4Element* elF;
    G4Material* fGlass;
    G4Material* fPlastic;
    G4Material* fBlackPlastic;
    G4Material* fLAr;
    G4Material* fAlum;
    G4Material* fSteel;
    G4Material* fVacuum;
    G4Material* fPTFE;
    G4Material* fTPBCoating;
    G4Material* fTPBFoil;

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
    G4MaterialPropertiesTable* fTPBCoating_mt = nullptr;
    G4MaterialPropertiesTable* fTPBFoil_mt = nullptr;

    // Sensitive Detector
    G4Cache<G4CCMPMTSD*> fPMT_SD;

};

#endif
