/*
  header file for the detector construction

initiates the methods used in G4CCMDetectorConstruction.cc and defines a number of global variables
Methods include: geometry modifications, pmt placement, material definition, and boolean swtichers
Global variables include most of the important volumes and materials.
*/
#ifndef G4CCMDetectorConstruction_H
#define G4CCMDetectorConstruction_H 1

#include "icetray/CCMPMTKey.h"
#include "dataclasses/I3Map.h"
#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"

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
class G4CCMScintSD;

class G4Box;
class G4Element;
class G4LogicalVolume;
class G4Material;
class G4MaterialPropertiesTable;
class G4Sphere;
class G4Tubs;
class G4VPhysicalVolume;

class G4CCMDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    // constructor and destructor
    G4CCMDetectorConstruction(G4double SingletTau, G4double TripletTau, G4bool UVAbsStatus, G4double Rayleigh128);
    ~G4CCMDetectorConstruction() override;

    // build the detector
    G4VPhysicalVolume* Construct() override;

    // Add SD
    void ConstructSDandField() override;

    //Functions to modify the geometry
    void SetDefaults();//Method to set the default values of all geometry modifications
    
    //Construct main volume
    void SetMainVolumeOn(G4bool b);
    G4bool GetMainVolumeOn() const { return fMainVolumeOn; }
   
    // set SD configuration
    void SetPMTSDStatus(bool PMTSDStatus) { PMTSDStatus_ = PMTSDStatus; }
    void SetLArSDStatus(bool LArSDStatus) { LArSDStatus_ = LArSDStatus; }
    
    // get SD configuration
    bool GetPMTSDStatus() { return PMTSDStatus_; }
    bool GetLArSDStatus() { return LArSDStatus_; }
 
    // set sodium source calibration run status
    void InitializeSodiumSourceRun(bool SodiumSourceRun, G4double SodiumSourceLocation){
        SodiumSourceRun_ = SodiumSourceRun; 
        SodiumSourceLocation_ = SodiumSourceLocation;
    }

    // set time cut
    void SetTimeCut(bool TimeCut){ TimeCut_ = TimeCut; }
    
    // set cerenkov control
    void SetCerenkovControl(bool CerenkovControl) { CerenkovControl_ = CerenkovControl; }    

    // get sodium source calibration run status
    bool GetSodiumSourceRun() { return SodiumSourceRun_; }
    G4double GetSodiumSourcePosition() { return SodiumSourceLocation_; }

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
    G4Material* fTPBFoil;
    G4Material* fTPBPMT;
    
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
    G4MaterialPropertiesTable* fTPBFoil_mt = nullptr;
    G4MaterialPropertiesTable* fTPBPMT_mt = nullptr;

    // Sensitive Detector
    G4Cache<G4CCMPMTSD*> fPMT_SD;
    G4Cache<G4CCMScintSD*> fScint_SD;

    // controls to turn SD on/off (set via G4Interface)
    bool PMTSDStatus_ = true; // turn PMT SD on/off
    bool LArSDStatus_ = true; // turn fiducial LAr SD on/off

    // controls to turn sodium source on/off
    bool SodiumSourceRun_ = false;
    G4double SodiumSourceLocation_ = 0.0; 

    G4double SingletTau_ = 8.2 * ns;
    G4double TripletTau_ = 743.0 * ns;
    G4double Rayleigh128_ = 95.0 * cm;
    G4bool UVAbsStatus_ = true;

    G4bool TimeCut_ = true;
    G4bool CerenkovControl_ = false;
};

#endif
