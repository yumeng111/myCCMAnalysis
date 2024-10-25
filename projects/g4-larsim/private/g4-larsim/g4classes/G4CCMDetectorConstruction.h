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
#include "g4-larsim/g4classes/G4CCMReadout.h"

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

class G4CCMDetectorConstruction : public G4VUserDetectorConstruction {
  public:

    // constructor and destructor
    G4CCMDetectorConstruction(G4double SingletTau, G4double TripletTau, G4double UVAbsLength,
                              G4double WLSNPhotonsEndCapFoil, G4double WLSNPhotonsSideFoil, G4double WLSNPhotonsPMT,
                              G4double EndCapFoilTPBThickness, G4double SideFoilTPBThickness, G4double PMTTPBThickness,
                              G4double Rayleigh128, G4double TPBAbsTau, G4double TPBAbsNorm, G4double TPBAbsScale,
                              G4double Mie_GG, G4double Mie_Ratio, G4double Normalization);
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

    // set SD configuration
    void SetPMTSDStatus(bool PMTSDStatus) { PMTSDStatus_ = PMTSDStatus; }
    void SetLArSDStatus(bool LArSDStatus) { LArSDStatus_ = LArSDStatus; }

    // get SD configuration
    bool GetPMTSDStatus() { return PMTSDStatus_; }
    bool GetLArSDStatus() { return LArSDStatus_; }

    // set sodium source calibration run status
    void InitializeSodiumSourceRun(bool SourceRodIn, G4double SourceRodLocation, bool CobaltSourceRun, bool SodiumSourceRun) {
        SourceRodIn_ = SourceRodIn;
        SourceRodLocation_ = SourceRodLocation;
        CobaltSourceRun_ = CobaltSourceRun;
        SodiumSourceRun_ = SodiumSourceRun;
    }

    // set time cut
    void SetTimeCut(bool TimeCut){ TimeCut_ = TimeCut; }

    // set cerenkov control
    void SetKillCherenkov(bool KillCherenkov) { KillCherenkov_ = KillCherenkov; }

    void SetPhotonTracking(bool FullPhotonTracking) { FullPhotonTracking_ = FullPhotonTracking; }

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
    G4Cache<G4CCMScintSD*> fScint_SD;

    // controls to turn SD on/off (set via G4Interface)
    bool PMTSDStatus_ = true; // turn PMT SD on/off
    bool LArSDStatus_ = true; // turn fiducial LAr SD on/off

    // controls to turn sodium source on/off
    bool SourceRodIn_ = false;
    G4double SourceRodLocation_ = 0.0 * cm;
    bool CobaltSourceRun_ = false;
    bool SodiumSourceRun_ = false;

    G4double SingletTau_ = 8.2 * ns;
    G4double TripletTau_ = 743.0 * ns;
    G4double Rayleigh128_ = 95.0 * cm;
    G4double UVAbsLength_ = 55.0 * cm;
    G4double WLSNPhotonsEndCapFoil_ = 0.605;
    G4double WLSNPhotonsSideFoil_ = 0.605;
    G4double WLSNPhotonsPMT_ = 0.605;
    G4double EndCapFoilTPBThickness_ = 0.00278035 * mm;
    G4double SideFoilTPBThickness_ = 0.00278035 * mm;
    G4double PMTTPBThickness_ = 0.00203892 * mm;
    G4double TPBAbsTau_ = 0.13457;
    G4double TPBAbsNorm_ = 8.13914e-21;
    G4double TPBAbsScale_ = 1.0;
    G4double Mie_GG_ = 0.99;
    G4double Mie_Ratio_ = 1.0;
    G4double Normalization_ = 1.0;

    G4bool TimeCut_ = true;
    G4bool KillCherenkov_ = false;
    G4bool FullPhotonTracking_ = true;
};

#endif
