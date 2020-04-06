/*
  header file for the detector construction

initiates the methods used in detectorConstruction.cc and defines a number of global variables
Methods include: geometry modifications, pmt placement, material definition, and boolean swtichers
Global variables include most of the important volumes and materials.
*/
#ifndef detectorConstruction_H
#define detectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class G4Sphere;

#include "G4Material.hh"
#include "detectorMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"

class detectorConstruction : public G4VUserDetectorConstruction
{
  public:

    detectorConstruction();
    virtual ~detectorConstruction();
  //constructor and deconstructor

    virtual G4VPhysicalVolume* Construct();
  //Deffines the method that will actually build the detector

    //Functions to modify the geometry
    void SetDefaults();//Method to set the default values of all geometry modifications
    void SetPMTRadius(G4double );//defines the radius of the PMTs (currently 4 inches ~ 10.2 cm)
    void SetPMTsOn(G4bool );//Turns the PMTs on or off
    void SetSodiumOn(G4bool );//Turns the Sodium on or off
    void SetMainVolOn(G4bool );//Sets the main volume on or off
    void SetCCM200(G4bool );//sets the numberr of pmts to 200 (true) or 120 (false)
    void SetRodin(G4bool );//Sets the calibration rod to true or false
    void SetTPBfoilOn(G4bool );//Turns the TPB coating on the foils on or off. Also turns on the triple cylinder for the LBOC template.
    void SetReflectorOn(G4bool );//Turns the ptfe foil reflector on or off
    void SetRodHeight( G4double );//Sets the height of the Laser calibration rod. also turns the laser to true
    

    //Get values
    G4double GetPMTRadius() const {return fOuterRadius_pmt;}//method to obtain the PMT radiu entered in the previous method (or the default, if unchanged)
    G4double GetRodHeight() const {return rodHeight;};//method to obtain the rodHeight defined in the previous method
    G4bool GetfLaser() const { return fLaser;};//method to obtain fLaser
    G4bool GetfSodium() const {return fSodium;};//method to obtain fSodium
    
  private:

    void DefineMaterials();
    void placePMT( G4String , G4double , G4double , G4double , G4bool );
    void placeTopBot( G4String , G4double , G4double , G4double , G4bool );
  //methods to define materials and place PMTS on the sides or on the top/bottom

    detectorMessenger* fDetectorMessenger;
  //links the detector messenger to this detector construction

    //world geometries (line 264)
    G4Box* wBox;
    G4LogicalVolume* wLog;
    G4VPhysicalVolume* wPhys;

    //Materials & Elements (from line 100)
    G4Element* fH;
    G4Element* fC;
    G4Element* fN;
    G4Element* fO;
    G4Material* lAr;
    G4Material* lAr1;
    G4Material* lAr2;
    G4Material* alum;
    G4Material* steel;
    G4Material* fVacuum;
    G4Material* ptfe;
    G4Material* fGlass;
    G4Material* fAir;
    G4Material* ice;
    G4Material* tPB;
  //G4Material* tPBfoil;
    G4Material* tPBhundred;

    //Material Property Tables (from line 154)
    G4MaterialPropertiesTable* lAr_mt;
    G4MaterialPropertiesTable* TPBProp; //line 237
    G4MaterialPropertiesTable* TPBsProp; //line 237
  
    //Internal geomtries (from line 280)
    G4Tubs* fCryoVessel;
    G4LogicalVolume* fLogicCryo;
    G4Tubs* fArgonOuter;
    G4LogicalVolume* fLogicArout;
    G4Tubs* fInnerFrame;
    G4LogicalVolume* fLogicFrame;
    G4Tubs* fTPBSides;
    G4LogicalVolume* fLogicTPB;
    G4Tubs* fTPBBottom;
    G4LogicalVolume* fLogicTPBb;
    G4Tubs* fFiducialAr;
    G4LogicalVolume* fLogicFiduc;
    G4VPhysicalVolume* lArFiducial;
    G4Tubs* fFiducialAr1;
    G4LogicalVolume* fLogicFiduc1;
    G4VPhysicalVolume* lArFiducial1;
    G4Tubs* fFiducialAr2;
    G4LogicalVolume* fLogicFiduc2;
    G4VPhysicalVolume* lArFiducial2;
    G4Tubs* fFiducialAr3;
    G4LogicalVolume* fLogicFiduc3;
    G4VPhysicalVolume* lArFiducial3;
    G4Tubs* fFiducialAr4;
    G4LogicalVolume* fLogicFiduc4;
    G4VPhysicalVolume* lArFiducial4;
    G4Tubs* fFiducialAr5;
    G4LogicalVolume* fLogicFiduc5;
    G4VPhysicalVolume* lArFiducial5;
   
    //Geometry and boolean switches
    G4double fOuterRadius_pmt;
    
    G4bool fPMTsOn;
    G4bool fRodin;
    G4bool fSodium;
    G4bool fTPBfoilOn;
    G4bool fReflectorOn;
    G4bool mainVolOn;
    G4bool cylinderOn;
    G4bool fLaser;
    G4bool fLayers;
    G4bool ccm200;
    G4double rodHeight;

};

#endif
