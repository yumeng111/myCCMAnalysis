/*
  Best Fit (to Laser, Cobalt, and Sodium) Detector Construction for CCM simulation

most of the acctual simulation is done here. 
This code builds the detector, including the optical physics and materials.
The code also includes a number of methods for modifying the detector from macros (see detectorMessenger.cc)
This code is very long, and should be split into a number of seperate portions each representing a different configuration.

This version is the current best fit for the calibration data I have. It includes the layered detector that clouds the outer edges of the upper portions of the Fiducial volume, the lower absorption lengths near the bottom, the bottom 'foil cone' that helps match data and represents the LED and foil distortions in the area, the less effective TPB for the foils as opposed to the PMTs. It also uses a number of tricks that the other versions don't to make it easier to change a number of variables without messing up their relations (as the best fit is not yet a perfect fit, I am continuing to improve it constantly and thus this version will not likely stay the best for long). 
 */

#include "CCMAnalysis/CCMDetectorSimulation/detectorConstruction.hh"
#include "CCMAnalysis/CCMDetectorSimulation/detectorMessenger.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//Constructor
detectorConstruction::detectorConstruction()
  : lAr_mt(nullptr), lAr2_mt(nullptr), TPBProp(nullptr), TPBsProp(nullptr)
{
  wBox = nullptr;
  wLog = nullptr;
  wPhys = nullptr;
  
  fH = fC = fN = fO = nullptr;
  alum = steel = fVacuum = fEj335 = nullptr;
  ptfe = fGlass = tPB = tPBhundred = fBlackPlastic = fPlastic = nullptr;
  lAr = lAr2 = nullptr;

  fCryoVessel = nullptr;
  fLogicCryo = nullptr;
  fArgonOuter = nullptr;
  fLogicArout = nullptr;
  fInnerFrame = nullptr;
  fLogicFrame = nullptr;
  fTPBSides = nullptr;
  fLogicTPB = nullptr;
  fFiducialAr = nullptr;
  fLogicFiduc = nullptr;
  lArFiducial = nullptr;
  fFiducialAr1 = nullptr;
  fLogicFiduc1 = nullptr;
  lArFiducial1 = nullptr;
  fFiducialAr2 = nullptr;
  fLogicFiduc2 = nullptr;
  lArFiducial2 = nullptr;
  fFiducialAr3 = nullptr;
  fLogicFiduc3 = nullptr;
  lArFiducial3 = nullptr;
  fFiducialAr4 = nullptr;
  fLogicFiduc4 = nullptr;
  lArFiducial4 = nullptr;
  fFiducialAr5 = nullptr;
  fLogicFiduc5 = nullptr;
  lArFiducial5 = nullptr;
  //A long series of nullptr's to clear all the pointers created in the header file.

  SetDefaults();
  //Sets a few defaults by calling a method that declears them below.

  DefineMaterials();

  fDetectorMessenger = new detectorMessenger(this);
  //Calls the methods to define the materials and link the detector messenger to this detector.
}

//Deconstructor
detectorConstruction::~detectorConstruction() {}

//Defines the materials and many of the physics
//Most notably defines TPB wavelength shifting and the lAR absorption lengths.
void detectorConstruction::DefineMaterials(){
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  //defines a few doubles for use later.

  //Define the elements used in compounds.
  fH = new G4Element("H", "H", z=1., a=1.01*g/mole);
  fC = new G4Element("C", "C", z=6., a=12.01*g/mole);
  fN = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  fO = new G4Element("O", "O", z=8., a= 16.00*g/mole);
  //G4Element* fGd = new G4Element("Gd", "Gd", z=64., a=155.0*g/mole);
  
  //Aluminum: for the frame
  alum = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);

  //Steel: for the cryogen
  G4NistManager* manager = G4NistManager::Instance();
  steel = manager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  //ptfe: for the reflector foils.
  ptfe = manager->FindOrBuildMaterial("G4_TEFLON");
  G4Element* fGd = manager->FindOrBuildElement(64);
  //Vacuum: for vacuum. 
  fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                          density=universe_mean_density,kStateGas,0.1*kelvin,
                          1.e-19*pascal);
  //Glass: for PMTs.
  fGlass = new G4Material("Glass", density=2.5*g/cm3,2);
  fGlass->AddElement(fC,91.533*perCent);
  fGlass->AddElement(fH,8.467*perCent);
  //TPB: TPB, wavelength shifting material (mostly defined in optical properties)
  tPB = new G4Material("tpb", density= 1.079*g/cm3, 2);
  tPB->AddElement(fC, 28);
  tPB->AddElement(fH, 22);
  //second kind of TPB for different optical properties.
  tPBhundred = new G4Material("t100", density= 1.079*g/cm3, 2);
  tPBhundred->AddElement(fC, 28);
  tPBhundred->AddElement(fH, 22);
  //plastic for PMT frill
  fPlastic = new G4Material("Plastic", density=1.20*g/cm3,3);
  fPlastic->AddElement(fC,15);
  fPlastic->AddElement(fH,16);
  fPlastic->AddElement(fO,2);
  //BlackPlastic for PMT frill
  fBlackPlastic = new G4Material("BlackPlastic", density=1.20*g/cm3,3);
  fBlackPlastic->AddElement(fC,15);
  fBlackPlastic->AddElement(fH,16);
  fBlackPlastic->AddElement(fO,2);

  //MARK FOR REACTOR
  //EJ-335 for liquid scintillator detector  
  fEj335 = new G4Material("ej-335",density=0.89*g/cm3,3);
  fEj335->AddElement(fH, 11.55*perCent);
  fEj335->AddElement(fC, 88.20*perCent);
  fEj335->AddElement(fGd, 0.25*perCent);

  //ej_wavelengths = {300, 380, 400, 420, 424, 436, 440, 446, 460, 480, 500, 520, 540, 600}
  G4double ej_energy[] =  {4.133*eV, 3.263*eV, 3.099*eV, 2.952*eV, 2.924*eV,
			   2.844*eV, 2.818*eV, 2.780*eV, 2.695*eV, 2.583*eV,
			   2.480*eV, 2.384*eV, 2.296*eV, 2.066*eV};
  G4double ej_scint_amp[] =  {0.0,  0.0,  0.0, 0.88, 1.0,
			      0.9,  0.75, 0.6, 0.44, 0.18,
			      0.06, 0.0,  0.0, 0.0};
  assert(sizeof(ej_energy) == sizeof(ej_scint_amp));
  G4double ej_rayleigh[] = {42.3*cm,  42.3*cm,  56.4*cm,  59.2*cm, 59.8*cm,
			    61.5*cm,  62.02*cm, 62.9*cm,  74.4*cm, 90.9*cm,
			    107.4*cm, 123.9*cm, 141.0*cm, 141.0*cm};
  assert(sizeof(ej_energy) == sizeof(ej_rayleigh));
  G4double ej_opAbsorb[] = {4.5*m, 4.5*m, 4.5*m, 4.5*m, 4.5*m,
			    4.5*m, 4.5*m, 4.5*m, 4.5*m, 4.5*m,
			    4.5*m, 4.5*m, 4.5*m, 4.5*m};
  assert(sizeof(ej_energy) == sizeof(ej_opAbsorb));
  G4double ej_rindex[] =  {1.49, 1.49, 1.49, 1.49, 1.49,
			   1.49, 1.49, 1.49, 1.49, 1.49,
			   1.49, 1.49, 1.49, 1.49};
  assert(sizeof(ej_energy) == sizeof(ej_rindex));

  const G4int ejnum = sizeof(ej_energy)/sizeof(G4double);

  G4MaterialPropertiesTable *ej335_mt = new G4MaterialPropertiesTable();
  ej335_mt->AddProperty("ABSLENGTH",ej_energy,ej_opAbsorb,ejnum);
  ej335_mt->AddProperty("RINDEX",ej_energy,ej_rindex,ejnum);
  ej335_mt->AddProperty("RAYLEIGH",ej_energy,ej_rayleigh,ejnum);
  ej335_mt->AddProperty("FASTCOMPONENT", ej_energy, ej_scint_amp, ejnum);
  G4double ej_scint_yeild=9.57/(1.0*keV); //scintillation yeild: 9.57/keV
  ej335_mt->AddConstProperty("SCINTILLATIONYIELD",ej_scint_yeild);
  ej335_mt->AddConstProperty("RESOLUTIONSCALE",0.11);
  ej335_mt->AddConstProperty("FASTTIMECONSTANT",3.8*ns);
  ej335_mt->AddConstProperty("YIELDRATIO",1.0);
  fEj335->SetMaterialPropertiesTable(ej335_mt);
  fEj335->GetIonisation()->SetBirksConstant(0.016*cm/MeV);
  //scintillation:
  /*light yeild = 55% anthracene
    anthracene = 17400 photons/MeV
    yield = 9.57/keV
    lifetime = 3.8 ns (singlet only)
  toFind = Rayleigh, Resolution Scale = 0.11, Birks Constant/ionization = 0.016cm/MeV*/
  //G4cout << "first definition of lAr below" << G4endl;

  //Define the Materials from liquid Argon to tpb.
  //Liquid Argon; three types defined for gradients of liquid argon absorption length
  lAr = new G4Material("lAr",z=18.,a=39.95*g/mole,density=1.396*g/cm3,kStateLiquid,88*kelvin);
  lAr2 = new G4Material("lAr2",z=18.,a=39.95*g/mole,density=1.396*g/cm3,kStateLiquid,88*kelvin);

  //G4cout << "first definition of lAr above" << G4endl;

  
  //Values for liquid Argon (all types)
  G4double lar_Energy_scin[] = { 3.87*eV , 4.51*eV , 4.74*eV , 5.03*eV , 5.36*eV , 5.55*eV , 5.82*eV , 6.06*eV , 6.54*eV , 6.79*eV , 7.03*eV , 7.36*eV , 7.76*eV , 7.98*eV , 8.33*eV , 8.71*eV , 8.96*eV , 9.33*eV , 9.91*eV , 10.31*eV , 10.61*eV , 10.88*eV , 11.27*eV , 11.81*eV , 12.40*eV , 13.05*eV , 13.78*eV , 14.59*eV , 15.50*eV }; //energies for scintillation spectrum
  
  G4double lar_Energy_rin[]    = { 1.771210*eV , 2.066412*eV , 2.479694*eV , 3.099618*eV ,
				   4.132823*eV , 6.199235*eV , 6.888039*eV , 7.749044*eV , 
				   8.856050*eV , 9.252590*eV , 9.686305*eV , 9.998766*eV , 
				   10.33206*eV , 11.27134*eV , 12.39847*eV }; //energies for refractive index and Rayleigh scattering lengths

  G4double lar_Energy_rs[]    = { 1.771210*eV , 2.066412*eV , 2.479694*eV , 3.099618*eV ,
				  4.132823*eV , 6.199235*eV , 6.888039*eV , 7.749044*eV , 
				  8.856050*eV , 9.252590*eV , 9.686305*eV , 9.998766*eV , 
				  10.33206*eV , 11.27134*eV , 12.39847*eV }; //energies for refractive index and Rayleigh scattering lengths*/

  const G4int larscin = sizeof(lar_Energy_scin)/sizeof(G4double);
  const G4int larrin =  sizeof(lar_Energy_rin)/sizeof(G4double);
  const G4int larrs  =  sizeof(lar_Energy_rs)/sizeof(G4double);
  //defines a few integers to make some tests easier.

  G4double lar_SCINT[] = { 0.00006, 0.00007, 0.00008, 0.00011, 0.00020, 0.00030, 0.00048, 0.00082, 0.00126, 0.00084, 0.00043, 0.00030, 0.00106, 0.00298, 0.00175, 0.00351, 0.01493, 0.12485, 0.49332, 0.20644, 0.07477, 0.04496, 0.01804, 0.00576, 0.00184, 0.00059, 0.00019, 0.00006, 0.00002 }; //liquid Argon scintillation spectrum. this one is centered at 128 nm (as it should be). 
  assert(sizeof(lar_SCINT) == sizeof(lar_Energy_scin)); //makes sure the energy and spectrum of the scintillation values are the same size.
  G4double lar_RIND[]  = { 1.22 , 1.222 , 1.225 , 1.23 ,
			   1.24 , 1.255 , 1.263 , 1.28 , 
			   1.315, 1.335 , 1.358 , 1.403, 
			   1.45 , 1.62  , 1.79 }; //index of refraction spectrum.
  G4double lar_RSL[]  = { 327028.6808*cm, 172560.2267*cm, 80456.5339*cm, 31177.44642*cm, 
			  8854.144327*cm, 1496.876298*cm, 906.5011168*cm, 480.2538294*cm, 
			  205.3758714*cm, 145.6326111*cm, 100.7813004*cm, 63.2898117*cm, 
			  40.07450411*cm, 11.43903548*cm, 3.626432195*cm }; //spectrum of rayleigh scattering lengths.*/ 
  assert(sizeof(lar_RIND) == sizeof(lar_Energy_rin)); 
  assert(sizeof(lar_RSL) == sizeof(lar_Energy_rs));

  //Takes the defined values above and uses them to define a materials properties table.
  lAr_mt = new G4MaterialPropertiesTable();//initiates the table
  lAr_mt->AddProperty("FASTCOMPONENT", lar_Energy_scin, lar_SCINT, larscin);//adds a property; this is the fast component of the scinitillation. 
  lAr_mt->AddProperty("SLOWCOMPONENT", lar_Energy_scin, lar_SCINT, larscin);
  lAr_mt->AddProperty("RINDEX",        lar_Energy_rin,  lar_RIND,  larrin);
  lAr_mt->AddProperty("RAYLEIGH",      lar_Energy_rs,   lar_RSL,   larrs);
  G4double scint_yeild=1.0/(19.5*eV); //scintillation yeild: 50 per keV.
  lAr_mt->AddConstProperty("SCINTILLATIONYIELD",scint_yeild);
  lAr_mt->AddConstProperty("RESOLUTIONSCALE",0.11);//numbers from https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf, slide on scintllation implementation.
  lAr_mt->AddConstProperty("FASTTIMECONSTANT",7.*ns);
  lAr_mt->AddConstProperty("SLOWTIMECONSTANT",1400.*ns);
  lAr_mt->AddConstProperty("YIELDRATIO",0.25);//for e/m scintillation
  //note: if you switch between nucleonic and e/m, make sure to change in the main file as well.
  lAr->SetMaterialPropertiesTable(lAr_mt);

  //Defines the MPT for the third type of liquid Argon. (only the absorption lenght is different). This kind is for the very cloudy bottom.
  lAr2_mt = new G4MaterialPropertiesTable();
  lAr2_mt->AddProperty("FASTCOMPONENT", lar_Energy_scin, lar_SCINT, larscin);
  lAr2_mt->AddProperty("SLOWCOMPONENT", lar_Energy_scin, lar_SCINT, larscin);
  lAr2_mt->AddProperty("RINDEX",        lar_Energy_rin,  lar_RIND,  larrin);
  lAr2_mt->AddProperty("RAYLEIGH",      lar_Energy_rs,   lar_RSL,   larrs);
  lAr2_mt->AddConstProperty("SCINTILLATIONYIELD",scint_yeild);
  lAr2_mt->AddConstProperty("RESOLUTIONSCALE",0.11);//numbers from https://indico.cern.ch/event/44566/contributions/1101918/attachments/943057/1337650/dipompeo.pdf, slide on scintllation implementation.
  lAr2_mt->AddConstProperty("FASTTIMECONSTANT",7.*ns);
  lAr2_mt->AddConstProperty("SLOWTIMECONSTANT",1400.*ns);
  lAr2_mt->AddConstProperty("YIELDRATIO",0.25);//for e/m scintillation
  lAr2->SetMaterialPropertiesTable(lAr2_mt);

  // Set the Birks Constant for the lAr scintillator (from the same paper as above)
  lAr->GetIonisation()->SetBirksConstant(0.0486*mm/MeV);
  lAr2->GetIonisation()->SetBirksConstant(0.0486*mm/MeV);
 
  //the following section defines the absorption lengths for the various kinds of liquid Argon.
  G4double base = 55.9506;//42.72;//Base absorption length for UV light
  G4double mult = 2800.0;//Absorption length for visible.
  DefineLAr(base,100.0,threehun,mult,1.358,1,100.7813004,1);//ultra is second


  //G4cout << "EdwardNote Defined Compounds" << G4endl;
  //***Material properties tables
  //these tables define the optical properties of the materials 
  //defines refractive index, scintillation properties, absorption length with respect to energy.
  //Definition of MPT for Glass
  G4double glass_Energy[] = { 7.0*eV, 7.07*eV, 7.14*eV };
  G4int glassnum = sizeof(glass_Energy) / sizeof(G4double);
  G4double glass_RIND[]={1.49,1.49,1.49};
  assert(sizeof(glass_RIND) == sizeof(glass_Energy));
  G4double glass_AbsLength[]={420.*cm,420.*cm,420.*cm};
  assert(sizeof(glass_AbsLength) == sizeof(glass_Energy));
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH",glass_Energy,glass_AbsLength,glassnum);
  glass_mt->AddProperty("RINDEX",glass_Energy,glass_RIND,glassnum);
  fGlass->SetMaterialPropertiesTable(glass_mt);


  //Definition of MPT for Plastic frills
  G4double plastic_Energy[] = { 1.0*eV,1.2*eV,2.5*eV,3.0*eV,3.4*eV,6.5*eV,10.0*eV,12.0*eV };
  //wavelengths: 1200, 1000, 500, 400, 350, 200, 120, 100
  G4int plasticnum = sizeof(plastic_Energy) / sizeof(G4double);
  G4double plastic_reflect[]={0.10, 0.10, 0.25, 0.30, 0.10, 0.05, 0.01, 0.01};
  assert(sizeof(plastic_RIND) == sizeof(plastic_Energy));
  G4double plastic_AbsLength[]={ 10.0*cm, 10.0*cm, 10.0*cm, 10.0*cm, 1.0e-3*cm, 1.0e-3*cm, 1.0e-3*cm, 1.0e-3*cm};
  G4double blackplastic_AbsLength[]={ 1e-9*cm, 1e-9*cm, 1e-9*cm, 1e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm};
  G4double plastic_RIND[] = { 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4 };
  assert(sizeof(plastic_AbsLength) == sizeof(plastic_Energy));
  G4MaterialPropertiesTable *plastic_mt = new G4MaterialPropertiesTable();
  plastic_mt->AddProperty("ABSLENGTH",plastic_Energy,plastic_AbsLength,plasticnum);
  plastic_mt->AddProperty("REFLECTIVITY",plastic_Energy,plastic_reflect,plasticnum);
  plastic_mt->AddProperty("RINDEX",plastic_Energy,plastic_RIND,plasticnum);
  fPlastic->SetMaterialPropertiesTable(plastic_mt);
  G4MaterialPropertiesTable *blackplastic_mt = new G4MaterialPropertiesTable();
  blackplastic_mt->AddProperty("ABSLENGTH",plastic_Energy,blackplastic_AbsLength,plasticnum);
  blackplastic_mt->AddProperty("REFLECTIVITY",plastic_Energy,plastic_reflect,plasticnum);
  blackplastic_mt->AddProperty("RINDEX",plastic_Energy,plastic_RIND,plasticnum);
  fBlackPlastic->SetMaterialPropertiesTable(blackplastic_mt);

  //Vacuum
  G4double vacuum_Energy[]={2.0*eV,7.0*eV,7.14*eV};
  const G4int vacnum = sizeof(vacuum_Energy)/sizeof(G4double);
  G4double vacuum_RIND[]={1.,1.,1.};
  assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND,vacnum);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  
  //Aluminum
  G4double alum_energy[40] = { 1.53353916299*eV, 1.55836885306*eV, 1.58716054918*eV, 1.61709864659*eV, 1.65046713174*eV,    
			       1.6817136839*eV,  1.71786713359*eV, 1.75433728722*eV, 1.79106381297*eV, 1.83220935967*eV, 
			       1.87670110428*eV, 1.92033477064*eV, 1.96440460941*eV, 2.01566223496*eV, 2.06248774244*eV, 
			       2.1191459895*eV,  2.17704608625*eV,  2.23198006977*eV, 2.29412813935*eV,2.37159279237*eV, 
			       2.43945693286*eV, 2.51927088373*eV, 2.59881298342*eV, 2.68963236717*eV, 2.78391786144*eV, 
			       2.88535232931*eV, 2.97955997826*eV, 3.08041053689*eV, 3.18832742402*eV, 3.26833943848*eV, 
			       3.36191419687*eV, 3.46107669506*eV, 3.52980225748*eV, 3.60084752681*eV, 3.70923697374*eV, 
			       3.75238245064*eV, 3.82680283055*eV, 3.9042347515*eV,  3.9111734511*eV, 20.6640330667*eV };
  G4double alum_reflect[40] = {0.8866682, 0.8975349, 0.9091666, 0.9176283, 0.9260717, 
			       0.9329575, 0.9398158, 0.9458907, 0.9511823, 0.9548705, 
			       0.9609272, 0.9662096, 0.9722937, 0.9775578, 0.9836511, 
			       0.9865376, 0.9886408, 0.9891864, 0.9905061, 0.9909877, 
			       0.9907316, 0.9904481, 0.9901829, 0.9891068, 0.9856623, 
			       0.9774627, 0.968507,  0.955589,  0.9426708, 0.9242781, 
			       0.9058671, 0.8866637, 0.862769,  0.8436296, 0.8228502, 
			       0.8005955, 0.7798801, 0.7591648, 0.001,     0.001 };
  assert(sizeof(alum_energy) == sizeof(alum_reflect));
  //Above: reflectivity of aluminum for frame. Below: absorption length for aluminum (very short).
  G4double alum_abslen[2] = { 1.0e-3*cm, 1.0e-3*cm};
  G4double alum_abseneg[2] = { 20.664*eV, 1.5498*eV };
  assert(sizeof(alum_abseneg) == sizeof(alum_abslen));
  G4MaterialPropertiesTable *alum_mt = new G4MaterialPropertiesTable();
  alum_mt->AddProperty("REFLECTIVITY", alum_energy, alum_reflect,40);
  alum_mt->AddProperty("ABSLENGTH", alum_abseneg, alum_abslen,2);
  alum->SetMaterialPropertiesTable(alum_mt);
  
  //Definition of Material properties for TPB. 
  const G4int nTPBEntries = 25;
  // Comment: These are different energies than pretty much all the other arrays (at least from 2.583->2.980eV)...                        
  G4double TPBEnergy[nTPBEntries] =
    { 0.602*eV/*(2066nm)*/, 0.689*eV/*(1799nm)*/, 1.030*eV/*(1204nm)*/, 1.926*eV/*(644nm)*/, 2.138*eV/* (580nm)*/, 
      2.250*eV/*(551nm)*/,  2.380*eV/*(521nm)*/,  2.480*eV/*(500nm)*/,  2.583*eV/*(480nm)*/, 2.800*eV/*(443nm)*/,
      2.880*eV/*(431nm)*/,  2.980*eV/*(416nm)*/,  3.124*eV/*(397nm)*/,  3.457*eV/*(359nm)*/, 3.643*eV/*(341nm)*/,
      3.812*eV/*(325nm)*/,  4.086*eV/*(304nm)*/,  4.511*eV/*(275nm)*/,  5.166*eV/*(240nm)*/, 5.821*eV/*(213nm)*/,
      6.526*eV/*(190nm)*/,  8.266*eV/*(150nm)*/,  9.686*eV/*(128nm)*/,  11.27*eV/*(110nm)*/, 12.60*eV/*(98nm)*/  };
                                                      
  //Emission spectrum for TPB; basically probabilities of photon being reemitted in the various energy bins.
  G4double TPBEmission[nTPBEntries] = //TPB Emission spectrum 
    {
      0.0000, 0.0000, 0.0000, 0.0000, 0.0005,
      0.0015, 0.0030, 0.0050, 0.0070, 0.0110,
      0.0110, 0.0060, 0.0020, 0.0000, 0.0000,
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
      0.0000, 0.0000, 0.0000, 0.0000, 0.0000
    };

  G4double TPBRIndex[nTPBEntries] =  //Refractive index of the TPB.
    {
      1.4, 1.4, 1.4, 1.4, 1.4,
      1.4, 1.4, 1.4, 1.4, 1.4,
      1.4, 1.4, 1.4, 1.4, 1.4,
      1.4, 1.4, 1.4, 1.4, 1.4,
      1.4, 1.4, 1.4, 1.4, 1.4
    };//Claim made for 1.7 from DEAP Thesis. Maybe try varying to check?. TPB rayleigh scattering length is also listed at 3 um, maybe try throwing that component in instead of randomization in stepping action. 

  //defining the MPTs for the various types of TPB. This version for foils
  TPBProp = new G4MaterialPropertiesTable();
  //TPBProp->AddProperty("RINDEX", TPBEnergy, TPBRIndex, nTPBEntries);
  TPBProp->AddProperty("WLSCOMPONENT", TPBEnergy, TPBEmission, nTPBEntries);
  TPBProp->AddConstProperty("WLSTIMECONSTANT", 1.7*ns);
  TPBProp->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.2);
  //allows geant to produce multiple photons per incoming, with the mean number being 1.2. if off WLSabsorbed to emitted is 1:1
  tPB->SetMaterialPropertiesTable(TPBProp);
  
  //MPT for second kind of TPB, this version for PMTs.
  TPBsProp = new G4MaterialPropertiesTable();
  //TPBsProp->AddProperty("RINDEX", TPBEnergy, TPBRIndex, nTPBEntries);
  TPBsProp->AddProperty("WLSCOMPONENT", TPBEnergy, TPBEmission, nTPBEntries);
  TPBsProp->AddConstProperty("WLSTIMECONSTANT", 1.7*ns);
  TPBsProp->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.2);
  tPBhundred->SetMaterialPropertiesTable(TPBsProp);

  DefineTpb(foilEff,tpbEff,tpbAbs,1.7,2.75,2.75);

  //Defines properties of the ptfe reflectors.
  const G4int nTefEntries = 25;
  G4double TefEnergy[nTefEntries] =
    {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV,
     2.138*eV, 2.25*eV,  2.38*eV,  2.48*eV, 
     2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV, 
     3.124*eV, 3.457*eV, 3.643*eV, 3.812*eV, 4.086*eV,
     4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
     7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
  G4double TefReflect[nTefEntries] =
    {1., 1., 1., 1., 
     1., 1., 1., 1., 
     1., 1., 1., 1.,
     1., 1., 1., 1., 1.,
     1., 1., 1., 1., 
     1., 1., 1., 0.99};
  G4double TefRIndex[nTefEntries] =
    {10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
     10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
     10., 10., 10., 10., 10.};
  G4MaterialPropertiesTable *TefProp = new G4MaterialPropertiesTable();
  TefProp->AddProperty("RINDEX", TefEnergy, TefRIndex, nTefEntries);
  TefProp->AddProperty("REFLECTIVITY", TefEnergy, TefReflect, nTefEntries);
  ptfe->SetMaterialPropertiesTable(TefProp);

  //G4cout << "EdwardNote Defined Material properties" << G4endl;
}

//method for more consistently performing the LAr Definition
void detectorConstruction::DefineLAr(G4double base, G4double uvlas, G4double three, G4double mult, G4double rindex, G4double scaler, G4double ray128, G4double scalerray) {
  //the following section defines the absorption lengths for the various kinds of liquid Argon.
  G4double time5 = 1.0;//five/100.;
  G4double row5 = time5*base;
  G4double uvlas5 = time5*uvlas;
  G4double three5 = three;
  G4double mult5 = mult;

  G4double lar_Energy_abs[47];
  G4double lar_wlv_abs[47];
  for (int i=0;i<47;++i){
    lar_wlv_abs[i] = i*20+80;
    lar_Energy_abs[i] = (1239.847/(i*20.0+80.0))*eV;
  }
  G4double lar_ABSL[47];
  G4double lar2_ABSL[47];
  for (int i=0;i<6;++i){
    lar_ABSL[i] = base*cm;
    lar2_ABSL[i] = row5*cm;
    //G4cout << lar_wlv_abs[i] << '\t' << lar_ABSL[i] << '\t' << lar2_ABSL[i] << G4endl; 
  }
  for (int i=6;i<11;++i){
    lar_ABSL[i] = uvlas*cm;
    lar2_ABSL[i] = uvlas5*cm;
    //G4cout << lar_wlv_abs[i] << '\t' << lar_ABSL[i] << '\t' << lar2_ABSL[i] << G4endl; 
  }
  for (int i=11;i<16;++i){
    lar_ABSL[i] = three*cm;
    lar2_ABSL[i] = three5*cm;
    //G4cout << lar_wlv_abs[i] << '\t' << lar_ABSL[i] << '\t' << lar2_ABSL[i] << G4endl; 
  }
  for (int i=16;i<47;++i){
    lar_ABSL[i] = mult*cm;
    lar2_ABSL[i] = mult5*cm;
    //G4cout << lar_wlv_abs[i] << '\t' << lar_ABSL[i] << '\t' << lar2_ABSL[i] << G4endl; 
  }
  const G4int larabs =  sizeof(lar_Energy_abs)/sizeof(G4double);
  assert(sizeof(lar_ABSL) == sizeof(lar_Energy_abs));

  //Wavelength list = 700 600 500 400 300 200 180 160 140 134 128 124 120 110 100
  G4double lar_Energy_rin[]    = { 1.771210*eV , 2.066412*eV , 2.479694*eV , 3.099618*eV ,
				   4.132823*eV , 6.199235*eV , 6.888039*eV , 7.749044*eV , 
				   8.856050*eV , 9.252590*eV , 9.686305*eV , 9.998766*eV , 
				   10.33206*eV , 11.27134*eV , 12.39847*eV }; //energies for refractive index and Rayleigh scattering lengths

  
  G4double set = rindex-1.24;
  //G4double scaler = 1.0;
  G4double scale[] = { 0.125, 0.196, 0.345, 0.637, 0.8065, 1.38, 1.78, 3.22, 4.66 };
  for (int i = 0; i < 9; ++i) {
    scale[i] = std::pow(scale[i],scaler)*set+1.24;
  }//*/
  G4double rayleigh = ray128*cm;//100.7813004*cm;
  G4double lar_RSL[] = { 3244.9341, 1712.2246, 798.328,   309.35745,
			 87.855032, 14.852719, 8.9947353, 4.7653069, 
			 2.0378371, 1.445036,  1,         0.6279916, 
			 0.3976383, 0.1135036, 0.035983 };
  for (int i = 0; i < 15; ++i) {
    /*if (lar_RSL[i] == 1) { lar_RSL[i] = lar_RSL[i]*rayleigh; }
    else if (lar_RSL[i] < 1) { lar_RSL[i] = lar_RSL[i]/scalerray*rayleigh; }
    else if (lar_RSL[i] > 1) { lar_RSL[i] = lar_RSL[i]*scalerray*rayleigh; }*/
    lar_RSL[i] = std::pow(lar_RSL[i],scalerray)*rayleigh;
  }//*/
    
    
  const G4int larrin =  sizeof(lar_Energy_rin)/sizeof(G4double);
  G4double lar_RIND[]  = { 1.22 ,     1.222 ,    1.225 ,    1.23 ,
			   1.24 ,     scale[0] , scale[1] , scale[2] , 
			   scale[3] , scale[4] , rindex ,   scale[5], 
			   scale[6] , scale[7] , scale[8] }; //index of refraction spectrum.
  /*G4double lar_RSL[]  = { 327028.6808*cm, 172560.2267*cm, 80456.5339*cm, 31177.44642*cm, 
			  8854.144327*cm, 1496.876298*cm, 906.5011168*cm, 480.2538294*cm, 
			  205.3758714*cm, 145.6326111*cm, 100.7813004*cm, 63.2898117*cm, 
			  40.07450411*cm, 11.43903548*cm, 3.626432195*cm }; //spectrum of rayleigh scattering lengths.*/
  assert(sizeof(lar_RIND) == sizeof(lar_Energy_rin)); 
  assert(sizeof(lar_RSL) == sizeof(lar_Energy_rin));

  lAr_mt->AddProperty("RINDEX",        lar_Energy_rin,  lar_RIND,  larrin);
  lAr_mt->AddProperty("RAYLEIGH",      lar_Energy_rin,  lar_RSL,   larrin);//*/
  
  //Takes the defined values above and uses them to define a materials properties table.
  lAr_mt->AddProperty("ABSLENGTH",     lar_Energy_abs,  lar_ABSL,  larabs);
  lAr2_mt->AddProperty("ABSLENGTH",     lar_Energy_abs,  lar2_ABSL,  larabs);
  //G4cout << "second definition of lAr above" << G4endl;
}


//method for performing the Plastic material Definition
void detectorConstruction::DefinePlastic(G4double abs, G4double refl, G4double rin) {
  //Defaults: absl = 10.0*cm, refl = 0.10; max refl = 0.30
  //Definition of MPT for Plastic frills
  G4double plastic_Energy[] = { 1.0*eV,1.2*eV,2.5*eV,3.0*eV,3.4*eV,6.5*eV,10.0*eV,12.0*eV };
  //wavelengths: 1200, 1000, 500, 400, 350, 200, 120, 100
  G4int plasticnum = sizeof(plastic_Energy) / sizeof(G4double);

  G4double uvr = refl/10.0;
  G4double refl2 = refl*2.5;
  G4double refl3 = refl*3.0;
  G4double plastic_reflect[]={refl, refl, refl2, refl3, refl, uvr, uvr, uvr};
  assert(sizeof(plastic_RIND) == sizeof(plastic_Energy));

  G4double absl = abs*cm;
  G4double plastic_AbsLength[]={ absl, absl, absl, absl, 1.0e-3*cm, 1.0e-3*cm, 1.0e-3*cm, 1.0e-3*cm};
  G4double blackplastic_AbsLength[]={ 1e-9*cm, 1e-9*cm, 1e-9*cm, 1e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm};

  G4double plastic_RIND[] = { rin, rin, rin, rin, rin, rin, rin, rin };
  assert(sizeof(plastic_AbsLength) == sizeof(plastic_Energy));

  G4MaterialPropertiesTable *plastic_mt = new G4MaterialPropertiesTable();
  plastic_mt->AddProperty("ABSLENGTH",plastic_Energy,plastic_AbsLength,plasticnum);
  plastic_mt->AddProperty("REFLECTIVITY",plastic_Energy,plastic_reflect,plasticnum);
  plastic_mt->AddProperty("RINDEX",plastic_Energy,plastic_RIND,plasticnum);
  fPlastic->SetMaterialPropertiesTable(plastic_mt);
  G4MaterialPropertiesTable *blackplastic_mt = new G4MaterialPropertiesTable();
  blackplastic_mt->AddProperty("ABSLENGTH",plastic_Energy,blackplastic_AbsLength,plasticnum);
  blackplastic_mt->AddProperty("REFLECTIVITY",plastic_Energy,plastic_reflect,plasticnum);
  blackplastic_mt->AddProperty("RINDEX",plastic_Energy,plastic_RIND,plasticnum);
  fBlackPlastic->SetMaterialPropertiesTable(blackplastic_mt);
}//*/

//method to more consistently perform the TPB definitions (pmt and foil)
void detectorConstruction::DefineTpb(G4double foil, G4double pmt, G4double abs, G4double rin, G4double ral, G4double ray) {
  const G4int nTPBEntries = 25;
  //Redefining TPB efficiency values for both foil and PMTs
  G4double TPBEnergy[nTPBEntries] =
    { 0.602*eV/*(2066nm)*/, 0.689*eV/*(1799nm)*/, 1.030*eV/*(1204nm)*/, 1.926*eV/*(644nm)*/, 2.138*eV/* (580nm)*/, 
      2.250*eV/*(551nm)*/,  2.380*eV/*(521nm)*/,  2.480*eV/*(500nm)*/,  2.583*eV/*(480nm)*/, 2.800*eV/*(443nm)*/,
      2.880*eV/*(431nm)*/,  2.980*eV/*(416nm)*/,  3.124*eV/*(397nm)*/,  3.457*eV/*(359nm)*/, 3.643*eV/*(341nm)*/,
      3.812*eV/*(325nm)*/,  4.086*eV/*(304nm)*/,  4.511*eV/*(275nm)*/,  5.166*eV/*(240nm)*/, 5.821*eV/*(213nm)*/,
      6.526*eV/*(190nm)*/,  8.266*eV/*(150nm)*/,  9.686*eV/*(128nm)*/,  11.27*eV/*(110nm)*/, 12.60*eV/*(98nm)*/  };

  G4double lengthconst = 0.0019/(std::log(1-foil));
  G4double wlAbf11 = lengthconst/(std::log(0.067));
  G4double wlAbf12 = lengthconst/(std::log(0.2));
  G4double wlAbf15 = lengthconst/(std::log(0.4));
  G4double wlAbf19 = lengthconst/(std::log(0.533));
  G4double wlAbf21 = lengthconst/(std::log(0.4));
  G4double wlAbf24 = lengthconst/(std::log(0.367));
  
  if (!ccm200) {
    wlAbf21 = -0.002/(std::log((1-foil)));
    wlAbf24 = -0.002/(std::log((1-(.222*foil/0.2))));
    wlAbf19 = -0.002/(std::log((1-(.156*foil/0.2))));
    wlAbf15 = -0.002/(std::log((1-(.114*foil/0.2))));
    wlAbf12 = -0.002/(std::log((1-(.191*foil/0.2))));
    wlAbf11 = wlAbf21/(std::log(5.797944));
  }    
    
  G4double TPBWLSAbsorption[nTPBEntries] =
    { 0.10000*m, 1000.000*m, 1000.000*m, 1000.000*m, 1000.000*m,
      1000.000*m, 1000.000*m, 1000.000*m, 1000.000*m, 1000.000*m,
      10000.000*m, 10000.000*m, 10000.000*m, 10000.000*m, 100000.0*m,
      100000.0*m, 100000.0*m, 100000.0*m, wlAbf24*mm, wlAbf21*mm, 
      wlAbf19*mm, wlAbf15*mm, wlAbf12*mm, wlAbf11*mm, wlAbf11*mm, };
    
  
  lengthconst = 0.0019/(std::log(1-pmt));
  G4double wlsAb110 = lengthconst/(std::log(0.067));
  G4double wlsAb128 = lengthconst/(std::log(0.2));
  G4double wlsAb190 = lengthconst/(std::log(0.533));
  G4double wlsAb213 = lengthconst/(std::log(0.4));
  G4double wlsAb240 = lengthconst/(std::log(0.367));
  G4double TPBWLSAbsorption100[nTPBEntries] = 
    {    0.10000*m, 1000.00*m, 1000.00*m, 1000.00*m, 1000.000*m,
	 1000.00*m, 1000.00*m, 1000.00*m, 1000.00*m, 1000.000*m,
	 10000.0*m, 10000.0*m, 10000.0*m, 10000.0*m, 100000.0*m, 
	 100000.0*m, 100000.0*m, 100000.0*m, wlsAb240*mm, wlsAb213*mm, 
	 wlsAb190*mm, wlsAb213*mm, wlsAb128*mm, wlsAb110*mm, wlsAb110*mm
    };

  std::ostringstream oss;

  oss << "TPBAbsorptionLengths: " << wlAbf11 << ',' << wlAbf12 << ',' << wlAbf15 << ',' << wlAbf19 << ',' << wlAbf21 << ',' << wlAbf24 << ',' << '\n' << wlsAb110 << ',' << wlsAb128 << ',' << wlsAb190 << ',' << wlsAb213 << ',' << wlsAb240 << ',' << G4endl;

  G4String tpbstring = oss.str();
  //G4cout << tpbstring << G4endl;

  G4double absltwo = -0.00211/(std::log(abs));
  G4double abslttw = absltwo*1.4;
  G4double abslttt = absltwo*1.55;
    
  G4double TPBAbsorption[nTPBEntries] = 
    {  0.02000*mm, absltwo*mm, absltwo*mm, absltwo*mm, absltwo*mm,
       abslttw*mm, abslttw*mm, abslttw*mm, abslttw*mm, abslttw*mm,
       abslttt*mm, abslttt*mm, abslttt*mm, 10.0000*mm, 100000.0*m,
       100000.0*m, 100000.0*m, 100000.0*m, 100000.0*m, 100000.0*m,
       100000.0*m, 100.0000*m, 100.0000*m, 100.0000*m, 100.0000*m
    };

  //G4double rin = 1.7;
  G4double rayl = ray*0.001*mm;
  G4double rayl2 = ral*0.001*mm;
  G4double TPBRIndex[nTPBEntries] =  //Refractive index of the TPB.
    {
      rin, rin, rin, rin, rin,
      rin, rin, rin, rin, rin,
      rin, rin, rin, rin, rin,
      rin, rin, rin, rin, rin,
      rin, rin, rin, rin, rin
    };//Claim made for 1.7 from DEAP Thesis. Maybe try varying to check?. TPB rayleigh scattering length is also listed at 2.75 um, maybe try throwing that component in instead of randomization in stepping action.
  
  G4double TPBRayleigh[nTPBEntries] =  //Rayleigh Scattering length for TPB
    {
      rayl, rayl, rayl, rayl, rayl,
      rayl, rayl, rayl, rayl, rayl,
      rayl, rayl, rayl, rayl, rayl,
      rayl, rayl, rayl, rayl, rayl,
      rayl, rayl, rayl, rayl, rayl
    };

  G4double TPBRayleigh2[nTPBEntries] =  //Rayleigh Scattering length for TPB
    {
      rayl2, rayl2, rayl2, rayl2, rayl2,
      rayl2, rayl2, rayl2, rayl2, rayl2,
      rayl2, rayl2, rayl2, rayl2, rayl2,
      rayl2, rayl2, rayl2, rayl2, rayl2,
      rayl2, rayl2, rayl2, rayl2, rayl2
    };

  //foil TPB properties
  //Split the rayleigh scattering lengths into 2 different variables: one for PMTs, one for Foils (foils look more scattered than PMTs). 
  TPBProp->AddProperty("RINDEX", TPBEnergy, TPBRIndex, nTPBEntries);
  TPBProp->AddProperty("RAYLEIGH", TPBEnergy, TPBRayleigh2, nTPBEntries);
  TPBProp->AddProperty("WLSABSLENGTH", TPBEnergy, TPBWLSAbsorption, nTPBEntries);
  TPBProp->AddProperty("ABSLENGTH", TPBEnergy, TPBAbsorption, nTPBEntries);

  //PMT TPB properties
  TPBsProp->AddProperty("RINDEX", TPBEnergy, TPBRIndex, nTPBEntries);
  TPBsProp->AddProperty("RAYLEIGH", TPBEnergy, TPBRayleigh, nTPBEntries);
  TPBsProp->AddProperty("WLSABSLENGTH", TPBEnergy, TPBWLSAbsorption100, nTPBEntries);
  TPBsProp->AddProperty("ABSLENGTH", TPBEnergy, TPBAbsorption, nTPBEntries);

}

//Construct method called in main code to actual build the detector using the materials defined above. 
//Also contains the if-else statements that build parts of the detector that are not always on.
G4VPhysicalVolume* detectorConstruction::Construct(){
  

  //Create a box 10x10x10 meters in which to place the detector. In theory is much larger than it needs to be.
  G4double expHall_x = 5*m;
  G4double expHall_y = 5*m;
  G4double expHall_z = 5*m;
  wBox = new G4Box("expHall",expHall_x,expHall_y,expHall_z);
  wLog = new G4LogicalVolume(wBox,fVacuum,"expHall",0,0,0);
  wPhys = new G4PVPlacement(0,G4ThreeVector(),wLog,"expHall",0,false,0);

  wLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  //create a slightly smaller box in order to kill things which go too far
  expHall_x = 4.9*m;
  expHall_y = 4.9*m;
  expHall_z = 4.9*m;
  G4Box* iBox = new G4Box("inner_box",expHall_x,expHall_y,expHall_z);
  G4LogicalVolume* iLog = new G4LogicalVolume(iBox,fVacuum,"inner_log");
  new G4PVPlacement(0,G4ThreeVector(0*cm, 0*cm, 0*cm),iLog,"inner",wLog,false,0,true);

  iLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  //makes the smaller box volume invisible to visualization software (as it should be)

  //Place the main volume if it is on. (i.e. the CCM detector)
  if(mainVolOn){
    //add concentric cylinders to form CCM detector volume: cryo (G4_STAINLESS-STEEL), LAr, support (Aluminum, G4_Al), inner support (ptfe, G4_TEFLON), TPB coating (TPB, get properties), fiducial Argon.

    //Outer cryogen
    fCryoVessel = new G4Tubs("Cryogen", 0*cm, 138*cm, 131*cm, 0*deg,360*deg);
    fLogicCryo  = new G4LogicalVolume(fCryoVessel, steel, "Cryogen");
    new G4PVPlacement(0,
		      G4ThreeVector(0*cm, 0*cm, 0*cm),
		      fLogicCryo,
		      "Cryogen",
		      iLog,
		      false,
		      0,
		      true);
    
    //Vacuum jacket
    G4Tubs* vacuum = new G4Tubs("Vacuum", 0*cm, 135*cm, 126*cm, 0*deg, 360*deg);
    G4LogicalVolume* fLogicVacuum = new G4LogicalVolume(vacuum, fVacuum, "Vacuum");
    new G4PVPlacement(0,
		      G4ThreeVector(0*cm, 0*cm, 0*cm),
		      fLogicVacuum,
		      "Vacuum",
		      fLogicCryo,
		      false,
		      0,
		      true);

    //inner cryogen
    G4Tubs* innerjack = new G4Tubs("innerjacket", 0*cm, 125*cm, 120*cm, 0*deg, 360*deg);
    G4LogicalVolume* fLogicJacket = new G4LogicalVolume(innerjack, steel, "innerjacket");
    new G4PVPlacement(0,
		      G4ThreeVector(0*cm, 0*cm, 0*cm),
		      fLogicJacket,
		      "innerjacket",
		      fLogicVacuum,
		      false,
		      0,
		      true);

    //Argon outside the fiducial volume
    fArgonOuter = new G4Tubs("liquidArgon", 0*cm, 120*cm, 115*cm, 0*deg, 360*deg);
    fLogicArout = new G4LogicalVolume(fArgonOuter, lAr, "liquidArgon");
    new G4PVPlacement(0,
		      G4ThreeVector(0*cm, 0*cm, 0*cm),
		      fLogicArout,
		      "liquidArgon",
		      fLogicJacket,
		      false,
		      0,
		      true);
    
    //Aluminum frame holding PMTs and instrumentation
    fInnerFrame = new G4Tubs("Frame",0*cm, 106*cm, 75*cm, 0*deg, 360*deg);
    fLogicFrame = new G4LogicalVolume(fInnerFrame,alum, "Frame");
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      fLogicFrame,
		      "Frame",
		      fLogicArout,
		      false,
		      0,
		      true);

    //Defines the thickness and placement of the TPB foils (slightly thicker on the bottom than the top). Can change the overall thickness without altering the ratios by just adjusting the first number, basethick
    const G4double basethick = 0.00019*cm;
    G4double thick = basethick*conewide;///fifth;//Define proportional foil thickness here: prel200. the best fit is half the normal thickness

    G4double deep = thick;//Defines the depth; used for making the bottom thicker
    G4double place = thick*0.0;//(1-randwide)/(1+randwide);//(thick-thick);//Defines the place; if the bottom is thicker the TPB cylinder needs to move slightly.

    //G4cout << "EdwardNote Defined Up to foils" << G4endl;
    
    //defines volumes for the TPB on the foils if it is on. Defines two seperate TPB volumes to distinguish the TPB on the top/bottom (unsmooth foil) from the TPB on the sides (smooth foil).
    G4double totalH = 136.0/2.0;
    if (ccm200) { totalH = 123.2/2.0; }

    if (fTPBfoilOn) {
      fTPBSides = new G4Tubs("TPBfoil", 0*cm, (96.0*cm+thick), totalH*cm, 0*deg, 360*deg);
      fLogicTPB = new G4LogicalVolume(fTPBSides,tPB,"TPBfoil");//*/
      
      fTPBBottom = new G4Tubs("TPBfoilb", 13.3*cm, (96.0*cm+thick), (totalH*cm+deep), 0*deg, 360*deg);
      
      G4SubtractionSolid* tpbBottom = new G4SubtractionSolid("TPBfoilsub", fTPBBottom, fTPBSides, 0, G4ThreeVector(0*cm, 0*cm, place));
      //G4cout << "TPBvolumahalfZlength = " << fTPBBottom->GetZHalfLength() << G4endl;
      fLogicTPBb = new G4LogicalVolume(tpbBottom,tPB,"TPBfoilb");//*/
    }
    
    //defines the optical surface of the TPB. Reflection, transmission, and efficiency of photons that intersect the surface.
    //currently set to 1's and 0's so all physics occurs in the volume itself.
    const G4int nAcTefEntries = 25;
    G4double TPBfoilOSEnergy[nAcTefEntries] =
      {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV, 2.38*eV,
       2.48*eV, 2.583*eV, 2.845*eV,  2.857*eV, 2.95*eV, 3.124*eV,  3.457*eV, 
       3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
       7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
    G4double TPBfoilOSTransmit[nAcTefEntries] =
      {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1., 1., 1., 1.}; //set to 1 and have all absorption in bulk 
    G4double TPBfoilOSReflect[nAcTefEntries] =
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0., 0., 0., 0., 0.}; //set to zero for the same
    G4double TPBfoilOSEff[nAcTefEntries] =
      {0., 0., 0., 0., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
       1.0, 1.0, 1.0, 1.0, 1., 1., 1., 1., 1., 1.,
       1., 1., 1., 1., 1.};

    G4OpticalSurface *TPBfoilOS = new G4OpticalSurface("TPB");
    TPBfoilOS->SetModel(unified); //Optical model    
    TPBfoilOS->SetType(dielectric_dielectric);
    TPBfoilOS->SetFinish(ground);
    TPBfoilOS->SetSigmaAlpha(0.05);
    
    G4MaterialPropertiesTable *TPBfoilMPT = new G4MaterialPropertiesTable();
    TPBfoilMPT->AddProperty("REFLECTIVITY", TPBfoilOSEnergy, TPBfoilOSReflect, nAcTefEntries);
    TPBfoilMPT->AddProperty("TRANSMITTANCE", TPBfoilOSEnergy, TPBfoilOSTransmit, nAcTefEntries);
    TPBfoilMPT->AddProperty("EFFICIENCY", TPBfoilOSEnergy, TPBfoilOSEff, nAcTefEntries);
    TPBfoilOS->SetMaterialPropertiesTable(TPBfoilMPT);

    //define the fiducal volumes of liquid Argon.
    fFiducialAr = new G4Tubs("Fiducial", 0*cm, 96*cm, totalH*cm, 0*deg, 360*deg);
    if (fInvBeta) {
      fLogicFiduc = new G4LogicalVolume(fFiducialAr,fEj335,"Fiducial");
    } else {
      fLogicFiduc = new G4LogicalVolume(fFiducialAr,lAr,"Fiducial");
    }
    //G4cout << "FiducialhalfZlength = " << fFiducialAr->GetZHalfLength() << G4endl;
    
    //defines layers for the lAr volumes, rather than having constant properties across all depths. 
    //These layers can be made as solid cylinders (0*cm to 96*cm) or as shells to modify the properties only near the outer edges (current: 70 to 96*cm)
    //each layer can have a different type of lAr (currently only two are available, but more can be added). 
    G4double innerradius = 96*cm - r5radius*cm;
    fFiducialAr2 = new G4Tubs("Fiducial2", 70*cm, 96*cm, 12.0*cm, 0*deg, 360*deg);
    fLogicFiduc2 = new G4LogicalVolume(fFiducialAr2,lAr,"Fiducial2");

    fFiducialAr3 = new G4Tubs("Fiducial3", 70*cm, 96*cm, 12.0*cm, 0*deg, 360*deg);
    fLogicFiduc3 = new G4LogicalVolume(fFiducialAr3,lAr,"Fiducial3");

    fFiducialAr4 = new G4Tubs("Fiducial4", 70*cm, 96*cm, 12.0*cm, 0*deg, 360*deg);
    fLogicFiduc4 = new G4LogicalVolume(fFiducialAr4,lAr,"Fiducial4");//*/

    G4double remainder = totalH - 36.0;
    remainder = remainder/2.0;
    G4double row1place = totalH-remainder;//134 -> 52
    //fFiducialAr1 = new G4Tubs("Fiducial1", innerradius, 96*cm, remainder*cm, 0*deg, 360*deg);
    //fLogicFiduc1 = new G4LogicalVolume(fFiducialAr1,lAr2,"Fiducial1");

    //The bottom layer is different from the prior (best fit so far). It has much cloudier lAR and a much larger contaminated ring.
    fFiducialAr5 = new G4Tubs("Fiducial5", 75*cm, 96*cm, remainder*cm, 0*deg, 360*deg);
    fLogicFiduc5 = new G4LogicalVolume(fFiducialAr5,lAr,"Fiducial5");

    //G4cout << "EdwardNote Defined TPB OS table" << G4endl;

    //defines the optical surface for the reflective foils. Does not use all the same values as the TPB in order to properly get reflectivity.
    G4OpticalSurface *reflfoilOS = new G4OpticalSurface("refl");
    
    reflfoilOS->SetModel(glisur); //Optical model                    
    reflfoilOS->SetType(dielectric_metal);
    reflfoilOS->SetFinish(ground);
    reflfoilOS->SetSigmaAlpha(1.0);
    
    G4MaterialPropertiesTable *reflfoilMPT = new G4MaterialPropertiesTable();
    G4double reflfoilOSEnergy[nAcTefEntries] =
      {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV,  2.38*eV,
       2.48*eV,  2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,  3.124*eV, 3.457*eV,
       3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
       7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
    G4double uvRf = 0.10; //change uvreflection
    G4double vsRf = 0.95; //change visible reflectivity
    //enter the uv and visible reflectivities defined above.
    G4double reflfoilOSReflect[nAcTefEntries] =
      {vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, 
       vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, 
       vsRf, vsRf, vsRf, vsRf, uvRf, uvRf, uvRf, 
       uvRf, uvRf, uvRf, uvRf};
    G4double reflfoilOSTransmit[nAcTefEntries] =
      {.02, .02, .02, .02, .02, .02, .02, .02, .02, .02,
       .02, .02, .02, .02, .02, .02, .02, .02, .02, .02,
       .02, .02, .02, .02, .02}; //minimal transmittance through foil
    G4double reflfoilOSEff[nAcTefEntries] =
      {0., 0., 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0, 0.0, 0.0, 0.0, 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0.};
    reflfoilMPT->AddProperty("TRANSMITTANCE", reflfoilOSEnergy, reflfoilOSTransmit, nAcTefEntries);
    reflfoilMPT->AddProperty("REFLECTIVITY", reflfoilOSEnergy, reflfoilOSReflect, nAcTefEntries);
    reflfoilMPT->AddProperty("EFFICIENCY", reflfoilOSEnergy, reflfoilOSEff, nAcTefEntries);
    reflfoilOS->SetMaterialPropertiesTable(reflfoilMPT);

    //if the reflector is on, place it.
    if (fReflectorOn) {
      //Places reflector if reflector is on
      G4Tubs* ffoilSides = new G4Tubs("ptfefoil", 0*cm, 96.05*cm, (totalH+0.05)*cm, 0*deg, 360*deg);
      G4LogicalVolume* fLogicfoil = new G4LogicalVolume(ffoilSides,ptfe,"ptfefoil");
      G4VPhysicalVolume* fPhysfoil = new G4PVPlacement(0,
						       G4ThreeVector(0,0,0),
						       fLogicfoil,
						       "ptfefoil",
						       fLogicFrame,
						       false,
						       0,
						       true);

      //if TPB foil is on place it and the fiducial volume within. 
      if (fTPBfoilOn){
	//TPB and Fiducial for both TPB foils and Reflector On
	
	fPhysTPBb = new G4PVPlacement(0,
				      G4ThreeVector(0,0,0),
				      fLogicTPBb,
				      "TPBfoilb",
				      fLogicfoil,
				      false,
				      0,
				      true);//*/
	

	fPhysTPB = new G4PVPlacement(0,
				     G4ThreeVector(0*cm,0*cm,place),
				     fLogicTPB,
				     "TPBfoil",
				     fLogicfoil,
				     false,
				     0,
				     true);//*/
	
	/*G4double sidedown = 10.0;
        G4double sideplace = (totalH-sidedown)*cm;
        G4double inrad = 96.0*cm+thick/4.0;
        G4Tubs* ffoilDown = new G4Tubs("ptfedown", inrad, (96.0*cm+thick), sidedown*cm, 0*deg, 360*deg);
        G4LogicalVolume* fLogicDown = new G4LogicalVolume(ffoilDown,ptfe,"ptfedown");
        G4VPhysicalVolume* fPhysDown = new G4PVPlacement(0,
                                                         G4ThreeVector(0,0,sideplace),
                                                         fLogicDown,
                                                         "ptfedown",
                                                         fLogicTPB,
                                                         false,
                                                         0,
                                                         true);
	new G4LogicalBorderSurface("refldown", fPhysTPB, fPhysDown, reflfoilOS);//*/

	
	lArFiducial = new G4PVPlacement(0,
					G4ThreeVector(0*cm, 0*cm, -0.0000*cm),
					fLogicFiduc,
					"Fiducial",
					fLogicTPB,
					false,
					0,
					true);

	//lAr to tpb optical surfaces, for both top/bottom and side layers
	new G4LogicalBorderSurface("refl", fPhysTPB, fPhysfoil, reflfoilOS);
 	new G4LogicalBorderSurface("refl2", fPhysTPBb, fPhysfoil, reflfoilOS);
 	new G4LogicalBorderSurface("reflhole", lArFiducial, fPhysfoil, reflfoilOS);

	new G4LogicalBorderSurface("TPB", lArFiducial, fPhysTPB, TPBfoilOS);
	new G4LogicalBorderSurface("TPB2", fPhysTPB, lArFiducial, TPBfoilOS);

	new G4LogicalBorderSurface("TPBb", lArFiducial, fPhysTPBb, TPBfoilOS);
	new G4LogicalBorderSurface("TPBb2", fPhysTPBb, lArFiducial, TPBfoilOS);//*/

	//adds in the layers of lAr if they are enabled. They are by default, but this can be changed.
	if (fLayers) {
	  lArFiducial1 = new G4PVPlacement(0,
					   G4ThreeVector(0*cm, 0*cm, row1place*cm),
					   fLogicFiduc1,
					   "Fiducial1",
					   fLogicFiduc,
					   false,
					   0,
					   true);
	  
	  lArFiducial2 = new G4PVPlacement(0,
					   G4ThreeVector(0*cm, 0*cm, 24.000*cm),
					   fLogicFiduc2,
					   "Fiducial2",
					   fLogicFiduc,
					   false,
					   0,
					   true);
	  
	  lArFiducial3 = new G4PVPlacement(0,
					   G4ThreeVector(0*cm, 0*cm, 0.000*cm),
					   fLogicFiduc3,
					   "Fiducial3",
					   fLogicFiduc,
					   false,
					   0,
					   true);
	  
	  lArFiducial4 = new G4PVPlacement(0,
					   G4ThreeVector(0*cm, 0*cm, -24.000*cm),
					   fLogicFiduc4,
					   "Fiducial4",
					   fLogicFiduc,
					   false,
					   0,
					   true);//*/
	  
	  lArFiducial5 = new G4PVPlacement(0, 
					   G4ThreeVector(0*cm, 0*cm, -1*row1place*cm),
					   fLogicFiduc5,
					   "Fiducial5",
					   fLogicFiduc,
					   false,
					   0,
					   true);//*/
	  
	  //Adds a small cone on the bottom near the center. makes the simulation better fit data and matches the LED in reality.
	  /*G4double bconewidth = conewide*cm;
	  G4double bconeheight = bconewidth*conehigh;
	  G4double bconeplace = bconeheight-totalH*cm;
	  
	  G4Cons* foilbcone = new G4Cons("ptfebcone", 0*cm, bconewidth, 0*cm, 0*cm, bconeheight, 0*deg, 360*deg);
	  G4LogicalVolume* fFoilBcone = new G4LogicalVolume(foilbcone, ptfe, "ptfebcone");
	  
	  G4Cons* tpbbcone = new G4Cons("tpbbcone", 0*cm, (bconewidth+thick*2.), 0*cm, 0*cm, (bconeheight+thick*2.), 0*deg, 360*deg);
	  G4LogicalVolume* fTPBBcone = new G4LogicalVolume(tpbbcone, tPB, "tpbbcone");
	  G4VPhysicalVolume* tpbBconePhys = new G4PVPlacement(0,
							      G4ThreeVector(0*cm, 0*cm, bconeplace),
							      fTPBBcone,
							      "tpbbcone",
							      fLogicFiduc5,
							      false,
							      0,
							      true);
	  
	  G4VPhysicalVolume* foilBconePhys = new G4PVPlacement(0,
							       G4ThreeVector(0*cm, 0*cm, -0.00015*cm),
							       fFoilBcone,
							       "ptfebcone",
							       fTPBBcone,
							       false,
							       0,
							       true);  
	  
	  //Layers require a lot of new border surfaces. Thus they are established.
	  new G4LogicalBorderSurface("reflbcone", tpbBconePhys, foilBconePhys, reflfoilOS);
	  
	  new G4LogicalBorderSurface("TPBbcone", lArFiducial5, tpbBconePhys, TPBfoilOS);
	  new G4LogicalBorderSurface("TPB2bcone", tpbBconePhys, lArFiducial5, TPBfoilOS);//*/
	  
	  new G4LogicalBorderSurface("1-TPB", lArFiducial1, fPhysTPB, TPBfoilOS);
	  new G4LogicalBorderSurface("1-TPB2", fPhysTPB, lArFiducial1, TPBfoilOS);
	  new G4LogicalBorderSurface("1-TPBb", lArFiducial1, fPhysTPBb, TPBfoilOS);
	  new G4LogicalBorderSurface("1-TPBb2", fPhysTPBb, lArFiducial1, TPBfoilOS);
	  new G4LogicalBorderSurface("2-TPB", lArFiducial2, fPhysTPB, TPBfoilOS);
	  new G4LogicalBorderSurface("2-TPB2", fPhysTPB, lArFiducial2, TPBfoilOS);
	  new G4LogicalBorderSurface("3-TPB", lArFiducial3, fPhysTPB, TPBfoilOS);
	  new G4LogicalBorderSurface("3-TPB2", fPhysTPB, lArFiducial3, TPBfoilOS);
	  new G4LogicalBorderSurface("4-TPB", lArFiducial4, fPhysTPB, TPBfoilOS);
	  new G4LogicalBorderSurface("4-TPB2", fPhysTPB, lArFiducial4, TPBfoilOS);//*/
	  new G4LogicalBorderSurface("5-TPB", lArFiducial5, fPhysTPB, TPBfoilOS);
	  new G4LogicalBorderSurface("5-TPB2", fPhysTPB, lArFiducial5, TPBfoilOS);//*/
	  new G4LogicalBorderSurface("5-TPBb", lArFiducial5, fPhysTPBb, TPBfoilOS);
	  new G4LogicalBorderSurface("5-TPBb2", fPhysTPBb, lArFiducial5, TPBfoilOS);//*/
	
	  //G4cout << "EdwardNote Defined Placed foil and fiducial" << G4endl;

	}//end fLayers if
	//G4cout << "Edwardnote Placed Optical Surfaces" << G4endl;

	//the below are a number of suboptions for if the reflector or tpb on the foils is turned off. Note from the format that they layered detector cannot be used if either tpb or reflector foils are missing. This is intentional; they make a best fit detector; neither of those things missing helps create a better fit.
      } else {
	//reflector On TPB off
	lArFiducial = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fLogicFiduc, "Fiducial",fLogicfoil,	false,	0,true);

        new G4LogicalBorderSurface("refl", lArFiducial, fPhysfoil, reflfoilOS);
      }

    } else if (fTPBfoilOn) {
      //Reflector Off TPB On
      fPhysTPB = new G4PVPlacement(0, G4ThreeVector(0,0,0), fLogicTPB, "TPBfoil",fLogicFrame, false, 0, true);

      lArFiducial = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fLogicFiduc, "Fiducial", fLogicTPB, false, 0, true);

      new G4LogicalBorderSurface("TPB", lArFiducial, fPhysTPB, TPBfoilOS);
      new G4LogicalBorderSurface("TPB2", fPhysTPB, lArFiducial, TPBfoilOS);
      
    } else {
      //places the lArFiducial Volume inside the Frame if both sets of foil are turned off.
      lArFiducial = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fLogicFiduc, "Fiducial", fLogicFrame, false, 0, true);
    }//end if-else for foils and fiducial volume



    //places the calibration rod (with extension) if its turned on.
    if (fRodin) {
      rodHeight = rodHeight + fifth*cm;//Shift the rod up or down.
      //G4cout << "Edwardnote Placing rod"  << G4endl;
      G4double rodLength = (totalH - 30.0*cm - rodHeight)/2.0;//remaining length of big rod in detector
      G4double mainplace = 30.0*cm+rodLength+rodHeight;//central position of big rod with respect to 0/center
      G4double knutplace = 32*cm+rodHeight;//position of the big knut around the end of the big rod
      G4double extendplace = 15*cm+rodHeight;//central position of the thin extension

      //check to make sure some of the rod is within the detector, otherwise do not place those volumes. 
      if (rodLength > 0) {
	//main rod
	G4Tubs* sourcerod = new G4Tubs("sourcerod", 0*cm, 2.5*cm, rodLength, 0*deg, 360*deg);
	G4LogicalVolume* fLogicRod = new G4LogicalVolume(sourcerod, steel, "sourcerod");
	new G4PVPlacement(0,
			  G4ThreeVector(0*cm, 0*cm, mainplace),//z=75*cm for extension,45 without
			  //G4ThreeVector(0*cm, 0*cm, 90*cm),//for upper position 18" above middle
			  fLogicRod,
			  "sourcerod",
			  fLogicFiduc,
			  false,
			  0,
			  true);

	/*//Clouding around the main source rod
	fFiducialAr1 = new G4Tubs("Fiducial1", 3*cm, 10*cm, rodLength, 0*deg, 360*deg);
	fLogicFiduc1 = new G4LogicalVolume(fFiducialAr1,lAr2,"Fiducial1");
	new G4PVPlacement(0,
			  G4ThreeVector(0*cm, 0*cm, mainplace),//z=75*cm for extension,45 without
			  //G4ThreeVector(0*cm, 0*cm, 90*cm),//for upper position 18" above middle
			  fLogicFiduc1,
			  "Fiducial1",
			  fLogicFiduc,
			  false,
			  0,
			  true);//*/
	
	//lAr inside the rod (it's hollow).
	G4Tubs* hollowrod = new G4Tubs("hollowrod", 0*cm, 2.*cm, 44*cm, 0*deg, 360*deg);
	G4LogicalVolume* fHoleRod = new G4LogicalVolume(hollowrod, lAr, "hollowrod");
	new G4PVPlacement(0,
			  G4ThreeVector(0*cm, 0*cm, 0*cm),
			  fHoleRod,
			  "hollowrod",
			  fLogicRod,
			  false,
			  0,
			  true);
	
	//endcap (that very large knut)
	G4Tubs* sourcecap = new G4Tubs("sourcecap", 2.5*cm, 3*cm, 2*cm, 0*deg, 360*deg);
	G4LogicalVolume* fLogicCap = new G4LogicalVolume(sourcecap, steel, "sourcecap");
	new G4PVPlacement(0,
			  G4ThreeVector(0*cm, 0*cm, knutplace),//z=32*cm for extension, 2 without
			  //G4ThreeVector(0*cm, 0*cm, 47*cm),//for upper position
			  fLogicCap,
			  "sourcecap",
			  fLogicFiduc,
			  false,
			  0,
			  true);
      }
      
      //the thinner extension
      G4Tubs* sourceextend = new G4Tubs("sourceextend", 0*cm, 0.5*cm, 15*cm, 0*deg, 360*deg);
      G4LogicalVolume* fLogicExtend = new G4LogicalVolume(sourceextend, steel, "sourceextend");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, extendplace),
			fLogicExtend,
			"sourceextend",
			fLogicFiduc,
			false,
			0,
			true);
      
      //extension is also hollow, or at least approximated as such.
      G4Tubs* hollowextend = new G4Tubs("hollowextend", 0*cm, 0.25*cm, 14.8*cm, 0*deg, 360*deg);
      G4LogicalVolume* fHoleExtend = new G4LogicalVolume(hollowextend, lAr, "hollowextend");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, 0.2*cm),
			fHoleExtend,
			"hollowextend",
			fLogicExtend,
			false,
			0,
			true);//*/

      //G4cout << "Edwardnote Placed rod"  << G4endl;
      /*//Clouding around the bottom
      fFiducialAr1 = new G4Tubs("Fiducial1", 1*cm, ultra*cm, 7.5*cm, 0*deg, 360*deg);
      fLogicFiduc1 = new G4LogicalVolume(fFiducialAr1,lAr2,"Fiducial1");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, -47.5*cm),//z=75*cm for extension,45 without
			//G4ThreeVector(0*cm, 0*cm, 90*cm),//for upper position 18" above middle
			fLogicFiduc1,
			"Fiducial1",
			fLogicFiduc,
			false,
			0,
			true);//*/

    } else if (fLaser) {
      //creates a cone of foil around the top of the laser rod, created by the loose foils there not being allowed to return to rest. Not significant in the source calibrations, though its inclusion won't hurt either.
      G4double coneheight = 6.7*cm;
      if (rodHeight > 30*cm) {
	coneheight = 3.87*cm;
      }
      G4double coneplace = totalH*cm-coneheight;
      G4double topwide = 15.2*cm;
      
      G4Cons* foilcone = new G4Cons("ptfecone", 1.8*cm, 1.81*cm, 1.8*cm, topwide, coneheight, 0*deg, 360*deg);
      G4LogicalVolume* fFoilCone = new G4LogicalVolume(foilcone, ptfe, "ptfecone");

      G4Cons* tpbcone = new G4Cons("tpbcone", 1.8*cm, (1.810*cm+(thick)), 1.8*cm, (15.200*cm+(thick)), (coneheight+(thick/2.)), 0*deg, 360*deg);
      G4LogicalVolume* fTPBCone = new G4LogicalVolume(tpbcone, tPB, "tpbcone");

      //place the cones of foil (ptfe) and tpb into the appropriate fiducial volume
      G4VPhysicalVolume* tpbConePhys;
      if (fLayers) {
	coneplace = coneplace - 52.0*cm;
	tpbConePhys = new G4PVPlacement(0,
					G4ThreeVector(0*cm, 0*cm, coneplace),
					fTPBCone,
					"tpbcone",
					fLogicFiduc1,
					false,
					0,
					true);
      } else {
	tpbConePhys = new G4PVPlacement(0,
					G4ThreeVector(0*cm, 0*cm, coneplace),
					fTPBCone,
					"tpbcone",
					fLogicFiduc,
					false,
					0,
					true);
      }

      G4VPhysicalVolume* foilConePhys = new G4PVPlacement(0,
							 G4ThreeVector(0*cm, 0*cm, 0*cm),
							 fFoilCone,
							 "tpbcone",
							 fTPBCone,
							 false,
							 0,
							 true);  

      //optical surfaces for the cone of foil
      new G4LogicalBorderSurface("reflcone", tpbConePhys, foilConePhys, reflfoilOS);
      
      if (fLayers) {
	new G4LogicalBorderSurface("TPBcone", lArFiducial1, tpbConePhys, TPBfoilOS);
	new G4LogicalBorderSurface("TPB2cone", tpbConePhys, lArFiducial1, TPBfoilOS);
      } else {
	new G4LogicalBorderSurface("TPBcone", lArFiducial, tpbConePhys, TPBfoilOS);
	new G4LogicalBorderSurface("TPB2cone", tpbConePhys, lArFiducial, TPBfoilOS);
      }

      //start creating and placing the laser rod, with variable height to be input by the user (default: 0*cm, center of detector).
      G4double rodplace = (61.79*cm+rodHeight);
      G4Tubs* sourcerod = new G4Tubs("sourcerod", 0*cm, 1.8*cm, 60*cm, 0*deg, 360*deg);
      G4LogicalVolume* fLogicRod = new G4LogicalVolume(sourcerod, steel, "sourcerod");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, rodplace),
			fLogicRod,
			"sourcerod",
			fLogicFiduc,
			false,
			0,
			true);

      //this rod is also hollow
      G4Tubs* hollowrod = new G4Tubs("hollowrod", 0*cm, 1.6*cm, 60*cm, 0*deg, 360*deg);
      G4LogicalVolume* fHoleRod = new G4LogicalVolume(hollowrod, lAr, "hollowrod");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, 0*cm),
			fHoleRod,
			"hollowrod",
			fLogicRod,
			false,
			0,
			true);
      
      //creates a cylinder with a sort of hollowed cone for the diffuser (geometry measurements from diffuser design diagram)
      G4Tubs* diffuserBody = new G4Tubs("diffusermain", 0*cm, 1.905*cm, 1.79*cm, 0*deg, 360*deg);
      G4Cons* diffuserCone = new G4Cons("diffusercone", 0*cm, 1.409*cm, 0*cm, .952*cm, .394*cm, 0*deg, 360*deg);
      G4SubtractionSolid* diffuser = new G4SubtractionSolid("diffuser", diffuserBody, diffuserCone, 0, G4ThreeVector(0*cm, 0*cm, -1.396*cm));
      G4LogicalVolume* fLogicDiff = new G4LogicalVolume(diffuser, steel, "diffuser");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, rodHeight),
			fLogicDiff,
			"diffuser",
			fLogicFiduc,
			false,
			0,
			true);
      
      //creates the diffuser plates out of glass (just in case you accidentally put the source of the photons inside them). 
      G4Tubs* plates = new G4Tubs("plates", 0*cm, 1.33*cm, .305*cm, 0*deg, 360*deg);
      G4LogicalVolume* logicPlates = new G4LogicalVolume(plates, fGlass, "plates");
      new G4PVPlacement(0,
			G4ThreeVector(0*cm, 0*cm, -.70*cm),
			logicPlates,
			"plates",
			fLogicDiff,
			false,
			0,
			true);
    }//end laser placement if statement.

    //G4cout << "Edwardnote after rod/laser ifelse"  << G4endl;

    //create a loop to place PMTs according to pattern
    G4double angle = 0.;
    G4double pmtxx;
    G4double pmtyy;
    G4double pmtzz;
    G4String pmtnm[120];//each pmt gets a unique name so they can be distinguished in the output.
    G4bool pmtcoat;

    G4double itang = 360.0/24;//interval angle; angle between two consecutive pmts
    G4double radius = 96.0;//radius: outer radius of Fiducial volume
    G4double pmtRad = 10.2;//radius of the PMT
    radius = radius + (pmtRad-topthick);
    G4int zs = 0;//counter

    //G4cout << "Placing PMTs" << G4endl;

    //begin loop for creation of PMT names
    for (G4int i=0; i<5; ++i) {
      for (G4int n=0; n<24; ++n) {
	zs = i*24+n;
	pmtnm[zs] = "C"+std::to_string(n+1)+"R"+std::to_string(i+1);
	//now all pmts are named according to column and row

	//G4cout << pmtnm[zs];
      }
    }
    zs = 0;

    for (G4int i=0; i<120; ++i) {
      if (i%24 == 0) {
	angle = 0;
      }
      //defines the angular position according to the interval angle
      G4double phi = angle*CLHEP::pi/180.0;
      pmtxx = radius*std::cos(phi);
      pmtyy = radius*std::sin(phi);

      //define the pmtzz value to translate row 1 (i/24=0) to the top and row 5 (i/24=4) to the bottom
      zs = 2-(i/24);
      pmtzz = zs*21.86+1.0;//totalH = 61.6, PMTrad = 10.2
      
      //iterate the angle according to the interval
      angle += itang;

      //define if PMTs are coated or uncoated.
      if (ccm200) {
	//start all pmts as coated and uncoated the appropriate ones
	pmtcoat = true;
	//Uncoated: Row 1: 2,6,10,14,18,22. Row2: 4,8,12,16,20,24. Row4: 1,5,9,13,17,21. Row5: 3,7,11,15,19,23.
	//Iterator %24 = column-1.
	if ((i%24)%4 == 1 && zs == 2)  { pmtcoat = false; }
	if ((i%24)%4 == 3 && zs == 1)  { pmtcoat = false; }
	if ((i%24)%4 == 0 && zs == -1) { pmtcoat = false; }
	if ((i%24)%4 == 2 && zs == -2) { pmtcoat = false; }
      }
      else {
	//start all pmts as coated and uncoated the selected.
	pmtcoat = true;
	//uncoated the odds on row 4 and the evens on row 2. Remember iterator is column-1. 
        if ((i%24)%2 != 0 && zs == 1)  {  pmtcoat = false;	}
	if ((i%24)%2 == 0 && zs == -1) {  pmtcoat = false;      }
      }
      //for universal, uncomment this and make all pmts whatever you want.
      //pmtcoat=false;
      if (cylinderOn || fInvBeta){//old flag for turning all pmts to uncoated for the LBOC design
	pmtcoat=false;
      }
      if (pmtcoat) {
	pmtnm[i] = pmtnm[i]+"_coated";
        //if coated, add coated to the name. Makes things much easier down the line.	
      }
      if (fPMTsOn){
	placePMT(pmtnm[i],pmtxx,pmtyy,pmtzz,phi,pmtcoat);
	//call the placePMT method to actually place the PMT at the chosen position with the chosen name and coating.
      }
    }

    //now do another set of loops for the top and bottom pmts if CCM200 is being used.
    if (ccm200) {      

      //define outer ring of pmts for top and bottom. This format was chosen to make copy-pasting easier
      G4String pmtnam;
      G4String pmtnam1;//new names
      G4int npmts = 20;//number of pmts in the ring
      itang = 360./npmts;//new iterating angle
      radius = 84.2;
      angle = 0;
      
      for (G4int n=0; n<npmts; ++n) {
	if (n%2 == 0) {
	  zs = 0;
	  pmtzz = +1.0*totalH;
	} else { 
	  zs = 6;
	  pmtzz = -1.0*totalH;//NEW VARIABLE: DISPLACEMENT
	}//let 6 be bottom and 0 be top, alternating around the circle in every other position from the design of CCM220
	
	//Uncoated on cross1, ring 4: 5, 15, 8, 18
	if (n==4 || n==14 || n==7 || n==17) { 
	  pmtcoat = false;
	} else { pmtcoat = true; }

	pmtnam = "C"+std::to_string(400+n+1)+"R"+std::to_string(zs);
	//define a new set of names, with the rows being R0x for the top (0 to indicate top, x for the circle there) and simuilarly R6x for the bottom.
	
	//define the position and coating and place the PMTs
	angle = n*itang;
	G4double phi = angle*CLHEP::pi/180.0;
	pmtxx = radius*std::cos(phi);
	pmtyy = radius*std::sin(phi);
	//pmtcoat = true;
	if (cylinderOn || fInvBeta) { pmtcoat=false; }//once more, turn all coatings off if using the LBOC design
	placeTopBot(pmtnam,pmtxx,pmtyy,pmtzz,pmtcoat);
      }

      //define the second ring of pmts on top and bottom. identical to the outer row, but with a slightly smaller circle
      npmts = 15;
      itang = 360./npmts;
      radius = 62.8;
      angle = 0;

      for (G4int n=0; n<npmts; ++n) {
	pmtnam = "C"+std::to_string(300+n+1)+"R6";
	pmtnam1 = "C"+std::to_string(300+n+1)+"R0";
	//G4cout << pmtnam << pmtnam1;
	
	angle = n*itang;
	G4double phi = angle*CLHEP::pi/180.0;
	pmtxx = radius*std::cos(phi);
	pmtyy = radius*std::sin(phi);
	pmtcoat = true;
	
	//if (n%3 == 0) { pmtcoat = false; }//every third pmt in the second row is uncoated. 10 total
	if (cylinderOn || fInvBeta) { pmtcoat=false; }

	pmtzz = -1*totalH;
	//ring3 uncoated: top 301, 308 bottom 303, 310.
	if (n==2 || n==9) { pmtcoat=false;}
	placeTopBot(pmtnam,pmtxx,pmtyy,pmtzz,pmtcoat);
	pmtzz = totalH;
	if (n==0 || n==7) { 
	  pmtcoat=false;
	} else { pmtcoat=true; }
	placeTopBot(pmtnam1,pmtxx,pmtyy,pmtzz,pmtcoat);
      }

      //third row on top and bottom. again smaller circle and npmts.
      npmts = 10;
      itang = 360./npmts;
      radius = 41.4;
      angle = 0;

      for (G4int n=0; n<npmts; ++n) {
	pmtnam = "C"+std::to_string(200+n+1)+"R6";
	pmtnam1 = "C"+std::to_string(200+n+1)+"R0";
	//G4cout << pmtnam << pmtnam1;
	
	angle = n*itang;
	G4double phi = angle*CLHEP::pi/180.0;
	pmtxx = radius*std::cos(phi);
	pmtyy = radius*std::sin(phi);
	pmtcoat = true;
	//if (n%3 == 2) { pmtcoat = false; }//every second pmt in the third row is uncoated, for another 10 total and thus all 20 uncoated on the top and bottom accounted for
	if (cylinderOn || fInvBeta) { pmtcoat=false; }

	//ring 2 uncoated: top 208, 203, bottom 205, 210
	pmtzz = -1*totalH;
	if (n==4 || n==9) { pmtcoat=false;}
	placeTopBot(pmtnam,pmtxx,pmtyy,pmtzz,pmtcoat);
	pmtzz = totalH;
	if (n==2 || n==7) { 
	  pmtcoat=false;
	} else { pmtcoat=true; }
	placeTopBot(pmtnam1,pmtxx,pmtyy,pmtzz,pmtcoat);
      }

      //fourth and final row on the top and bottom
      npmts = 5;
      itang = 360./npmts;
      radius = 20.0;
      angle = 0;

      for (G4int n=0; n<npmts; ++n) {
	pmtnam = "C"+std::to_string(100+n+1)+"R6";
	pmtnam1 = "C"+std::to_string(100+n+1)+"R0";
	//G4cout << pmtnam << pmtnam1;
	
	angle = n*itang;
	G4double phi = angle*CLHEP::pi/180.0;
	pmtxx = radius*std::cos(phi);
	pmtyy = radius*std::sin(phi);
	pmtcoat = true;
	if (cylinderOn || fInvBeta) { pmtcoat=false; }

	//ring 1 uncoated: top 101, 103, bottom 102, 104
	pmtzz = -1*totalH;
	if (n==1 || n==3) { pmtcoat=false;}
	placeTopBot(pmtnam,pmtxx,pmtyy,pmtzz,pmtcoat);
	pmtzz = totalH;
	if (n==0 || n==2) { 
	  pmtcoat=false;
	} else { pmtcoat=true; }
	placeTopBot(pmtnam1,pmtxx,pmtyy,pmtzz,pmtcoat);
      }
    }//*/
    
    //G4cout << "Edwardnote Placed PMTs" << G4endl;

  }//end place main volumes

  //G4cout << "Edwardnote Main Volume Done" << G4endl;

  return wPhys;
}

//Class for creating and placing PMT half sphere facing inwards (for the side 120 pmts)

void detectorConstruction::placePMT(G4String name, 
				    G4double pmt_x, G4double pmt_y, G4double pmt_z, 
				    G4double phi, G4bool coated) {

  //defines total detector height
  G4double totalH = 136.0/2.0;
  if (ccm200) { totalH = 123.2/2.0; }

  G4double remainder = totalH/2.0-0.005;
  G4double row1place = (totalH+36.0)-remainder;//134 -> 52; 123.2 -> 

  //reset the angles to 0 for each new PMT
  G4double init_angle = 0;
  G4double fin_angle = 0;
  G4double vert_angle = 0;

  //sets the pmts to half spheres; 180 degrees from initial angle to final angle horizontally and vertically.
  if (mainVolOn){
    fin_angle = 180*deg;
    vert_angle= 180*deg;
  } else {
    fin_angle = 360*deg;
    vert_angle= 90*deg;
  }

  //set the outer radii of the glass and the thickness of the TPB.
  G4double radout = 10.2;
  G4double tpbout = 0.00009*cm;//In the best fit, the TPB on the PMTs is slightly thinner than normal. Normal:0.00019*cm.
  G4double frillin = std::sqrt(radout*radout-(radout-topthick)*(radout-topthick))*cm;
  //G4cout << "Frill inner radius in mm: " << frillin << "  details: " << radout << '\t' << topthick << '\t' << radout*radout << '\t' << (radout-topthick)*(radout*topthick) << '\t' << radout*radout-(radout-topthick)*(radout-topthick) << G4endl;

  tpbout = tpbout*conehigh;///randwide;
  radout = radout*cm;

  //define the initial angle so the pmt is facing inwards according to the position
  init_angle = (std::atan(pmt_y/pmt_x)*180/CLHEP::pi);
  if (pmt_x < 0) {
    init_angle+=180;
  }
  init_angle = init_angle*deg;
  
  //defines the name of the PMT so it contains PMT
  G4String namePMT = "PMT_"+name;
  G4String frilName = "frilblack_"+name;
  G4String frillName = "frill_"+name;
  
  //defines a three vector for the position.
  G4ThreeVector trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);

  G4double frillx = 95.9*std::cos(phi);
  G4double frilly = 95.9*std::sin(phi);
  G4ThreeVector transf = G4ThreeVector(frillx*cm,frilly*cm,pmt_z*cm);
  G4ThreeVector xfrill = G4ThreeVector(0,0,1);
  G4ThreeVector yfrill = G4ThreeVector(0,1,0);
  G4ThreeVector zfrill = G4ThreeVector(-1,0,0);
  G4RotationMatrix* rotationMatrix = new G4RotationMatrix(xfrill,yfrill,zfrill);
  rotationMatrix->rotateZ(init_angle);

  init_angle = init_angle + 90*deg;
  //creates the glass sphere properly oriented for the position
  G4Sphere* pmt = new G4Sphere(namePMT, 0*cm, radout, init_angle, fin_angle, 0*deg, vert_angle);
  G4LogicalVolume* pmtLogic = new G4LogicalVolume(pmt, fGlass, namePMT);

  //create the frill around the base of the PMT, in 2 parts: 1 wider in thin plastic, the other thin but thick black plastic.
  G4double ft = 0.2*cm;
  G4double ftt = 1.0*cm*1.0+ft;
  G4double fh = 1.0*cm*0.2;
  G4double fhh = 0.1*cm;
    
  G4Tubs* fril = new G4Tubs(frilName, frillin+tpbout, frillin+ft, fh, 0*deg, 360*deg);
  G4Tubs* frill = new G4Tubs(frillName, frillin+ft, frillin+ftt, fhh, 0*deg, 360*deg);
  G4LogicalVolume* frilLogic = new G4LogicalVolume(fril,fBlackPlastic,frilName);//*/
  G4LogicalVolume* frillLogic = new G4LogicalVolume(frill,fPlastic,frillName);//*/

  //begin defining the values for the tpb or ice optical surface (for some reason it won't pull it from the Construct method all the time).
  const G4int nAcTefEntries = 15;
  G4double TPBOSEnergy[nAcTefEntries] =
    {0.602*eV, 1.03*eV,  2.25*eV,  2.48*eV,  2.95*eV,
     3.457*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV,
     6.262*eV, 7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
  G4double TPBOSTransmit[nAcTefEntries] =
    {1., 1., 1., 1., 1.,
     1., 1., 1., 1., 1.,
     1., 1., 1., 1., 1.}; //set to 1 and have all absorption in bulk 
  G4double TPBOSReflect[nAcTefEntries] =
    {0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0};
  G4double ref = 0.1;//conehigh;
  G4double FrillOSReflect[nAcTefEntries] =
    {ref, ref, ref, ref, ref,
     ref, ref, ref, ref, ref,
     ref, ref, ref, ref, ref};//*/
  G4double TPBOSEff[nAcTefEntries] =
    {1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0};

  G4OpticalSurface* frillCoat = new G4OpticalSurface("FRILL Coat");
  frillCoat->SetModel(unified);
  frillCoat->SetType(dielectric_dielectric);
  frillCoat->SetFinish(ground);
  frillCoat->SetSigmaAlpha(0.05);
  
  G4MaterialPropertiesTable *FRILLMPT = new G4MaterialPropertiesTable();
  FRILLMPT->AddProperty("REFLECTIVITY", TPBOSEnergy, FrillOSReflect, nAcTefEntries);
  FRILLMPT->AddProperty("TRANSMITTANCE", TPBOSEnergy, TPBOSTransmit, nAcTefEntries);
  FRILLMPT->AddProperty("EFFICIENCY", TPBOSEnergy, TPBOSEff, nAcTefEntries);
  frillCoat->SetMaterialPropertiesTable(FRILLMPT);//*/

  if (coated) {
    //rename with coated if the pmt is coated and create a larger semisphere for the TPB
    namePMT = "PMT_"+name+"_coated";
    G4String namecoat = "tpbcoat_"+name;
    G4Sphere* tpbcoating = new G4Sphere("tpbcoat", 0*cm, (radout+tpbout), init_angle, fin_angle, 0*deg, vert_angle);//define pmttpb coating thickness
    G4LogicalVolume* tpbLogic = new G4LogicalVolume(tpbcoating, tPBhundred, namecoat);

    //define the TPB optical properties table
    G4OpticalSurface* tpbCoat = new G4OpticalSurface("TPB Coat");
    tpbCoat->SetModel(unified);
    tpbCoat->SetType(dielectric_dielectric);
    tpbCoat->SetFinish(ground);
    tpbCoat->SetSigmaAlpha(0.05);

    G4MaterialPropertiesTable *TPBMPT = new G4MaterialPropertiesTable();
    TPBMPT->AddProperty("REFLECTIVITY", TPBOSEnergy, TPBOSReflect, nAcTefEntries);
    TPBMPT->AddProperty("TRANSMITTANCE", TPBOSEnergy, TPBOSTransmit, nAcTefEntries);
    TPBMPT->AddProperty("EFFICIENCY", TPBOSEnergy, TPBOSEff, nAcTefEntries);
    tpbCoat->SetMaterialPropertiesTable(TPBMPT);//*/

    //if Layered detector, place the TPB spheres in the proper fiducial volume layer
    if (fLayers) {
      if (pmt_z > 40) {
	pmt_z = pmt_z-row1place;
	trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);
	
	G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
							   trans,
							   tpbLogic,
							   namecoat,
							   fLogicFiduc1,
							   false,
							   0,
							   true);
	
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial1,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial1,
				   tpbCoat);
	
      } else if (pmt_z > 20) {
	pmt_z = pmt_z-24.0;
	trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);
	
	G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
							   trans,
							   tpbLogic,
							   namecoat,
							   fLogicFiduc2,
							   false,
							   0,
							   true);
	
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial2,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial2,
				   tpbCoat);
	
      } else if (pmt_z > -40 && pmt_z < -20) {
	pmt_z = pmt_z+24.0;
	trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);
	
	G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
							   trans,
							   tpbLogic,
							   namecoat,
							   fLogicFiduc4,
							   false,
							   0,
							   true);
	
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial4,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial4,
				   tpbCoat);
	
      } else if (pmt_z == 0) {
	G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
							   trans,
							   tpbLogic,
							   namecoat,
							   fLogicFiduc3,
							   false,
							   0,
							   true);
      
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial3,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial3,
				   tpbCoat);
	//*/
      } else {//*/
	pmt_z = pmt_z+row1place;
	trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);

	G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
							   trans,
							   tpbLogic,
							   namecoat,
							   fLogicFiduc5,
							   false,
							   0,
							   true);
	
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial5,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial5,
				   tpbCoat);
      }//*/
    } else {//*/
      //place the tpb half sphere and the optical surfaces with the Fiducial volume
      G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
							 trans, //G4Transform3D(*rotationMatrix,trans),
							 tpbLogic,
							 namecoat,
							 fLogicFiduc,
							 false,
							 0,
							 false);
      
      new G4LogicalBorderSurface("TPB Coat", 
				 lArFiducial,
				 tpbCoating,
				 tpbCoat);
      new G4LogicalBorderSurface("TPB2 Coat", 
				 tpbCoating,
				 lArFiducial,
				 tpbCoat);//*/
    }//*/
    //place the glass PMT inside the tpb semi-sphere
    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      pmtLogic,
		      namePMT,
		      tpbLogic,
		      false,
		      0,
		      false);

    G4VPhysicalVolume*  frillPhys = new G4PVPlacement(G4Transform3D(*rotationMatrix,transf),
						      frillLogic,
						      frillName,
						      fLogicFiduc,
						      false,
						      0,
						      false);
		      
    G4VPhysicalVolume*  frilPhys = new G4PVPlacement(G4Transform3D(*rotationMatrix,transf),
						     frilLogic,
						     frilName,
						     fLogicFiduc,
						     false,
						     0,
						     false);
    
    new G4LogicalBorderSurface("Frill Coat", frillPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", lArFiducial, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Frill Coat", frillPhys, fPhysTPB, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", fPhysTPB, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Fril Coat", frilPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Fril2 Coat", lArFiducial, frilPhys, frillCoat);
		      
    
  } else if (!coated && fLayers) {
    //now place the uncoated PMTs with a layered detector
    if (pmt_z > 40) {
      pmt_z = pmt_z-52.0;
      trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);

      new G4PVPlacement(0,
			trans,
			pmtLogic,
			namePMT,
			fLogicFiduc1,
			false,
			0,
			true);
      
    } else if (pmt_z > 20) {
      pmt_z = pmt_z-24.0;
      trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);

      new G4PVPlacement(0,
			trans,
			pmtLogic,
			namePMT,
			fLogicFiduc2,
			false,
			0,
			true);
      
    } else if (pmt_z > -40 && pmt_z < -20) {
      pmt_z = pmt_z+24.0;
      trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);

      new G4PVPlacement(0,
			trans,
			pmtLogic,
			namePMT,
			fLogicFiduc4,
			false,
			0,
			true);
      
    } else if (pmt_z == 0) {
      new G4PVPlacement(0,
			trans,
			pmtLogic,
			namePMT,
			fLogicFiduc3,
			false,
			0,
			true);
      
    } else {//*/
      pmt_z = pmt_z+51.9;
      trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);

      new G4PVPlacement(0,
			trans,
			pmtLogic,
			namePMT,
			fLogicFiduc5,
			false,
			0,
			true);
      
    }//*/
  } else if (!coated && !fLayers) {
    //place the uncoated PMT in the fiducial volume if uncoated and there are no layers
    new G4PVPlacement(0,
		      trans, //G4Transform3D(*rotationMatrix,trans),
		      pmtLogic,
		      namePMT,
		      fLogicFiduc,
		      false,
		      0,
		      false);
    
    G4VPhysicalVolume*  frillPhys = new G4PVPlacement(G4Transform3D(*rotationMatrix,transf),
						      frillLogic,
						      frillName,
						      fLogicFiduc,
						      false,
						      0,
						      false);
		      
    G4VPhysicalVolume*  frilPhys = new G4PVPlacement(G4Transform3D(*rotationMatrix,transf),
						     frilLogic,
						     frilName,
						     fLogicFiduc,
						     false,
						     0,
						     false);
    
    new G4LogicalBorderSurface("Frill Coat", frillPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", lArFiducial, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Frill Coat", frillPhys, fPhysTPB, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", fPhysTPB, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Fril Coat", frilPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Fril2 Coat", lArFiducial, frilPhys, frillCoat);
		      
  }//*/

  //G4cout << "PMTplacement made: " << namePMT << G4endl;

}

//create and place the PMT half spheres on the top and bottom.

void detectorConstruction::placeTopBot(G4String name, 
				       G4double pmt_x, G4double pmt_y, G4double pmt_z, 
				       G4bool coated) {
  //reset the angles to 0 for each new PMT
  G4double init_angle = 0;
  G4double fin_angle = 0;
  G4double vert_angle = 0;

  //define the vertical angle so the PMT is facing either up or down.
  if (pmt_z > 50) {
    vert_angle = 90*deg;
  } else {
    vert_angle = 0*deg;
  }

  G4double radius = pmt_x*pmt_x+pmt_y*pmt_y;

  fin_angle = 360*deg;//all top/bottom pmts should be full circles in the horizontal.

  //set the radii of the glass and the tpb thickness
  G4double radout = 10.2;
  G4double tpbout = 0.00009*cm;
  G4double frillin = std::sqrt(radout*radout-(radout-topthick)*(radout-topthick))*cm;

  tpbout = tpbout*conehigh;//
  G4double newz = pmt_z + (radout-topthick);
  if (pmt_z < -50) {
    newz = pmt_z - (radout-topthick);
    tpbout = tpbout;//
    pmt_z = pmt_z + 0.1*cm;
  }
  else {  pmt_z = pmt_z - 0.1*cm; }
  radout = radout*cm;

  //add PMT to the name
  G4String namePMT = "PMT_"+name;
  G4String frillName = "frill_"+name;
  G4String frilName = "frilblack_"+name;
  
  //define the position three-vector
  G4ThreeVector transf = G4ThreeVector(pmt_x*cm,pmt_y*cm,pmt_z*cm);
  G4ThreeVector trans = G4ThreeVector(pmt_x*cm,pmt_y*cm,newz*cm);
  //G4cout << "topbottom PMT position = " << pmt_x << '\t' << pmt_y << '\t' << pmt_z << G4endl;

  //create the frill around the base of the PMT, in 2 parts: 1 wider in thin plastic, the other thin but thick black plastic.
  G4double ft = 0.2*cm;
  G4double ftt = 1.0*cm*1.0+ft;
  G4double fh = 1.0*cm*0.2;
  G4double fhh = 0.1*cm;
    
  G4Tubs* fril = new G4Tubs(frilName, frillin+tpbout, frillin+ft, fh, 0*deg, 360*deg);
  G4Tubs* frill = new G4Tubs(frillName, frillin+ft, frillin+ftt, fhh, 0*deg, 360*deg);
  G4LogicalVolume* frilLogic = new G4LogicalVolume(fril,fBlackPlastic,frilName);//*/
  G4LogicalVolume* frillLogic = new G4LogicalVolume(frill,fPlastic,frillName);//*/

  //redo the tpb optical surface properties.
  const G4int nAcTefEntries = 15;
  G4double TPBOSEnergy[nAcTefEntries] =
    {0.602*eV, 1.03*eV,  2.25*eV,  2.48*eV,  2.95*eV,
     3.457*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV,
     6.262*eV, 7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
  G4double TPBOSTransmit[nAcTefEntries] =
    {1., 1., 1., 1., 1.,
     1., 1., 1., 1., 1.,
     1., 1., 1., 1., 1.}; //set to 1 and have all absorption in bulk 
  G4double TPBOSReflect[nAcTefEntries] =
    {0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0, 0.0};
  G4double ref = 0.1;//conehigh;
  G4double FrillOSReflect[nAcTefEntries] =
    {ref, ref, ref, ref, ref,
     ref, ref, ref, ref, ref,
     ref, ref, ref, ref, ref};//*/
  G4double TPBOSEff[nAcTefEntries] =
    {0.0, 0.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0,
     1.0, 1.0, 1.0, 1.0, 1.0};

  G4OpticalSurface* frillCoat = new G4OpticalSurface("FRILL Coat");
  frillCoat->SetModel(unified);
  frillCoat->SetType(dielectric_dielectric);
  frillCoat->SetFinish(ground);
  frillCoat->SetSigmaAlpha(0.05);
  
  G4MaterialPropertiesTable *FRILLMPT = new G4MaterialPropertiesTable();
  FRILLMPT->AddProperty("REFLECTIVITY", TPBOSEnergy, FrillOSReflect, nAcTefEntries);
  FRILLMPT->AddProperty("TRANSMITTANCE", TPBOSEnergy, TPBOSTransmit, nAcTefEntries);
  FRILLMPT->AddProperty("EFFICIENCY", TPBOSEnergy, TPBOSEff, nAcTefEntries);
  frillCoat->SetMaterialPropertiesTable(FRILLMPT);//*/

  //create a properly oriented glass half sphere for the PMT
  G4Sphere* pmt = new G4Sphere(namePMT, 0*cm, radout, init_angle, fin_angle, vert_angle, 90*deg);
  G4LogicalVolume* pmtLogic = new G4LogicalVolume(pmt, fGlass, namePMT);

  //if coated, add the coating and place the PMT
  if (coated) {
    namePMT = "PMT_"+name+"_coated";
    G4String namecoat = "tpbcoat_"+name;
    G4Sphere* tpbcoating = new G4Sphere(namecoat, 0*cm, (radout+tpbout), init_angle, fin_angle, vert_angle, (vert_angle+90*deg));//define pmttpb coating thickness
    G4LogicalVolume* tpbLogic = new G4LogicalVolume(tpbcoating, tPBhundred, namecoat);
    G4OpticalSurface* tpbCoat = new G4OpticalSurface("TPB Coat");
    tpbCoat->SetModel(unified);
    tpbCoat->SetType(dielectric_dielectric);
    tpbCoat->SetFinish(ground);
    tpbCoat->SetSigmaAlpha(0.05);

    G4MaterialPropertiesTable *TPBMPT = new G4MaterialPropertiesTable();
    TPBMPT->AddProperty("REFLECTIVITY", TPBOSEnergy, TPBOSReflect, nAcTefEntries);
    TPBMPT->AddProperty("TRANSMITTANCE", TPBOSEnergy, TPBOSTransmit, nAcTefEntries);
    TPBMPT->AddProperty("EFFICIENCY", TPBOSEnergy, TPBOSEff, nAcTefEntries);
    tpbCoat->SetMaterialPropertiesTable(TPBMPT);//*/

    //place the PMT in the fiducial volume even if Layers are on. Note: this works, but with some errors due to overlaps of daughter volumes with the same mother. Thus it is not recommended to use the layered detector and ccm200 at the same time.
    G4VPhysicalVolume*  tpbCoating = new G4PVPlacement(0,
						       trans,
						       tpbLogic,
						       namecoat,
						       fLogicFiduc,
						       false,
						       0,
						       true);
	
    G4VPhysicalVolume*  frillPhys = new G4PVPlacement(0,
						      transf,
						      frillLogic,
						      frillName,
						      fLogicFiduc,
						      false,
						      0,
						      false);
		      
    G4VPhysicalVolume*  frilPhys = new G4PVPlacement(0,
						     transf,
						     frilLogic,
						     frilName,
						     fLogicFiduc,
						     false,
						     0,
						     false);

    new G4LogicalBorderSurface("Frill Coat", frillPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", lArFiducial, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Frill Coat", frillPhys, fPhysTPBb, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", fPhysTPBb, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Fril Coat", frilPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Fril2 Coat", lArFiducial, frilPhys, frillCoat);      

    new G4PVPlacement(0,
		      G4ThreeVector(0,0,0),
		      pmtLogic,
		      namePMT,
		      tpbLogic,
		      false,
		      0,
		      true);
    
    new G4LogicalBorderSurface("TPB Coat", lArFiducial, tpbCoating, tpbCoat);
    new G4LogicalBorderSurface("TPB2 Coat", tpbCoating, lArFiducial,tpbCoat);
    
    if (fLayers) {
      if (pmt_z < -51.0 && radius > 5625.0) {
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial5,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial5,
				   tpbCoat);
      }      
      
      if (pmt_z > 51.0 && radius > 4160.0) {
	new G4LogicalBorderSurface("TPB Coat", 
				   lArFiducial1,
				   tpbCoating,
				   tpbCoat);
	new G4LogicalBorderSurface("TPB2 Coat", 
				   tpbCoating,
				   lArFiducial1,
				   tpbCoat);
      }      
    }
    
  } else {
    new G4PVPlacement(0,
		      trans,
		      pmtLogic,
		      namePMT,
		      fLogicFiduc,
		      false,
		      0,
		      true);
    
    G4VPhysicalVolume*  frillPhys = new G4PVPlacement(0,
						      transf,
						      frillLogic,
						      frillName,
						      fLogicFiduc,
						      false,
						      0,
						      false);
		      
    G4VPhysicalVolume*  frilPhys = new G4PVPlacement(0,
						     transf,
						     frilLogic,
						     frilName,
						     fLogicFiduc,
						     false,
						     0,
						     false);
    
    new G4LogicalBorderSurface("Frill Coat", frillPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", lArFiducial, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Frill Coat", frillPhys, fPhysTPBb, frillCoat);
    new G4LogicalBorderSurface("Frill2 Coat", fPhysTPBb, frillPhys, frillCoat);
    new G4LogicalBorderSurface("Fril Coat", frilPhys, lArFiducial, frillCoat);
    new G4LogicalBorderSurface("Fril2 Coat", lArFiducial, frilPhys, frillCoat);
		      
  }//*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void detectorConstruction::SetDefaults() {
  //Resets to default values
  fOuterRadius_pmt = 10.2*cm;

  //G4cout << "Edwardnote: setting defaults" << G4endl;
  
  fTPBfoilOn = true;
  fPMTsOn = true;
  fReflectorOn = true;
  mainVolOn = true;
  fRodin = false;
  cylinderOn = false;
  fLaser=false;
  fSodium=false;
  fAr39=false;
  fCosmic=false;
  fInvBeta=false;
  darkMatter=false;
  alp=false;
  fLayers=false;//true;//better results when false according to simulation set. 
  ccm200 = false;
  conewide = 1.0;//7.555;
  conehigh = 1.0;//0.590;
  ultra = 37.55;
  threehun = 1310.0;
  fifth = 12.2318;
  r5radius = 31.948;
  topthick = 2.12;
  randwide = 0.80;//2.922;
  tpbEff = 0.96274;
  foilEff = .45548;
  tpbAbs = 0.8735;
  randomized = false;
  rootset = false;
  variableString = "Random: All Mean Values";
  rootfile = "defaultfile.root";

}

//method to randomize all variables without correlations
void detectorConstruction::SetRandoms() {
  /*
  //CCM120 Variable Set
  conewide = G4RandFlat::shoot(3.5,9.0);
  conehigh = G4RandFlat::shoot(0.45,0.65);
  ultra = G4RandFlat::shoot(20.0,80.0);
  threehun = G4RandFlat::shoot(1000.0,1600.);
  randwide = G4RandFlat::shoot(0.5,2.5);
  topthick = G4RandFlat::shoot(1.0,50.0);
  fifth = G4RandFlat::shoot(10.0,50.0);
  r5radius = G4RandFlat::shoot(20.0,60.0);
  foilEff = G4RandFlat::shoot(0.30,0.60);*/

  /*
    CCM200 Simulation Notes as of 0608:
    LAr Variables: Abs length = 46.5 - 100 - 1600 - 2400
        Rindex = 1.358 + steepness = 1
	Rayleigh = 100.58 @ 128 + steepness = 1

    PMT Variables: Height = 5.5-6.0
        Frill: Absl = 0.5, refl = 0.17, wide = 0.6, high = 0.27
	RodShift = -1.7 or -2.0
	foilRefl = 0.9 vis, 0.1 uv

    TPB Variables: 
        foil: Eff = 0.76, scatter = <1.0, thick = 0.7 => Set scatter to 0.3, thick to 0.7
	PMT: Eff = 0.95, scatter = 2.5-3.0, thick = 3.5
	Other: Abs = 0.90, rin = 1.67

   */
  
  //ccm200 variable set:
  
  //ccm200 variable set

  //LiquidArgon Variables:
  G4double base = G4RandGauss::shoot(38.5,2.29);//AbsLAr < 200 nm
  while (base < 10) {  base = G4RandGauss::shoot(38.5,2.29); }
  ultra = 100.0;//G4RandFlat::shoot(50.0,200.0);//67.7;//G4RandFlat::shoot(20.0,200.0);//AbsLAr 200-300 nm
  threehun = 1400.;//G4RandFlat::shoot(1000.0,1800.);//AbsLAr 300-400 nm
  G4double mult = 2400.0;//G4RandFlat::shoot(500,3000.0);//AbsLAr 400+ nm
  G4double larRin = 1.69;//G4RandGauss::shoot(1.69,0.3);//1.69;//LAr rindex @ 128
  G4double scaleRin = 1.0;//G4RandFlat::shoot(0.5,2.0);//1.24;//scaling slope of rindex
  G4double ray128 = G4RandGauss::shoot(99.49,6.0);//100.7813004;//LAr rayleigh length @128
  while (ray128 < 10) {  ray128 = G4RandGauss::shoot(94.49,11.0); }
  G4double scaleRay = G4RandGauss::shoot(1.90,0.2);//1.31;//scaling slope of rayleigh
  while (scaleRay < 1.0) {  scaleRay = G4RandGauss::shoot(1.77,0.2); }

  DefineLAr(base,ultra,threehun,mult,larRin,scaleRin,ray128,scaleRay);
  
  //PMT Variables:
  fifth = -0.92;//G4RandGauss::shoot(-0.92,0.25);//vertical shift of the entire rod (missed the cent)
  topthick = G4RandGauss::shoot(7.23,0.939);//5.56;//pmtHeight
  while (topthick < 1.1) { topthick = G4RandGauss::shoot(5.83,0.5); }
  G4double plasAbs = 2.4;//G4RandGauss::shoot(2.25,0.25);//2.25;//
  G4double plasRefl = 0.27;//G4RandFlat::shoot(0.01,0.5);
  G4double plasRin = 1.6;//G4RandGauss::shoot(1.4,0.2);
  //conewide = G4RandFlat::shoot(0.2,4.0);//frill Thick
  //conehigh = G4RandFlat::shoot(0.2,4.0);//frill Height
  
  DefinePlastic(plasAbs,plasRefl,plasRin);//frill absorption length, reflection %, rindex
 
  //TPB Variables:
  tpbEff = G4RandGauss::shoot(0.92,0.11);//0.95;//pmt tpb Eff
  while (tpbEff < 0.6 || tpbEff > 0.999) {  tpbEff = G4RandGauss::shoot(0.90,0.04); }
  foilEff = 0.8;//G4RandGauss::shoot(0.894,0.21);//0.80;//foil tpb eff
  tpbAbs = G4RandGauss::shoot(0.72,0.15);//all tpb abs (visible light)
  while (tpbAbs < 0.5 || tpbAbs > 0.99) {  tpbAbs = G4RandGauss::shoot(0.83,0.8); }
  randwide = 1.0;//G4RandFlat::shoot(0.3,3.0);//1.0;//probability of scattering in foil TPB
  G4double scatter = G4RandGauss::shoot(1.03,0.3);//probability of scattering in pmt TPB
  while (scatter < 0.5) {  scatter = G4RandGauss::shoot(2.0,0.4); }
  conewide = G4RandGauss::shoot(0.93,0.33);//1.0;//foilThick
  while (conewide < 0.1 || conewide > 3.0) { conewide = G4RandGauss::shoot(1.0,0.5); }
  conehigh = G4RandGauss::shoot(0.65,0.33);//pmtThick
  while (conehigh < 0.1 || conehigh > 2.0) { conehigh = G4RandGauss::shoot(0.25,0.2); }
  G4double tpbRin = G4RandGauss::shoot(1.71,0.3);//tpb index of refraction
  while (tpbRin < 1.3 || tpbRin > 2.5) { tpbRin = G4RandGauss::shoot(1.71,0.5); }
  
  DefineTpb(foilEff, tpbEff, tpbAbs, tpbRin, randwide, scatter);
  
  G4RunManager::GetRunManager()->ReinitializeGeometry();

  std::ostringstream oss;

  //RandomSet07/20s: oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t ray128 \t" << ray128 << "\t rodShift \t" << fifth << "\t pmtHeight \t" << topthick << "\t tpbAbs \t" << tpbAbs << "\t plastAbs \t" << plasAbs << "\t pmtThick \t" << conehigh << "\t tpbRin \t" << tpbRin << "\t pmtScatter \t" << scatter << "\t abs400s \t" << mult << "\n";
  //oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t ray128 \t" << ray128 << "\t rodShift \t" << fifth << "\t pmtHeight \t" << topthick << "\t tpbAbs \t" << tpbAbs << "\t foilThick \t" << conewide << "\t pmtThick \t" << conehigh << "\t tpbRin \t" << tpbRin << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\n";

  //oss << "Randoms: cone\t" << conewide << "\t high\t" << conehigh << "\t ultra\t" << ultra << "\t threehun\t" << threehun << "\t unsmooth\t" << randwide << "\t top\t" << topthick << "\t fifth\t" << fifth << "\t rad\t" << r5radius << "\t foil\t" << foilEff << "\t VUVabsorb\t" << base;
  //oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t ray128 \t" << ray128 << "\t rodShift \t" << fifth << "\t pmtHeight \t" << topthick << "\t tpbAbs \t" << tpbAbs << "\t foilThick \t" << conewide << "\t pmtThick \t" << conehigh << "\t pmtEff \t" << tpbEff << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\n";
  //oss << "Randoms_" << rootfile << "\t pmtEff \t" << tpbEff << "\t foilEff \t" << foilEff << "\t tpbAbs \t" << tpbAbs << "\t tpbRin \t" << tpbRin << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\t pmtThick \t" << conehigh << "\t foilThick \t" << conewide << "\t abs100s \t" << base << "\t rodShift \t" << fifth << "\n";
  oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t scaleRay \t" << scaleRay << "\t ray128 \t" << ray128 << "\t rodShift \t" << fifth << "\t pmtEff \t" << tpbEff << "\t tpbAbs \t" << tpbAbs << "\t tpbRin \t" << tpbRin << "\t pmtThick \t" << conehigh << "\t plasAbs \t" << plasAbs << "\t abs200s \t" << ultra << "\t abs300s \t" << threehun << "\t abs400s \t" << mult << "\t scaleRin \t" << scaleRin << "\t larRin \t" << larRin << "\t pmtHeight \t" << topthick  << "\t foilEff \t" << foilEff << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\t foilThick \t" << conewide << "\t plasRin \t" << plasRin << "\t plasRefl \t" << plasRefl << "\n";

  randomized = true;
  
  variableString = oss.str();
  G4cout << variableString << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......                                          
//methods to set the PMT radius as different (does not currently work) and the rodheight for the laser
void detectorConstruction::SetPMTRadius(G4double outerRadius_pmt) {
  fOuterRadius_pmt=outerRadius_pmt;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void detectorConstruction::SetLaserRod(G4double rHeight) {
  rodHeight=rHeight;
  fLaser=true;
  fRodin=false;//turns the calibration rod off. Running both of these at the same time works, but really shouldn't.
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//turns on the calibration rod. Also turns the laser rod off, helpfully.
void detectorConstruction::SetRodin(G4double rod) {
  fRodin=true;
  rodHeight=rod;
  fLaser=false;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//turn PMTs on or off. turning them off results in no output, though, so be warned.
void detectorConstruction::SetPMTsOn(G4bool b) {
  fPMTsOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//turn Sodium source on or off.
//one should always turn on the calibraiton rod as well (see SetRodIn
void detectorConstruction::SetSodiumOn(G4bool b) {
  fSodium=b;
  G4cout << "Set Sodium to " << fSodium << G4endl;
  //G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//turn Argon39 source on or off.
void detectorConstruction::SetAr39(G4bool b) {
  fAr39=b;
  G4cout << "Set Ar39 to " << fAr39 << G4endl;
}
//turn Cosmic source on or off.
void detectorConstruction::SetCosmic(G4bool b) {
  fCosmic=b;
  G4cout << "Set Cosmic to " << fCosmic << G4endl;
}
//MARK FOR REACTOR: set inverse beta event
void detectorConstruction::SetInverseBeta(G4bool b) {
  fInvBeta=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  G4cout << "Set InverseBeta to " << fInvBeta << G4endl;
}
//sets the dark matter and changes the scintillation to nuclear recoil or electronic
void detectorConstruction::SetDarkMatter(G4bool b) {
  darkMatter=b;
  if (b) {
      lAr_mt->AddConstProperty("YIELDRATIO",1.0);//all in singlet for nuclear recoil
      lAr2_mt->AddConstProperty("YIELDRATIO",1.0);//all in singlet for nuclear recoil
  } else {
    lAr_mt->AddConstProperty("YIELDRATIO",0.25);//twentyfive percent in singlet for ee
    lAr2_mt->AddConstProperty("YIELDRATIO",0.25);//twentyfive percent in singlet for ee
  }
  G4cout << "Setting Dark Matter to " << darkMatter << G4endl;  
}
//sets the alp and changes the scintillation to nuclear recoil or electronic
void detectorConstruction::SetALP(G4bool b) {
  alp=b;
  if (b) {
      lAr_mt->AddConstProperty("YIELDRATIO",0.25);//twentyfive percent in singlet for ee
      lAr2_mt->AddConstProperty("YIELDRATIO",0.25);//twentyfive percent in singlet for ee
  } else {
    lAr_mt->AddConstProperty("YIELDRATIO",1.0);//all in singlet for nuclear recoil
    lAr2_mt->AddConstProperty("YIELDRATIO",1.0);//all in singlet for nuclear recoil
  }
  G4cout << "Setting ALP to " << alp << G4endl;  
}
//sets ccm200 on or off. off defaults to ccm120
void detectorConstruction::SetCCM200(G4bool b) {
  ccm200=b;
  if (ccm200) {
    //LAR variables
    G4double base = 39.56;//AbsLAr < 200 nm
    ultra = 100.0;//AbsLAr 200-300 nm
    threehun = 1400.;//AbsLAr 300-400 nm
    G4double mult = 2400.;//AbsLAr 400+ nm
    G4double larRin = 1.69;//LAr rindex @ 128
    G4double scaleRin = 1.0;//scaling slope of rindex
    G4double ray128 = 94.096;//LAr rayleigh length @128
    G4double scaleRay = 1.73;//scaling slope of rayleigh

    //PMT Variables:
    fifth = -0.92;//vertical shift of the entire rod (missed the center)
    topthick = 7.052; //pmtHeight
    G4double plasAbs = 2.4;//plastic absorbtion length cm
    G4double plasRefl = 0.27;//plastic absorbtion length cm
    G4double plasRin = 1.6;//plastic absorbtion length cm
  
    //TPB Variables:
    tpbEff = 0.878;//pmt tpb Eff
    foilEff = 0.80;//foil tpb eff
    tpbAbs = 0.7215;//all tpb abs (visible light)
    randwide = 1.0;//scattering length in foil TPB
    G4double scatter = 1.13;//scattering length in pmt TPB
    conewide = 0.897;//foilThick
    conehigh = 0.601;//pmtThick
    G4double tpbRin = 1.592;//tpb index of refraction

    DefineLAr(base,ultra,threehun,mult,larRin,scaleRin,ray128,scaleRay);
    DefineTpb(foilEff, tpbEff, tpbAbs, tpbRin, randwide, scatter);
    DefinePlastic(plasAbs,plasRefl,plasRin);
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();

}
//turn the tpb on the foils on or off. The tpb on the pmts will remain even if this is false
void detectorConstruction::SetTPBfoilOn(G4bool b) {
  fTPBfoilOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//turns the reflective foils on or off.
void detectorConstruction::SetReflectorOn(G4bool b) {
  fReflectorOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//method to set the rootfilename
void detectorConstruction::SetRootFile (G4String b) {
  rootfile = b;
  rootset = false;
}

//method to modulate a single random variable a defined sigma away from mean
int detectorConstruction::ModulateRandom(G4int var, G4double sigma) {
  G4double variables[] = {7.555, 0.590, 37.55, 1310.0, 2.922, 26.12, 
			  12.2318, 31.948, .45548, 55.9506};
  G4double errors[] = {1.488, 0.075, 18.17, 172.0, 14.17, 0.480,
                       5.92, 11.08, .0797, 6.923};
  
  if (var < 10) {
    variables[var] = variables[var]+errors[var]*sigma;
  } else {
    G4cout << "invalid variable changed!" << G4endl; 
    return 0;
  }   
  conewide = variables[0];
  conehigh = variables[1];
  ultra = variables[2];
  threehun = variables[3];
  randwide = variables[4];
  topthick = variables[5];
  
  fifth = variables[6];
  r5radius = variables[7];
  foilEff = variables[8];
  G4double base = variables[9];

  G4double mult = 2800.0;

  DefineLAr(base,100.0,threehun,mult,1.358,1,100.7813004,1);//ultra is second
  //DefineLAr(base,ultra,fifth,threehun,mult);
  DefineTpb(foilEff, tpbEff, tpbAbs, 1.7, 2.75, 2.75);

  G4RunManager::GetRunManager()->ReinitializeGeometry();

  std::ostringstream oss;
  oss << "Randoms: cone\t" << conewide << "\t high\t" << conehigh << "\t ultra\t" << ultra << "\t threehun\t" << threehun << "\t unsmooth\t" << randwide << "\t top\t" << topthick << "\t fifth\t" << fifth << "\t rad\t" << r5radius << "\t foil\t" << foilEff << "\t VUVabsorb\t" << base;  

  randomized = true;
  
  variableString = oss.str();
  G4cout << variableString << G4endl;
  
  return 1;
}

//throws a set of randomized correlated OM parameters
void detectorConstruction::CorrelateRandom() {
  G4double variables[] = {7.555, 0.590, 37.55, 1310.0, 2.922, 26.12, 12.2318, 31.948, .45548, 55.9506};
  if (ccm200) {
    variables[0] = 39.556;
    variables[1] = 1.73;
    variables[2] = 94.096;
    variables[3] = 0.878;
    variables[4] = 0.7215;
    variables[5] = 1.592;
    variables[6] = 0.601;
    variables[7] = 7.052;
    variables[8] = 1.13;
    variables[9] = 0.537;
  }
  
  G4double errors[10];
  G4double throws[10];

  for (int i = 0; i < 10; ++i){
    throws[i] = G4RandGauss::shoot(0.0,1.0);
  }
  
  //3-18 laser best fit
  errors[0] = throws[0]*1.435616295;
  errors[1] = throws[0]*-0.002265963+throws[1]*0.07467269;
  errors[2] = throws[0]*6.351105428+throws[1]*-4.33286749+throws[2]*15.73753036;
  errors[3] = throws[0]*18.56168772+throws[1]*3.123084468+throws[2]*-8.391898059+throws[3]*169.9346988;
  errors[4] = throws[0]*0.126009605+throws[1]*-0.017048481+throws[2]*0.086180583+throws[3]*-0.050823025+throws[4]*0.428396012;
  errors[5] = throws[0]*-1.737026377+throws[1]*-3.235993852+throws[2]*2.61676093+throws[3]*0.712776544+throws[4]*0.88608069+throws[5]*13.13058164;
  /*3-14 sodium best fit, replaced with 3-21
  errors[6] = throws[6]*7.002392289;
  errors[7] = throws[6]*10.8252071+throws[7]*6.300770135;
  errors[8] = throws[6]*-0.019606209+throws[7]*0.024887974+throws[8]*0.025921294;
  errors[9] = throws[6]*0.290075114+throws[7]*0.41847509+throws[8]*-0.784878399+throws[9]*2.110057621;//*/
  //3-21 sodium best fit
  errors[6] = throws[6]*6.568616069;
  errors[7] = throws[6]*8.262090646+throws[7]*7.799671268;
  errors[8] = throws[6]*-0.035238137+throws[7]*0.034674956+throws[8]*0.051884376;
  errors[9] = throws[6]*-0.645903006+throws[7]*0.551776675+throws[8]*-1.960392196+throws[9]*6.116201031;

  if (ccm200) {
    G4double cov[10][10] =
      {{3.8959, 0.0108, -2.3116, -0.0234, -0.0600, 0.0294, -0.0613, -0.1600, -0.0180, 0.0391 },
       { 0.0108, 0.0433, 0.0468, 0.0044, 0.0000, 0.0032, -0.0016, 0.0112, 0.0078, -0.0015 },
       { -2.3116, 0.0468, 34.1389, 0.0365, 0.0710, 0.1831, -0.0773, 0.2002, -0.1160, -0.0886 },
       { -0.0234, 0.0044, 0.0365, 0.0074, -0.0031, 0.0050, -0.0024, -0.0031, 0.0056, -0.0028 },
       { -0.0600, 0.0000, 0.0710, -0.0031, 0.0095, 0.0048, -0.0005, -0.0020, 0.0004, 0.0029 },
       { 0.0294, 0.0032, 0.1831, 0.0050, 0.0048, 0.0414, -0.0155, -0.0622, 0.0947, -0.0333 },
       { -0.0613, -0.0016, -0.0773, -0.0024, -0.0005, -0.0155, 0.0224, 0.0630, -0.0458, 0.0254 },
       { -0.1600, 0.0112, 0.2002, -0.0031, -0.0020, -0.0622, 0.0630, 1.2100, -0.4828, 0.1828 },
       { -0.0180, 0.0078, -0.1160, 0.0056, 0.0004, 0.0947, -0.0458, -0.4828, 0.5947, -0.1561 },
       { 0.0391, -0.0015, -0.0886, -0.0028, 0.0029, -0.0333, 0.0254, 0.1828, -0.1561, 0.0989 }};
    for (int i = 0; i < 10; ++i) {
      errors[i] = 0.0;
      for (int j = i; j < 10; ++j) {
	errors[i] += throws[j]*cov[i][j];
      }
    } 
  }

  for (int i = 0; i < 10; ++i){
    variables[i] = variables[i] + errors[i];
    if (variables[i] < 0) {
      variables[i] = -1*variables[i];
    }
  }
  
  std::ostringstream oss;
  if (ccm200) {
    //LAR variables
    G4double base = variables[0];//AbsLAr < 200 nm
    ultra = 100.0;//AbsLAr 200-300 nm
    threehun = 1400.;//AbsLAr 300-400 nm
    G4double mult = 2400.;//AbsLAr 400+ nm
    G4double larRin = 1.69;//LAr rindex @ 128
    G4double scaleRin = 1.0;//scaling slope of rindex
    G4double ray128 = variables[2];//LAr rayleigh length @128
    G4double scaleRay = variables[1];//scaling slope of rayleigh

    //PMT Variables:
    fifth = -0.92;//vertical shift of the entire rod (missed the center)
    topthick = variables[7]; //pmtHeight
    G4double plasAbs = 2.4;//plastic absorbtion length cm
    G4double plasRefl = 0.27;//plastic absorbtion length cm
    G4double plasRin = 1.6;//plastic absorbtion length cm
  
    //TPB Variables:
    tpbEff = variables[3];//pmt tpb Eff
    while (tpbEff > 0.999) { tpbEff = 0.999; }
    foilEff = 0.80;//foil tpb eff
    tpbAbs = variables[4];//all tpb abs (visible light)
    while (tpbAbs > 0.999) { tpbAbs = 0.999; }
    randwide = 1.0;//scattering length in foil TPB
    G4double scatter = variables[8];//scattering length in pmt TPB
    conewide = variables[9];//foilThick
    conehigh = variables[6];//pmtThick
    G4double tpbRin = variables[5];//tpb index of refraction
    while (tpbRin < 1.001) { tpbRin = 1.001; }

    DefineLAr(base,ultra,threehun,mult,larRin,scaleRin,ray128,scaleRay);
    DefineTpb(foilEff, tpbEff, tpbAbs, tpbRin, randwide, scatter);
    DefinePlastic(plasAbs,plasRefl,plasRin);

    oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t scaleRay \t" << scaleRay << "\t ray128 \t" << ray128 << "\t rodShift \t" << fifth << "\t pmtEff \t" << tpbEff << "\t tpbAbs \t" << tpbAbs << "\t tpbRin \t" << tpbRin << "\t pmtThick \t" << conehigh << "\t plasAbs \t" << plasAbs << "\t abs200s \t" << ultra << "\t abs300s \t" << threehun << "\t abs400s \t" << mult << "\t scaleRin \t" << scaleRin << "\t larRin \t" << larRin << "\t pmtHeight \t" << topthick  << "\t foilEff \t" << foilEff << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\t foilThick \t" << conewide << "\t plasRin \t" << plasRin << "\t plasRefl \t" << plasRefl << "\n";
  } else {
    conewide = variables[0];
    conehigh = variables[1];
    ultra = variables[2];
    threehun = variables[3];
    randwide = variables[4];
    topthick = variables[5];
    
    fifth = variables[6];
    r5radius = variables[7];
    foilEff = variables[8];
    G4double base = variables[9];
    
    G4double mult = 2800.0;
    
    DefineLAr(base,100.0,threehun,mult,1.358,1,100.7813004,1);//ultra is second
    //DefineLAr(base,ultra,fifth,threehun,mult);
    DefineTpb(foilEff, tpbEff, tpbAbs, 1.7, 2.75, 2.75);
    
    oss << "Randoms: cone\t" << conewide << "\t high\t" << conehigh << "\t ultra\t" << ultra << "\t threehun\t" << threehun << "\t unsmooth\t" << randwide << "\t top\t" << topthick << "\t fifth\t" << fifth << "\t rad\t" << r5radius << "\t foil\t" << foilEff << "\t VUVabsorb\t" << base;  
  }
  G4RunManager::GetRunManager()->ReinitializeGeometry();
    
  randomized = true;
  
  variableString = oss.str();
  G4cout << variableString << G4endl;
  
  //  return 1;
}

//method for randomizing one parameter at a time
void detectorConstruction::OneRandom(G4int var) {
  /*// ccm120 Variable set
  G4double base = 55.9506;

  if (var > 11) {
    return;
  } else if (var == 0) {
    base = 55.9506;
  } else if (var == 1) {
    conewide = G4RandFlat::shoot(3.5,9.0);
  } else if (var == 2) {
    conehigh = G4RandFlat::shoot(0.45,0.65);
  } else if (var == 3) {
    ultra = G4RandFlat::shoot(20.0,80.0);
  } else if (var == 4) {
    threehun = G4RandFlat::shoot(1000.0,1600.);
  } else if (var == 5) {
    randwide = G4RandFlat::shoot(0.5,2.5);
  } else if (var == 6) {
    topthick = G4RandFlat::shoot(1.0,50.0);
  } else if (var == 7) {
    fifth = G4RandFlat::shoot(10.0,50.0);
  } else if (var == 8) {
    r5radius = G4RandFlat::shoot(20.0,60.0);
  } else if (var == 9) {
    foilEff = G4RandFlat::shoot(0.30,0.60);
  } else if (var == 10) {
    base = G4RandFlat::shoot(20.0,70.0);
    }*/

  //ccm200 variable set
  //LAR variables
  G4double base = 39.25;//AbsLAr < 200 nm
  ultra = 100.0;//AbsLAr 200-300 nm
  threehun = 1400.;//AbsLAr 300-400 nm
  G4double mult = 2400.;//AbsLAr 400+ nm
  G4double larRin = 1.69;//LAr rindex @ 128
  G4double scaleRin = 1.0;//scaling slope of rindex
  G4double ray128 = 94.49;//100.78;//LAr rayleigh length @128
  G4double scaleRay = 1.77;//scaling slope of rayleigh

  //PMT Variables:
  fifth = -0.92;//vertical shift of the entire rod (missed the center)
  topthick = 5.83; //pmtHeight
  //conewide = G4RandFlat::shoot(0.2,4.0);//frill Thick
  //conehigh = G4RandFlat::shoot(0.2,4.0);//frill Height
  G4double plasAbs = 2.4;//plastic absorbtion length cm
  G4double plasRefl = 0.27;//plastic absorbtion length cm
  G4double plasRin = 1.6;//plastic absorbtion length cm
  
  //TPB Variables:
  tpbEff = 0.96;//G4RandGauss::shoot(0.9,0.11);//pmt tpb Eff
  foilEff = 0.80;//G4RandGauss::shoot(0.894,0.21);//foil tpb eff
  tpbAbs = 0.83;//G4RandGauss::shoot(0.913,0.172);//all tpb abs (visible light)
  randwide = 1.0;//G4RandFlat::shoot(0.2,5.0);//probability of scattering in foil TPB
  G4double scatter = 2.0;//G4RandFlat::shoot(0.2,5.0);//probability of scattering in pmt TPB
  conewide = 1.0;//G4RandFlat::shoot(0.2,4.0);//foilThick
  conehigh = 0.30;//G4RandFlat::shoot(0.2,4.0);//pmtThick
  G4double tpbRin = 1.62;//1.38;//G4RandFlat::shoot(1.1,2.0);//tpb index of refraction
  if (var > 11) {
    base = 42.5;
    //return;
  } else if (var == 0) {
    base = G4RandFlat::shoot(30.0,50.0);
  } else if (var == 1) {
    fifth = G4RandFlat::shoot(-5.0,2.0);//vertical shift of the entire rod
  } else if (var == 2) {
    topthick = G4RandFlat::shoot(4.0,8.0);//AbsLAr 200-300 nm
  } else if (var == 3) {
    scaleRay = G4RandFlat::shoot(1.0,2.5);//AbsLAr 300-400 nm
  } else if (var == 4) {
    tpbRin = G4RandFlat::shoot(1.2,2.4);//AbsLAr 400+ nm
  } else if (var == 5) {
    tpbAbs = G4RandFlat::shoot(0.60,0.999);//LAr rindex @ 128
  } else if (var == 6) {
    plasAbs = G4RandFlat::shoot(1.5,3.5);
    //randwide = G4RandFlat::shoot(0.2,3.0);//FoilScattering
  } else if (var == 7) {
    ray128 = G4RandFlat::shoot(70.0,120.0);//100.7813004;//LAr rayleigh length @128
  } else if (var == 8) {
    tpbEff = G4RandFlat::shoot(0.6,0.999);//scaling slope of rayleigh
  } else if (var == 9) {
    conehigh = G4RandFlat::shoot(0.15,1.0);//pmtHeight
  } else if (var == 10) {
    fifth = G4RandFlat::shoot(-6.0,3.0);//vertical shift of the entire rod
  }
  DefineLAr(base,ultra,threehun,mult,larRin,scaleRin,ray128,scaleRay);
  //DefineLAr(base,100.0,threehun,mult,1.358,1,100.7813004,1);//ultra is second
  DefineTpb(foilEff, tpbEff, tpbAbs, tpbRin, randwide, scatter);
  DefinePlastic(plasAbs,plasRefl,plasRin);//absorption length, reflection %, rindex
  
  G4RunManager::GetRunManager()->ReinitializeGeometry();

  std::ostringstream oss;

  oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t larRin \t" << larRin << "\t scaleRay \t" << scaleRay << "\t rodShift \t" << fifth << "\t pmtEff \t" << tpbEff << "\t tpbAbs \t" << tpbAbs << "\t tpbRin \t" << tpbRin << "\t pmtHeight \t" << topthick  << "\t plasRin \t" << plasRin << "\t plasRefl \t" << plasRefl << "\t abs200s \t" << ultra << "\t abs300s \t" << threehun << "\t abs400s \t" << mult << "\t scaleRin \t" << scaleRin << "\t ray128 \t" << ray128 << "\t pmtThick \t" << conehigh << "\t foilEff \t" << foilEff << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\t foilThick \t" << conewide << "\t plasAbs \t" << plasAbs << "\n";

  //oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t ray128 \t" << ray128 << "\t rodShift \t" << fifth << "\t pmtHeight \t" << topthick << "\t tpbAbs \t" << tpbAbs << "\t plastAbs \t" << plastAbs << "\t pmtThick \t" << conehigh << "\t tpbRin \t" << tpbRin << "\t pmtScatter \t" << scatter << "\t scaleRay \t" << scaleRay << "\n";
  //oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t abs200s \t" << ultra << "\t abs300s \t" << threehun << "\t abs400s \t" << mult << "\t larRin \t" << larRin << "\t scaleRin \t" << scaleRin << "\t ray128 \t" << ray128 << "\t scaleRay \t" << scaleRay << "\t rodShift \t" << fifth << "\t pmtHeight \t" << topthick << "\n";
  //oss << "Randoms_" << rootfile << "\t pmtEff \t" << tpbEff << "\t foilEff \t" << foilEff << "\t tpbAbs \t" << tpbAbs << "\t tpbRin \t" << tpbRin << "\t pmtScatter \t" << scatter << "\t foilScatter \t" << randwide << "\t pmtThick \t" << conehigh << "\t foilThick \t" << conewide << "\t abs100s \t" << base << "\t larRin \t" << ultra << "\n";
  //oss << "Randoms_" << rootfile << "\t abs100s \t" << base << "\t rodShift \t" << fifth << "\t tpbScatter \t" << randwide << "\t pmtHeight \t" << topthick << "\t tpbAbs \t" << tpbAbs  << "\t frillAbsl \t" << conewide << "\t frillRefl \t" << conehigh << "\t pmtEff \t" << tpbEff << "\t foilEff \t" << foilEff << "\t larRin \t" << ultra << "\n";

  randomized = true;
  
  variableString = oss.str();
  G4cout << variableString << G4endl;
}

//method for reverting to 'clean' argon (low contamination levels).
void detectorConstruction::CleanArgon() {
  ultra = 1000.;
  fifth = 100.0;
  
  DefineLAr(ultra,ultra,2800.,2800.,1.358,1,100.7813004,1);//ultra is second
  //DefineLAr(ultra, ultra, fifth, 2800., 2800.);

  std::ostringstream oss;

    oss << "Randoms: cone\t" << conewide << "\t high\t" << conehigh << "\t ultra\t" << ultra << "\t threehun\t" << threehun << "\t unsmooth\t" << randwide << "\t top\t" << topthick << "\t fifth\t" << fifth << "\t rad\t" << r5radius << "\t foil\t" << foilEff << "\t VUVabsorb\t" << ultra;  
  
  randomized = true;
  
  variableString = oss.str();
  G4cout << variableString << G4endl;

}

