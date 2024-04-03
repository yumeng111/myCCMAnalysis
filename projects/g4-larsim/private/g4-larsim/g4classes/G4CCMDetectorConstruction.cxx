

#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4CCMScintSD.h"
#include "g4-larsim/g4classes/G4CCMMainVolume.h"
#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

#include <sstream>

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Tubs.hh>
#include <globals.hh>
#include <G4Sphere.hh>
#include <G4Material.hh>
#include <G4UImanager.hh>
#include <G4SDManager.hh>
#include <G4SolidStore.hh>
#include <G4RunManager.hh>
#include <G4PVPlacement.hh>
#include <G4ThreeVector.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include <G4SystemOfUnits.hh>
#include <G4LogicalVolume.hh>
#include <G4MaterialTable.hh>
#include <G4OpticalSurface.hh>
#include <G4GeometryManager.hh>
#include <G4SubtractionSolid.hh>
#include <G4PhysicalConstants.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalBorderSurface.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4CCMDetectorConstruction::G4CCMDetectorConstruction()
{
  SetDefaults();
  DefineMaterials();
  fDetectorMessenger = new G4CCMDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMDetectorConstruction::~G4CCMDetectorConstruction()
{
  delete fMainVolume;
  delete fDetectorMessenger;
  delete fH;
  delete fC;
  delete fN;
  delete fO;
  delete elFe;
  delete elCr;
  delete elNi;
  delete elC;
  delete elMn;
  delete elF;
  delete fGlass;
  delete fPlastic;
  delete fBlackPlastic;
  delete fLAr;
  delete fAlum;
  delete fSteel;
  delete fVacuum;
  delete fPTFE;
  delete fTPBCoating;
  delete fTPBFoil;
  delete fGlass_mt;
  delete fPlastic_mt;
  delete fBlackPlastic_mt;
  delete fLAr_mt;
  delete fAlum_mt;
  delete fSteel_mt;
  delete fVacuum_mt;
  delete fPTFE_mt;
  delete fTPBCoating_mt;
  delete fTPBFoil_mt;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMDetectorConstruction::DefineMaterials() {
    G4double a;  // atomic mass
    G4double z;  // atomic number
    G4double density;

    //***Elements
    fH = new G4Element("H", "H", z = 1., a = 1.01 * g / mole);
    fC = new G4Element("C", "C", z = 6., a = 12.01 * g / mole);
    fN = new G4Element("N", "N", z = 7., a = 14.01 * g / mole);
    fO = new G4Element("O", "O", z = 8., a = 16.00 * g / mole);

    //***Materials
    // Liquid Argon
    fLAr = new G4Material("LAr", z=18., a=39.95*g/mole, density = 1.396*g/cm3, kStateLiquid,88*kelvin);
    
    // Aluminum
    fAlum = new G4Material("Alum", z=13., a=26.98*g/mole, density=2.7*g/cm3);
    
    // Vacuum
    fVacuum = new G4Material("Vacuum", z=1., a=1.01*g/mole, density=universe_mean_density, kStateGas, 0.1*kelvin, 1.e-19*pascal);

    // Steel
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Element* elFe = nistManager->FindOrBuildElement("Fe");
    G4Element* elCr = nistManager->FindOrBuildElement("Cr");
    G4Element* elNi = nistManager->FindOrBuildElement("Ni");
    G4Element* elC = nistManager->FindOrBuildElement("C");
    G4Element* elMn = nistManager->FindOrBuildElement("Mn");
    G4Material* fSteel = new G4Material("Steel", 7.8*g/cm3, 5); // 5 elements in stainless steel
    fSteel->AddElement(elFe, 0.70); // 70% Iron
    fSteel->AddElement(elCr, 0.18); // 18% Chromium
    fSteel->AddElement(elNi, 0.10); // 10% Nickel
    fSteel->AddElement(elC,  0.02); // 2% Carbon
    fSteel->AddElement(elMn, 0.005); // 0.5% Manganese

    // PTFE (for reflector foils) 
    G4Element* elF = nistManager->FindOrBuildElement("F");
    G4Material* fPTFE = new G4Material("PTFE", 2.2*g/cm3, 2); // 2 elements in Teflon
    fPTFE->AddElement(elC, 0.240183); // Fraction of Carbon in fPTFE
    fPTFE->AddElement(elF, 0.759817); // Fraction of Fluorine in fPTFE

    // Glass
    fGlass = new G4Material("Glass", density=1.032*g/cm3, 2);
    fGlass->AddElement(fC, 91.533 * perCent);
    fGlass->AddElement(fH, 8.467 * perCent);

    // TPB (one for coating on PMTs and one for foils)
    fTPBCoating = new G4Material("TPBCoating", density= 1.079*g/cm3, 2);
    fTPBCoating->AddElement(fC, 28);
    fTPBCoating->AddElement(fH, 22);
    
    fTPBFoil = new G4Material("TPBFoil", density= 1.079*g/cm3, 2);
    fTPBFoil->AddElement(fC, 28);
    fTPBFoil->AddElement(fH, 22);

    // Plastic for PMT frill
    fPlastic = new G4Material("Plastic", density=1.20*g/cm3,3);
    fPlastic->AddElement(fC,15);
    fPlastic->AddElement(fH,16);
    fPlastic->AddElement(fO,2);
    
    // Black plastic for PMT frill
    fBlackPlastic = new G4Material("BlackPlastic", density=1.20*g/cm3,3);
    fBlackPlastic->AddElement(fC,15);
    fBlackPlastic->AddElement(fH,16);
    fBlackPlastic->AddElement(fO,2);

    //***Material properties tables

    std::vector<G4double> LAr_Energy_Scint = { 3.87*eV , 4.51*eV , 4.74*eV , 5.03*eV , 5.36*eV , 5.55*eV , 5.82*eV , 6.06*eV , 6.54*eV , 6.79*eV ,
                                               7.03*eV , 7.36*eV , 7.76*eV , 7.98*eV , 8.33*eV , 8.71*eV , 8.96*eV , 9.33*eV , 9.91*eV , 10.31*eV ,
                                               10.61*eV , 10.88*eV , 11.27*eV , 11.81*eV , 12.40*eV , 13.05*eV , 13.78*eV , 14.59*eV , 15.50*eV }; // energies for scintillation spectrum

    std::vector<G4double> LAr_Energy_RIn = {1.771210*eV , 2.066412*eV , 2.479694*eV , 3.099618*eV ,
                                            4.132823*eV , 6.199235*eV , 6.888039*eV , 7.749044*eV , 
                                            8.856050*eV , 9.252590*eV , 9.686305*eV , 9.998766*eV , 
                                            10.33206*eV , 11.27134*eV , 12.39847*eV }; // energies for refractive index and Rayleigh scattering lengths

    std::vector<G4double> LAr_SCINT = {0.00006, 0.00007, 0.00008, 0.00011, 0.00020, 0.00030, 0.00048, 0.00082, 0.00126, 0.00084, 0.00043, 0.00030,
                                       0.00106, 0.00298, 0.00175, 0.00351, 0.01493, 0.12485, 0.49332, 0.20644, 0.07477, 0.04496, 0.01804, 0.00576,
                                       0.00184, 0.00059, 0.00019, 0.00006, 0.00002 }; // liquid Argon scintillation spectrum. this one is centered at 128 nm (as it should be). 

    std::vector<G4double> LAr_RIND = {1.22 , 1.222 , 1.225 , 1.23 ,
                                      1.24 , 1.255 , 1.263 , 1.28 , 
                                      1.315, 1.335 , 1.358 , 1.403, 
                                      1.45 , 1.62  , 1.79 }; // index of refraction spectrum.

    std::vector<G4double> LAr_RSL = {327028.6808*cm, 172560.2267*cm, 80456.5339*cm, 31177.44642*cm, 
                                     8854.144327*cm, 1496.876298*cm, 906.5011168*cm, 480.2538294*cm, 
                                     205.3758714*cm, 145.6326111*cm, 100.7813004*cm, 63.2898117*cm, 
                                     40.07450411*cm, 11.43903548*cm, 3.626432195*cm }; // spectrum of rayleigh scattering lengths.

    // LAr absorption length...plz change at some point
    std::vector<G4double> LAr_Energy_Abs = {1.239847*eV, 1.26515*eV, 1.2915072916666666*eV, 1.3189861702127659*eV, 1.3476597826086956*eV, 1.3776077777777778*eV,
                                            1.4089170454545454*eV, 1.441682558139535*eV, 1.4760083333333334*eV, 1.512008536585366*eV, 1.54980875*eV, 1.589547435897436*eV,
                                            1.6313776315789474*eV, 1.675468918918919*eV, 1.7220097222222222*eV, 1.77121*eV, 1.823304411764706*eV, 1.8785560606060605*eV,
                                            1.9372609375*eV, 1.9997532258064517*eV, 2.0664116666666668*eV, 2.1376672413793103*eV, 2.2140125*eV, 2.296012962962963*eV,
                                            2.384321153846154*eV, 2.479694*eV, 2.583014583333333*eV, 2.695319565217391*eV, 2.817834090909091*eV, 2.9520166666666667*eV,
                                            3.0996175*eV, 3.262755263157895*eV, 3.4440194444444443*eV, 3.646608823529412*eV, 3.874521875*eV, 4.1328233333333335*eV, 4.428025*eV,
                                            4.768642307692308*eV, 5.166029166666666*eV, 5.635668181818182*eV, 6.199235*eV, 6.888038888888889*eV, 7.74904375*eV, 8.85605*eV,
                                            10.332058333333332*eV, 12.39847*eV, 15.4980875*eV};
        
    std::vector<G4double> LAr_ABS = {2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm,
                                     2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm,
                                     2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 2800.0*cm, 1310.0*cm, 1310.0*cm, 1310.0*cm, 1310.0*cm, 1310.0*cm, 100.0*cm, 100.0*cm, 100.0*cm,
                                     100.0*cm, 100.0*cm, 55.9506*cm, 55.9506*cm, 55.9506*cm, 55.9506*cm, 55.9506*cm, 55.9506*cm};

    // Takes the defined values above and uses them to define a materials properties table.
    
    G4MaterialPropertiesTable* fLAr_mt = new G4MaterialPropertiesTable();
    fLAr_mt->AddProperty("SCINTILLATIONCOMPONENT1", LAr_Energy_Scint, LAr_SCINT); 
    fLAr_mt->AddProperty("SCINTILLATIONCOMPONENT2", LAr_Energy_Scint, LAr_SCINT);
    fLAr_mt->AddProperty("RINDEX",        LAr_Energy_RIn,   LAr_RIND);
    fLAr_mt->AddProperty("RAYLEIGH",      LAr_Energy_RIn,  LAr_RSL);
    fLAr_mt->AddProperty("ABSLENGTH", LAr_Energy_Abs, LAr_ABS);
    G4double scint_yeild=1.0/(19.5*eV); // scintillation yield: 50 per keV.
    fLAr_mt->AddConstProperty("SCINTILLATIONYIELD", scint_yeild);
    fLAr_mt->AddConstProperty("RESOLUTIONSCALE",0.11); 
    fLAr_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",7.*ns);
    fLAr_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",1400.*ns);
    fLAr_mt->AddConstProperty("SCINTILLATIONYIELD1",0.25); // for e/m scintillation
    fLAr->SetMaterialPropertiesTable(fLAr_mt);


    // Set the Birks Constant for the LAr scintillator
    fLAr->GetIonisation()->SetBirksConstant(0.0486*mm/MeV);

    // Set PMT glass constants
    std::vector<G4double> glass_energy = {7.0*eV, 7.07*eV, 7.14*eV};
    std::vector<G4double> glass_AbsLength = { 420. * cm, 420. * cm, 420. * cm };
    std::vector<G4double> glass_RIND = {1.49,1.49,1.49};

    G4MaterialPropertiesTable* fGlass_mt = new G4MaterialPropertiesTable();
    fGlass_mt->AddProperty("ABSLENGTH", glass_energy, glass_AbsLength);
    fGlass_mt->AddProperty("RINDEX", glass_energy, glass_RIND);
    fGlass->SetMaterialPropertiesTable(fGlass_mt);

    std::vector<G4double> Vacuum_Energy = {2.0*eV,7.0*eV,7.14*eV};
    std::vector<G4double> Vacuum_RIND={1.,1.,1.};
   
    G4MaterialPropertiesTable* fVacuum_mt = new G4MaterialPropertiesTable(); 
    fVacuum_mt->AddProperty("RINDEX", Vacuum_Energy, Vacuum_RIND);
    fVacuum->SetMaterialPropertiesTable(fVacuum_mt);

    // Definition of MPT for Plastic frills
    std::vector<G4double> plastic_Energy = { 1.0*eV,1.2*eV,2.5*eV,3.0*eV,3.4*eV,6.5*eV,10.0*eV,12.0*eV };
    std::vector<G4double> plastic_reflect = {0.10, 0.10, 0.25, 0.30, 0.10, 0.05, 0.01, 0.01};
    std::vector<G4double> plastic_AbsLength = { 10.0*cm, 10.0*cm, 10.0*cm, 10.0*cm, 1.0e-3*cm, 1.0e-3*cm, 1.0e-3*cm, 1.0e-3*cm};
    std::vector<G4double> blackplastic_AbsLength = { 1e-9*cm, 1e-9*cm, 1e-9*cm, 1e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm, 1.0e-9*cm};
    std::vector<G4double> plastic_RIND = { 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4 };
   
    G4MaterialPropertiesTable* fPlastic_mt = new G4MaterialPropertiesTable(); 
    fPlastic_mt->AddProperty("ABSLENGTH", plastic_Energy, plastic_AbsLength);
    fPlastic_mt->AddProperty("REFLECTIVITY", plastic_Energy, plastic_reflect);
    fPlastic_mt->AddProperty("RINDEX", plastic_Energy, plastic_RIND);
    fPlastic->SetMaterialPropertiesTable(fPlastic_mt);
   
    G4MaterialPropertiesTable* fBlackPlastic_mt = new G4MaterialPropertiesTable();  
    fBlackPlastic_mt->AddProperty("ABSLENGTH", plastic_Energy, blackplastic_AbsLength);
    fBlackPlastic_mt->AddProperty("REFLECTIVITY", plastic_Energy, plastic_reflect);
    fBlackPlastic_mt->AddProperty("RINDEX", plastic_Energy, plastic_RIND);
    fBlackPlastic->SetMaterialPropertiesTable(fBlackPlastic_mt);

    // Aluminum
    std::vector<G4double> alum_energy = { 1.53353916299*eV, 1.55836885306*eV, 1.58716054918*eV, 1.61709864659*eV, 1.65046713174*eV,    
                                          1.6817136839*eV,  1.71786713359*eV, 1.75433728722*eV, 1.79106381297*eV, 1.83220935967*eV, 
                                          1.87670110428*eV, 1.92033477064*eV, 1.96440460941*eV, 2.01566223496*eV, 2.06248774244*eV, 
                                          2.1191459895*eV,  2.17704608625*eV,  2.23198006977*eV, 2.29412813935*eV,2.37159279237*eV, 
                                          2.43945693286*eV, 2.51927088373*eV, 2.59881298342*eV, 2.68963236717*eV, 2.78391786144*eV, 
                                          2.88535232931*eV, 2.97955997826*eV, 3.08041053689*eV, 3.18832742402*eV, 3.26833943848*eV, 
                                          3.36191419687*eV, 3.46107669506*eV, 3.52980225748*eV, 3.60084752681*eV, 3.70923697374*eV, 
                                          3.75238245064*eV, 3.82680283055*eV, 3.9042347515*eV,  3.9111734511*eV, 20.6640330667*eV };
    std::vector<G4double> alum_reflect = {0.8866682, 0.8975349, 0.9091666, 0.9176283, 0.9260717, 
                                          0.9329575, 0.9398158, 0.9458907, 0.9511823, 0.9548705, 
                                          0.9609272, 0.9662096, 0.9722937, 0.9775578, 0.9836511, 
                                          0.9865376, 0.9886408, 0.9891864, 0.9905061, 0.9909877, 
                                          0.9907316, 0.9904481, 0.9901829, 0.9891068, 0.9856623, 
                                          0.9774627, 0.968507,  0.955589,  0.9426708, 0.9242781, 
                                          0.9058671, 0.8866637, 0.862769,  0.8436296, 0.8228502, 
                                          0.8005955, 0.7798801, 0.7591648, 0.001,     0.001 };    

    std::vector<G4double> alum_abslen = { 1.0e-3*cm, 1.0e-3*cm};
    std::vector<G4double> alum_abseneg = { 1.5498*eV, 20.664*eV};

    G4MaterialPropertiesTable* fAlum_mt = new G4MaterialPropertiesTable();
    fAlum_mt->AddProperty("REFLECTIVITY", alum_energy, alum_reflect);
    fAlum_mt->AddProperty("ABSLENGTH", alum_abseneg, alum_abslen);
    fAlum->SetMaterialPropertiesTable(fAlum_mt);

    // TPB
    std::vector<G4double> TPBEnergy = { 0.602*eV/*(2066nm)*/, 0.689*eV/*(1799nm)*/, 1.030*eV/*(1204nm)*/, 1.926*eV/*(644nm)*/, 2.138*eV/* (580nm)*/,
                                        2.250*eV/*(551nm)*/,  2.380*eV/*(521nm)*/,  2.480*eV/*(500nm)*/,  2.583*eV/*(480nm)*/, 2.800*eV/*(443nm)*/,
                                        2.880*eV/*(431nm)*/,  2.980*eV/*(416nm)*/,  3.124*eV/*(397nm)*/,  3.457*eV/*(359nm)*/, 3.643*eV/*(341nm)*/,
                                        3.812*eV/*(325nm)*/,  4.086*eV/*(304nm)*/,  4.511*eV/*(275nm)*/,  5.166*eV/*(240nm)*/, 5.821*eV/*(213nm)*/,
                                        6.526*eV/*(190nm)*/,  8.266*eV/*(150nm)*/,  9.686*eV/*(128nm)*/,  11.27*eV/*(110nm)*/, 12.60*eV/*(98nm)*/  };

    // Emission spectrum for TPB; basically probabilities of photon being reemitted in the various energy bins.
    std::vector<G4double> TPBEmission = { 0.0000, 0.0000, 0.0000, 0.0000, 0.0005,
                                          0.0015, 0.0030, 0.0050, 0.0070, 0.0110,
                                          0.0110, 0.0060, 0.0020, 0.0000, 0.0000,
                                          0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                                          0.0000, 0.0000, 0.0000, 0.0000, 0.0000 };

    // Refractive index of the TPB.
    std::vector<G4double> TPBRIndex = { 1.4, 1.4, 1.4, 1.4, 1.4,
                                        1.4, 1.4, 1.4, 1.4, 1.4,
                                        1.4, 1.4, 1.4, 1.4, 1.4,
                                        1.4, 1.4, 1.4, 1.4, 1.4,
                                        1.4, 1.4, 1.4, 1.4, 1.4 };


    G4MaterialPropertiesTable* fTPBCoating_mt = new G4MaterialPropertiesTable();
    //fTPBCoating_mt->AddProperty("RINDEX", TPBEnergy, TPBRIndex);
    fTPBCoating_mt->AddProperty("WLSCOMPONENT", TPBEnergy, TPBEmission);
    fTPBCoating_mt->AddConstProperty("WLSTIMECONSTANT", 1.7*ns);
    fTPBCoating_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.2);
    fTPBCoating->SetMaterialPropertiesTable(fTPBCoating_mt);
    
    G4MaterialPropertiesTable* fTPBFoil_mt = new G4MaterialPropertiesTable();
    //fTPBFoil_mt->AddProperty("RINDEX", TPBEnergy, TPBRIndex);
    fTPBFoil_mt->AddProperty("WLSCOMPONENT", TPBEnergy, TPBEmission);
    fTPBFoil_mt->AddConstProperty("WLSTIMECONSTANT", 1.7*ns);
    fTPBFoil_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", 1.2);
    fTPBFoil->SetMaterialPropertiesTable(fTPBFoil_mt);

    // Defines properties of the PTFE reflectors.
    std::vector<G4double> TefEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV,
                                       2.138*eV, 2.25*eV,  2.38*eV,  2.48*eV, 
                                       2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV, 
                                       3.124*eV, 3.457*eV, 3.643*eV, 3.812*eV, 4.086*eV,
                                       4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                       7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
    std::vector<G4double> TefReflect = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                                        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.99};
    
    std::vector<G4double> TefRIndex = {10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
                                       10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
                                        10., 10., 10., 10., 10.};
    
    
    G4MaterialPropertiesTable* fPTFE_mt = new G4MaterialPropertiesTable();
    fPTFE_mt->AddProperty("RINDEX", TefEnergy, TefRIndex);
    fPTFE_mt->AddProperty("REFLECTIVITY", TefEnergy, TefReflect);
    fPTFE->SetMaterialPropertiesTable(fPTFE_mt);



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* G4CCMDetectorConstruction::Construct() {

    // Create a box 5x5x5 meters in which to place the detector. In theory is much larger than it needs to be.
    G4double expHall_x = 5.*m;
    G4double expHall_y = 5.*m;
    G4double expHall_z = 5.*m;
    fExperimentalHall_box = new G4Box("expHall",expHall_x,expHall_y,expHall_z);
    fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box,fVacuum,"expHall_log");
    fExperimentalHall_phys = new G4PVPlacement(nullptr, G4ThreeVector(), fExperimentalHall_log, "expHall", nullptr, false, 0);
    
    fExperimentalHall_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Place the main volume
    if(fMainVolumeOn) {
        fMainVolume = new G4CCMMainVolume(nullptr, G4ThreeVector(), fExperimentalHall_log, false, 0, this);
    }

    return fExperimentalHall_phys;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// construct our pmts SD
void G4CCMDetectorConstruction::ConstructSDandField(){
    if(!fMainVolume)
        return;

    if (PMTSDStatus_){
        // PMT SD
        G4CCMPMTSD* pmt = fPMT_SD.Get();
        if(!pmt) {
            // Created here so it exists as pmts are being placed
            G4cout << "Construction /LAr/pmtSD" << G4endl;
            auto pmt_SD = new G4CCMPMTSD("/LAr/pmtSD");
            fPMT_SD.Put(pmt_SD);

            pmt_SD->InitPMTs();
            pmt_SD->SetPmtPositions(fMainVolume->GetPMTPositions());
        }
        else {
            pmt->InitPMTs();
            pmt->SetPmtPositions(fMainVolume->GetPMTPositions());
        }
        G4SDManager::GetSDMpointer()->AddNewDetector(fPMT_SD.Get());
        SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPMT_SD.Get());
    }
    
    if (LArSDStatus_){
        // Scint SD
        if(!fScint_SD.Get()) {
            G4cout << "Construction /LAr/scintSD" << G4endl;
            auto scint_SD = new G4CCMScintSD("/LAr/scintSD");
            scint_SD->SetPMTSDStatus(PMTSDStatus_);
            fScint_SD.Put(scint_SD);
        }
        G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
        SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMDetectorConstruction::SetDefaults() {
    //Resets to default values
    fMainVolumeOn = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMDetectorConstruction::SetMainVolumeOn(G4bool b) {
    fMainVolumeOn = b;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}



