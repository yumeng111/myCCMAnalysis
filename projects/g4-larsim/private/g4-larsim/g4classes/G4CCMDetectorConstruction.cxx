

#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4CCMScintSD.h"
#include "g4-larsim/g4classes/G4CCMMainVolume.h"
#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

#include <sstream>
#include <cmath>

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Tubs.hh>
#include <globals.hh>
#include <G4Sphere.hh>
#include <G4Material.hh>
#include <G4UImanager.hh>
#include <G4SDManager.hh>
#include <G4SolidStore.hh>
#include <G4MTRunManager.hh>
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
G4CCMDetectorConstruction::G4CCMDetectorConstruction(G4double SingletTau, G4double TripletTau, G4double UVAbsLength, G4double WLSNPhotonsEndCapFoil,
                                                     G4double WLSNPhotonsSideFoil, G4double WLSNPhotonsPMT,
                                                     G4double EndCapFoilTPBThickness, G4double SideFoilTPBThickness, G4double PMTTPBThickness,
                                                     G4double Rayleigh128, G4double TPBAbsTau, G4double TPBAbsNorm, G4double TPBAbsScale, G4double Mie_GG, G4double Mie_Ratio) {
    SingletTau_ = SingletTau;
    TripletTau_ = TripletTau;
    UVAbsLength_ = UVAbsLength;
    WLSNPhotonsEndCapFoil_ = WLSNPhotonsEndCapFoil;
    WLSNPhotonsSideFoil_ = WLSNPhotonsSideFoil;
    WLSNPhotonsPMT_ = WLSNPhotonsPMT;
    EndCapFoilTPBThickness_ = EndCapFoilTPBThickness;
    SideFoilTPBThickness_ = SideFoilTPBThickness;
    PMTTPBThickness_ = PMTTPBThickness;
    Rayleigh128_ = Rayleigh128;
    TPBAbsTau_ = TPBAbsTau;
    TPBAbsNorm_ = TPBAbsNorm;
    TPBAbsScale_ = TPBAbsScale;
    Mie_GG_ = Mie_GG;
    Mie_Ratio_ = Mie_Ratio;
    SetDefaults();
    DefineMaterials();
    fDetectorMessenger = new G4CCMDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMDetectorConstruction::~G4CCMDetectorConstruction()
{
  if(fMainVolume != nullptr) delete fMainVolume;
  if(fDetectorMessenger != nullptr) delete fDetectorMessenger;
  if(fH != nullptr) delete fH;
  if(fC != nullptr) delete fC;
  if(fN != nullptr) delete fN;
  if(fO != nullptr) delete fO;
  if(elFe != nullptr) delete elFe;
  if(elCr != nullptr) delete elCr;
  if(elNi != nullptr) delete elNi;
  if(elC != nullptr) delete elC;
  if(elMn != nullptr) delete elMn;
  if(elF != nullptr) delete elF;
  if(fGlass != nullptr) delete fGlass;
  if(fPlastic != nullptr) delete fPlastic;
  if(fBlackPlastic != nullptr) delete fBlackPlastic;
  if(fLAr != nullptr) delete fLAr;
  if(fAlum != nullptr) delete fAlum;
  if(fSteel != nullptr) delete fSteel;
  if(fVacuum != nullptr) delete fVacuum;
  if(fPTFE != nullptr) delete fPTFE;
  if(fTPBFoilSides != nullptr) delete fTPBFoilSides;
  if(fTPBFoilTopBottom != nullptr) delete fTPBFoilTopBottom;
  if(fGlass_mt != nullptr) delete fGlass_mt;
  if(fPlastic_mt != nullptr) delete fPlastic_mt;
  if(fBlackPlastic_mt != nullptr) delete fBlackPlastic_mt;
  if(fLAr_mt != nullptr) delete fLAr_mt;
  if(fAlum_mt != nullptr) delete fAlum_mt;
  if(fSteel_mt != nullptr) delete fSteel_mt;
  if(fVacuum_mt != nullptr) delete fVacuum_mt;
  if(fPTFE_mt != nullptr) delete fPTFE_mt;
  if(fTPBFoilSides_mt != nullptr) delete fTPBFoilSides_mt;
  if(fTPBFoilTopBottom_mt != nullptr) delete fTPBFoilTopBottom_mt;
  if(fTPBPMT_mt != nullptr) delete fTPBPMT_mt;
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
    G4Material* fSteel = new G4Material("Steel", 7.8*g/cm3, 4); // 5 elements in stainless steel
    fSteel->AddElement(elFe, 0.70); // 70% Iron
    fSteel->AddElement(elCr, 0.20); // 20% Chromium
    fSteel->AddElement(elNi, 0.08); // 8% Nickel
    //fSteel->AddElement(elC,  0.008); // 0.08% Carbon
    fSteel->AddElement(elMn, 0.02); // 2% Manganese

    // PTFE (for reflector foils)
    G4Element* elF = nistManager->FindOrBuildElement("F");
    G4Material* fPTFE = new G4Material("PTFE", 2.2*g/cm3, 2); // 2 elements in Teflon
    fPTFE->AddElement(elC, 0.240183); // Fraction of Carbon in fPTFE
    fPTFE->AddElement(elF, 0.759817); // Fraction of Fluorine in fPTFE

    // Glass
    fGlass = new G4Material("Glass", density=1.032*g/cm3, 2);
    fGlass->AddElement(fC, 91.533 * perCent);
    fGlass->AddElement(fH, 8.467 * perCent);

    // TPB Foil
    fTPBFoilSides = new G4Material("TPBFoilSides", density= 1.079*g/cm3, 2);
    fTPBFoilSides->AddElement(fC, 28);
    fTPBFoilSides->AddElement(fH, 22);

    fTPBFoilTopBottom = new G4Material("TPBFoilTopBottom", density= 1.079*g/cm3, 2);
    fTPBFoilTopBottom->AddElement(fC, 28);
    fTPBFoilTopBottom->AddElement(fH, 22);

    // TPB PMT
    fTPBPMT = new G4Material("TPBPMT", density= 1.079*g/cm3, 2);
    fTPBPMT->AddElement(fC, 28);
    fTPBPMT->AddElement(fH, 22);


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
    //fLAr_mt->AddProperty("SCINTILLATIONCOMPONENT2", LAr_Energy_Scint, LAr_SCINT);

    // let's do index of refraction and rayleigh now
    G4double rindex = 1.358;//LAr rindex @ 128
    G4double scaler = 1.0;//scaling slope of rindex
    G4double scalerray = 1.73;//scaling slope of rayleigh
    std::vector<G4double> lar_Energy_rin = { 1.771210*eV , 2.066412*eV , 2.479694*eV , 3.099618*eV , 4.132823*eV , 6.199235*eV , 6.888039*eV , 7.749044*eV ,
                                             8.856050*eV , 9.252590*eV , 9.686305*eV , 9.998766*eV , 10.33206*eV , 11.27134*eV , 12.39847*eV }; //energies for refractive index and Rayleigh scattering lengths

    G4double set = rindex-1.24;
    G4double scaleRin = 1.0;//scaling slope of rindex
    std::vector<G4double> scale = { 0.125, 0.196, 0.345, 0.637, 0.8065, 1.38, 1.78, 3.22, 4.66 };
    for (int i = 0; i < 9; ++i) {
        scale[i] = std::pow(scale[i],scaler)*set+1.24;
    }
    std::vector<G4double> lar_RSL = { 3244.9341, 1712.2246, 798.328,   309.35745, 87.855032, 14.852719, 8.9947353, 4.7653069,
                                      2.0378371, 1.445036,  1,         0.6279916, 0.3976383, 0.1135036, 0.035983 };
    for (int i = 0; i < 15; ++i) {
        lar_RSL[i] = std::pow(lar_RSL[i],scalerray) * Rayleigh128_;
    }

    std::cout << "using " << Rayleigh128_ / cm<< " for Rayleigh scattering length at 128nm" << std::endl;

    const G4int larrin =  sizeof(lar_Energy_rin)/sizeof(G4double);
    std::vector<G4double> lar_RIND  = { 1.22 ,     1.222 ,    1.225 ,    1.23 , 1.24 ,     scale[0] , scale[1] , scale[2] , scale[3] ,
                                        scale[4] , rindex ,   scale[5], scale[6] , scale[7] , scale[8] }; //index of refraction spectrum.

    fLAr_mt->AddProperty("RINDEX", lar_Energy_rin,  lar_RIND, larrin);
    fLAr_mt->AddProperty("RAYLEIGH", lar_Energy_rin,  lar_RSL, larrin);

    // now add absorption length
    std::vector<G4double> flat_abs = {300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm,
                                      300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm,
                                      300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm, 300000*cm,
                                      300000*cm, 300000*cm, 300000*cm, 300000*cm,
                                      UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_,
                                      UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_, UVAbsLength_};

    //fLAr_mt->AddProperty("ABSLENGTH", LAr_Energy_Abs, LAr_ABS);
    //fLAr_mt->AddProperty("ABSLENGTH", LAr_Energy_Abs, flat_abs);

    G4double scint_yeild=1.0/(19.5*eV); // scintillation yield: 50 per keV.
    fLAr_mt->AddConstProperty("SCINTILLATIONYIELD", scint_yeild);
    fLAr_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
    fLAr_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", SingletTau_);
    //fLAr_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT2", TripletTau_);
    //fLAr_mt->AddConstProperty("SCINTILLATIONYIELD1",0.25); // for e/m scintillation
    fLAr_mt->AddConstProperty("SCINTILLATIONYIELD1",1.0); // for e/m scintillation
    fLAr->SetMaterialPropertiesTable(fLAr_mt);

    // Set the Birks Constant for the LAr scintillator
    fLAr->GetIonisation()->SetBirksConstant(0.0486*mm/MeV);

    // Set PMT glass constants
    std::vector<G4double> glass_energy = {7.0*eV, 7.07*eV, 7.14*eV};
    std::vector<G4double> glass_AbsLength = {1e-12*mm, 1e-12*mm, 1e-12*mm};
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
    fPlastic_mt->AddProperty("REFLECTIVITY", plastic_Energy, plastic_reflect); // pretty sure this doesnt do anything
    fPlastic_mt->AddProperty("RINDEX", plastic_Energy, plastic_RIND);
    fPlastic->SetMaterialPropertiesTable(fPlastic_mt);

    G4MaterialPropertiesTable* fBlackPlastic_mt = new G4MaterialPropertiesTable();
    fBlackPlastic_mt->AddProperty("ABSLENGTH", plastic_Energy, blackplastic_AbsLength); // pretty sure this doesnt do anything
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

    std::vector<G4double> alum_abslen = { 1.0e-10*mm, 1.0e-10*mm};
    std::vector<G4double> alum_abseneg = { 1.5498*eV, 20.664*eV};

    G4MaterialPropertiesTable* fAlum_mt = new G4MaterialPropertiesTable();
    //fAlum_mt->AddProperty("REFLECTIVITY", alum_energy, alum_reflect); this doesnt do anything...
    fAlum_mt->AddProperty("ABSLENGTH", alum_abseneg, alum_abslen);
    fAlum->SetMaterialPropertiesTable(fAlum_mt);

    // now time to define foil + pmt TPB
    // emession spectrum for TPB + PTFE (https://arxiv.org/pdf/2211.05024)
    std::vector<G4double> TPB_PTFE_Emission_Energy = {2.0685664749347477*eV, 2.0742413589180417*eV, 2.0805207528485083*eV, 2.0868375149681437*eV, 2.0931932082425835*eV, 2.099587647474017*eV,
                2.106021274818547*eV, 2.112494131525724*eV, 2.119007045952072*eV, 2.12555991971599*eV, 2.132153329185213*eV, 2.138788039429448*eV, 2.145462998960346*eV, 2.152179994420006*eV,
                2.1589388770455353*eV, 2.1657412326952583*eV, 2.1725864662565497*eV, 2.179474860522975*eV, 2.1864067314622337*eV, 2.193383621332569*eV, 2.200405401685091*eV, 2.2074731426719727*eV,
                2.2145859216890855*eV, 2.2217433336672547*eV, 2.228947970252437*eV, 2.236198669829078*eV, 2.243496400573736*eV, 2.250843075164009*eV, 2.258236327514419*eV, 2.265677270871168*eV,
                2.27316704186735*eV, 2.2807053764506375*eV, 2.288294318910927*eV, 2.2959357211202307*eV, 2.3031132208657867*eV, 2.3064249200502367*eV, 2.3089987782740105*eV, 2.310753056791676*eV,
                2.3123937390541376*eV, 2.3137607497117663*eV, 2.315207040598785*eV, 2.315908550533234*eV, 2.3174715527811522*eV, 2.319511938132526*eV, 2.320396040998038*eV, 2.3215564607402435*eV,
                2.3231261582527347*eV, 2.3239333030607687*eV, 2.3252510550513774*eV, 2.326606009528413*eV, 2.3288955123996273*eV, 2.3296218044811994*eV, 2.3315105459000147*eV, 2.3336095670498143*eV,
                2.3392267743396813*eV, 2.3434320787205567*eV, 2.3448964788261755*eV, 2.34672623048496*eV, 2.347830664304569*eV, 2.34948351875669*eV, 2.350956442682144*eV, 2.352616275497997*eV,
                2.3542788612441328*eV, 2.355429742544156*eV, 2.3569122356677346*eV, 2.3583965961108655*eV, 2.359863137579771*eV, 2.3618825640867933*eV, 2.3621648981516272*eV, 2.3643916173329202*eV,
                2.3669920554642006*eV, 2.3715481052664877*eV, 2.37976314520725*eV, 2.3880295328504966*eV, 2.396355234293263*eV, 2.4047384763570903*eV, 2.4131801233476007*eV, 2.4216787608940966*eV,
                2.4302382790609602*eV, 2.438858674901718*eV, 2.4475401705062194*eV, 2.456284361987954*eV, 2.4650891967604536*eV, 2.4739551874273493*eV, 2.4828854247031*eV, 2.4918785039180884*eV,
                2.500937455744144*eV, 2.5092286092893503*eV, 2.515479330062037*eV, 2.5209205592845185*eV, 2.5268086088717574*eV, 2.533148080032263*eV, 2.5399467855740334*eV, 2.5463522450031966*eV,
                2.5540832764731194*eV, 2.5636053856961847*eV, 2.5732054967073124*eV, 2.5828790399538524*eV, 2.592617130429985*eV, 2.6024231434079277*eV, 2.610951434653585*eV, 2.618174289740801*eV,
                2.6249816924571023*eV, 2.6313664255135407*eV, 2.638599130699968*eV, 2.6446894110191903*eV, 2.649376829706567*eV, 2.6562778089959154*eV, 2.662345774917437*eV, 2.6671572148443854*eV,
                2.67209571453112*eV, 2.676393339739669*eV, 2.680419665707901*eV, 2.6853475764163504*eV, 2.6897411247096445*eV, 2.6940419882848023*eV, 2.6978756231800833*eV, 2.7002736900568207*eV,
                2.7075145707266524*eV, 2.7157706136036364*eV, 2.723095384003621*eV, 2.7299663336229294*eV, 2.7368719679359534*eV, 2.743812627384684*eV, 2.751880520718895*eV, 2.7562917843098393*eV,
                2.7628272350469816*eV, 2.7693929790893814*eV, 2.775554342547651*eV, 2.782618920365959*eV, 2.7914414251133404*eV, 2.797005036527513*eV, 2.80177236298375*eV, 2.806323149863457*eV,
                2.8104072536417566*eV, 2.814246190862789*eV, 2.818835698864991*eV, 2.823030428106631*eV, 2.8272374453588065*eV, 2.8314564452626763*eV, 2.83568884883466*eV, 2.83993392437019*eV,
                2.8441149617296557*eV, 2.8473914164715097*eV, 2.851137565797999*eV, 2.855970198617966*eV, 2.86054991169904*eV, 2.867315503133049*eV, 2.876574331652757*eV, 2.8864450385590703*eV,
                2.896942495109509*eV, 2.909198597614656*eV, 2.92044497894063*eV, 2.9275282613965077*eV, 2.9319253376658416*eV, 2.936941004157542*eV, 2.9401399085894635*eV, 2.9424427536811075*eV,
                2.944995247603916*eV, 2.947224026420486*eV, 2.950629931104428*eV, 2.9530961862054896*eV, 2.9551739231653156*eV, 2.957171754757947*eV, 2.9595033311512613*eV, 2.9618384297880627*eV,
                2.9641772949799576*eV, 2.9665198569390037*eV, 2.968866124436694*eV, 2.9712173727466893*eV, 2.9738656346526167*eV, 2.976225923589673*eV, 2.9788203597532155*eV, 2.981323058454508*eV,
                2.9839269958496897*eV, 2.986008450509771*eV, 2.9880922773480423*eV, 2.990476655383541*eV, 2.992863958299538*eV, 2.9952551562778336*eV, 2.997650017146839*eV, 3.000048468614792*eV,
                3.0024511653467503*eV, 3.005089915986259*eV, 3.007310787083497*eV, 3.009681671633881*eV, 3.0121001549438677*eV, 3.0145226911680596*eV, 3.0169492088763286*eV, 3.0193797996297134*eV,
                3.0218142280550877*eV, 3.0242518471992907*eV, 3.02673463613197*eV, 3.02883033610366*eV, 3.031276903173656*eV, 3.0337272608748043*eV, 3.0365176160389673*eV, 3.0402728354689623*eV,
                3.0426348341018628*eV, 3.0463740356874967*eV, 3.048804342523251*eV, 3.0525584735233675*eV, 3.055677528940057*eV, 3.060893002720675*eV, 3.065252695835357*eV, 3.069623671241484*eV,
                3.0758832078926326*eV, 3.085943494297697*eV, 3.095434887250183*eV, 3.1075345808794155*eV, 3.121657771067854*eV, 3.1359095373649297*eV, 3.150290025719692*eV, 3.164802163179325*eV,
                3.17944809538817*eV, 3.1942272194483787*eV, 3.209145791920373*eV, 3.224201524990315*eV, 3.2393993289846046*eV, 3.254738392193289*eV};

    std::vector<G4double> TPB_PTFE_Emission = {0.012512215916579318, 0.013288259609399031, 0.01232653203599243, 0.012821661395969356, 0.012453468128755732, 0.012247147854140541,
                0.012040827579525168, 0.01242804161110332, 0.012545467321684265, 0.013256427338458736, 0.014183218012030568, 0.0146243897078082, 0.017169910307361626, 0.019283769593320142,
                0.021937205521272526, 0.023025869187442795, 0.024330363510410607, 0.02606651914697333, 0.028396209089729577, 0.029376957427500704, 0.029980002115876404, 0.029126189870868208,
                0.029135700253050564, 0.03141143253160715, 0.03233822320517935, 0.03461395548373612, 0.03737530674008709, 0.038248139749459814, 0.041872813633001334, 0.04717017510672393,
                0.05306107088663986, 0.06073256958513601, 0.0677026186490395, 0.07186686917457506, 0.08342269090982822, 0.09857089576818341, 0.1124891479629985, 0.12805509788700778,
                0.14111486997705686, 0.15509758142888513, 0.17079504360751582, 0.18266659431158058, 0.1964045883967792, 0.2138909081763603, 0.2305787881427304, 0.24618168171585192,
                0.26136111625894864, 0.2750981840082206, 0.2898533566968226, 0.30865732757661113, 0.3237457650538275, 0.33657829664840033, 0.3517439699566726, 0.36473359783181697,
                0.3713539971895155, 0.3559271074424945, 0.3397656255860199, 0.32455847044031233, 0.30761106545597855, 0.2908044489513116, 0.2732227242907315, 0.2525696757406509,
                0.23130825452672196, 0.21043202540567976, 0.18971940727993006, 0.16900678915418035, 0.1477006131427772, 0.1346452163525201, 0.12633606028338817, 0.11075837848640897,
                0.09659077295673812, 0.08431395519706662, 0.07898165682350908, 0.08201279640085744, 0.08261584108923295, 0.08424408139739709, 0.08651981367595367, 0.09230279412747083,
                0.0969526636308008, 0.1013867024773332, 0.1061984449732614, 0.11009290717779996, 0.11679316792070678, 0.12646110019458012, 0.13580528648325715, 0.14763152532510634,
                0.1588102721965627, 0.17107929836385274, 0.18526465780773635, 0.19926072908215398, 0.211900582803943, 0.22512380429512846, 0.23737236600195252, 0.250839360511753,
                0.2648375932366665, 0.2705126583597849, 0.26777032786779825, 0.2634632251140293, 0.26956995155074265, 0.28274513199757656, 0.29464426575740227, 0.30748636177041966,
                0.320020659513345, 0.3327453038503046, 0.34722305608318593, 0.3612757504988503, 0.37625981050235424, 0.39254603121084447, 0.40673712090015407, 0.42163927278389784,
                0.43661406818054277, 0.45007145109738383, 0.4635549539426627, 0.4775800394306586, 0.4910440652053573, 0.504981173654138, 0.5184133456525551, 0.530780805148745,
                0.5441418113895824, 0.556341344194323, 0.5692254032248266, 0.5826283724831498, 0.5961161323566437, 0.6096038922301378, 0.6234625906936828, 0.636775959875559,
                0.6494719083840292, 0.6630016312750092, 0.6753439150192319, 0.689647722808013, 0.7036774061815361, 0.7175114242489115, 0.7304610536249093, 0.7447921068402795,
                0.7596550587139702, 0.7737994329394016, 0.7901336644381743, 0.8062367408144614, 0.822562392555571, 0.8394815786028738, 0.8555846549791609, 0.871687731355448,
                0.8881405072635847, 0.9019640330326955, 0.9143570877051752, 0.9286206400223621, 0.9415930390785291, 0.9566699745513206, 0.9710747699389533, 0.9844964264807075,
                0.9941091042200839, 1.0, 0.9916452078510781, 0.9767380620342111, 0.9636918890701992, 0.9474394384525532, 0.930950072275306, 0.9173005123932566, 0.9034250294712715,
                0.8894079818344215, 0.8734063006447936, 0.8571354056495819, 0.8415072128130362, 0.8260768029894932, 0.8105724484005896, 0.7952164773882346, 0.7797863145876053,
                0.764356151786976, 0.7489259889863468, 0.7323087575733308, 0.7148756326343488, 0.6972197161854946, 0.6775113283738716, 0.6562118668320194, 0.6359429188423249,
                0.6182620556523296, 0.6010758043841623, 0.5820844357463731, 0.5639091767795998, 0.5456597260245524, 0.5275586588460532, 0.5096801670323765, 0.491430716277329,
                0.4732153237112231, 0.4549848397795798, 0.43735009010665404, 0.41880387219850984, 0.40010927071381736, 0.38134047744085064, 0.3624233005913358, 0.34358031553009494,
                0.3254050565633217, 0.3090740243543703, 0.2937284051459977, 0.2777047080391751, 0.2618293945089009, 0.24285175108176107, 0.22676602034570434, 0.21168201384200458,
                0.1961467582304747, 0.18082255282094437, 0.16551744316487726, 0.1504087180486383, 0.13312566261759137, 0.1188591347035548, 0.10560656122926523, 0.0922563622409077,
                0.07806085084144082, 0.06481880157254234, 0.053941056332199755, 0.046180663069669545, 0.038744015792335655, 0.032980056105183016, 0.027917546052622498, 0.02328669731365707,
                0.021083943463664306, 0.017748078665484338, 0.016678435763678782, 0.015500877533474266, 0.016427668207046282};

    std::vector<G4double> TPB_Emission_Energy = {2.067768179441798*eV, 2.0711968936209915*eV, 2.0770256356539587*eV, 2.0828860695279343*eV, 2.0887801678705675*eV, 2.0947076470021497*eV,
                2.1006678523807305*eV, 2.1066628714606246*eV, 2.1126923525118997*eV, 2.118755810633005*eV, 2.124855895993544*eV, 2.130989524274953*eV, 2.137160160957121*eV,
                2.143366186498423*eV, 2.149607830530791*eV, 2.1558859330094995*eV, 2.1622011203809492*eV, 2.168553799156252*eV, 2.1749408782209616*eV, 2.1813703262619417*eV,
                2.187834433136552*eV, 2.1943390922112487*eV, 2.2008819646150983*eV, 2.207464696419716*eV, 2.214085997362992*eV, 2.2207459466791235*eV, 2.2274478952350636*eV,
                2.2341894644722657*eV, 2.2409727628071736*eV, 2.2477960957054406*eV, 2.2546617273091436*eV, 2.261570917109327*eV, 2.2690790780362478*eV, 2.275918360488418*eV,
                2.2825494435685756*eV, 2.289631000619393*eV, 2.303927281249834*eV, 2.311143279273195*eV, 2.3184047089117805*eV, 2.325712708799065*eV, 2.333065857382064*eV,
                2.3404658292859812*eV, 2.347913614534111*eV, 2.3554084066762226*eV, 2.362952022880465*eV, 2.3705443898139276*eV, 2.378186166708452*eV, 2.385877557921112*eV,
                2.39361904820099*eV, 2.401411694943464*eV, 2.4092551517643392*eV, 2.417150873493553*eV, 2.425098902855604*eV, 2.43309985823439*eV, 2.4411542700042*eV, 2.449263166570215*eV,
                2.457426608281983*eV, 2.4656452468071404*eV, 2.4739206455057636*eV, 2.48225258717489*eV, 2.4906410440138873*eV, 2.499087410736674*eV, 2.5075915684931345*eV, 2.516154734144066*eV,
                2.524777002173984*eV, 2.5334591966532396*eV, 2.5422021559960752*eV, 2.5510061368014743*eV, 2.5598747617369026*eV, 2.568805416910652*eV, 2.5778014308884507*eV,
                2.586448802813165*eV, 2.5947403760514343*eV, 2.60308674603487*eV, 2.6110657010388754*eV, 2.6186715711698088*eV, 2.6258957646546524*eV, 2.6331602742724325*eV,
                2.6404664670251035*eV, 2.647380461440109*eV, 2.654330759152101*eV, 2.6613174874481835*eV, 2.668341734678859*eV, 2.675405092457309*eV, 2.682505782323715*eV,
                2.6891971411609465*eV, 2.6954737815622134*eV, 2.7017794146330925*eV, 2.708115750686812*eV, 2.7144816876207662*eV, 2.720877052523635*eV, 2.727302815076956*eV,
                2.7337590002976824*eV, 2.7402446659873796*eV, 2.7472273584415867*eV, 2.7542457284081454*eV, 2.761299020465093*eV, 2.768861193875619*eV, 2.776940442151452*eV,
                2.7850656971082666*eV, 2.7937193127095*eV, 2.8029112074486324*eV, 2.812650284678237*eV, 2.82343910796717*eV, 2.834310754990234*eV, 2.845262870646866*eV,
                2.856296349903497*eV, 2.867410217057209*eV, 2.8786039935020593*eV, 2.889876350742903*eV, 2.901227832437383*eV, 2.9126571798818017*eV, 2.923639997514584*eV,
                2.9331100815988824*eV, 2.9410433360216066*eV, 2.94794805046771*eV, 2.9543467728840396*eV, 2.960233620679133*eV, 2.965604256124921*eV, 2.9709909183145484*eV,
                2.975854365324016*eV, 2.980733761153649*eV, 2.9856262956544306*eV, 2.989986205611953*eV, 2.9943600784321798*eV, 2.998741909609724*eV, 3.0031365839873723*eV,
                3.0075433438647448*eV, 3.0114068881054026*eV, 3.0147270690390626*eV, 3.0180529393848383*eV, 3.0213850604061303*eV, 3.0247239982518894*eV, 3.0280714245291565*eV,
                3.0314251651541078*eV, 3.0365693520234185*eV, 3.0405792321553435*eV, 3.047060787862289*eV, 3.0517296494432844*eV, 3.0558209563241894*eV, 3.0582096913284937*eV,
                3.0625349984167625*eV, 3.066538957624989*eV, 3.0699762348999835*eV, 3.0734217933695915*eV, 3.0768790717507026*eV, 3.080344706475243*eV, 3.084399397230364*eV,
                3.08963270184922*eV, 3.0948868131593352*eV, 3.100162631540684*eV, 3.106048991810478*eV, 3.112555022893491*eV, 3.1196879734878573*eV, 3.12745887942332*eV,
                3.136480280348734*eV, 3.1473828561184702*eV, 3.160208944447646*eV, 3.1739727659683186*eV, 3.1853374734348137*eV, 3.2000033206912755*eV, 3.213926754660582*eV,
                3.227974780964764*eV, 3.2421463258878163*eV, 3.2564438935098026*eV, 3.270870748178594*eV, 3.2854275907735015*eV, 3.300108826645378*eV, 3.314929516983376*eV,
                3.329882241236572*eV, 3.344970105618417*eV, 3.3601940264973393*eV, 3.375558090309753*eV, 3.391063864615245*eV, 3.4067127492143876*eV, 3.42250673454375*eV,
                3.439807832718487*eV, 3.4545383507470926*eV, 3.4707797607643514*eV, 3.4871748082163454*eV, 3.503726487179263*eV, 3.5204366472344355*eV, 3.5373071658780413*eV,
                3.554338294847484*eV, 3.5715342173450098*eV, 3.588896493574171*eV, 3.606428826421865*eV, 3.6241326517720798*eV, 3.64354981506753*eV};

    std::vector<G4double> TPB_Emission = {0.006779687128117056, 0.0074420158190842516, 0.008922515245951421, 0.009078357290884739, 0.009779646493085329, 0.010403014672818593,
                0.009935488538018645, 0.01032509365035207, 0.010870540807619343, 0.01074153574520864, 0.012428961256953306, 0.012351040234486514, 0.013831539661353953,
                0.014844512953421176, 0.015312039088221392, 0.015779565223021607, 0.01655877544768872, 0.017727590784689527, 0.015838005989871833, 0.018584722031823434,
                0.017883432829622844, 0.01928601123402349, 0.02011847121186657, 0.02165914174347512, 0.022300259840934332, 0.02179031251911761, 0.02301929896085953,
                0.02333920680715724, 0.02441495601369449, 0.02428378523805227, 0.024733561455127076, 0.026569468563742243, 0.025800594399193195, 0.02881693952655687,
                0.027078113480694107, 0.029026139042362026, 0.03338971630049755, 0.036194873109299114, 0.03907795094056693, 0.04266231797403587, 0.04531163273790385,
                0.048116789546704876, 0.051545314535240236, 0.05450631338897511, 0.05816860144491031, 0.062064652568245886, 0.06635030880391461, 0.07079180708451693,
                0.07538914741005256, 0.08060985591532199, 0.08575264339812463, 0.09159672008312786, 0.09775248085799745, 0.10429784674520128, 0.11123281774473798, 0.11894699896894208,
                0.12705078530547959, 0.13562209777681758, 0.14544014660762264, 0.15588156361816125, 0.16647882267363343, 0.17785529195377245, 0.18946552430131183, 0.20177704585105183,
                0.2144002514906582, 0.22749098326506503, 0.24120508321920567, 0.2552620356721996, 0.2718280450486214, 0.2885031438564966, 0.30720418924850657, 0.32597536356073636,
                0.345175103496533, 0.36540340092889023, 0.3846317119062579, 0.40396478336938646, 0.42227405917621624, 0.4408214269961391, 0.46030925826777785, 0.47895186289293756,
                0.49759446751809705, 0.5161299307373652, 0.5350939595802, 0.5553436852937356, 0.5754862696013793, 0.5944656043593414, 0.6129551498332273, 0.631199800665075,
                0.6501791354230371, 0.6690360228599803, 0.6875255683338662, 0.706137561128771, 0.724749553923676, 0.7426268627924663, 0.7614837502294093, 0.7803406376663525,
                0.7985546766679454, 0.8171972812931049, 0.8363398791457594, 0.8546967733552079, 0.8731679517309406, 0.892053410209455, 0.9106245872306868, 0.9294814746676301,
                0.9481825200596397, 0.9647796978450484, 0.9792730080238559, 0.9905715562815285, 0.9978961323933989, 1.0, 0.9967273170563984, 0.9869092682255932, 0.9712237664030451,
                0.9523025941225669, 0.9335222362612609, 0.9141959674270778, 0.894339093535145, 0.8744536486016414, 0.8555967611646982, 0.8348541849840606, 0.8154830187988373,
                0.7961118526136141, 0.7751978501835499, 0.7552695486876895, 0.7359840956271794, 0.714127248825268, 0.6922704020233568, 0.6699849895978784, 0.648913846439173,
                0.631199800665075, 0.6126286236438432, 0.5934860257911887, 0.574057717522823, 0.55520083008588, 0.5357725218175142, 0.5147523983759001, 0.4885486716778104,
                0.4599878340501025, 0.4319780093669811, 0.4100599389045599, 0.3874071845160244, 0.3665605281125206, 0.3443465432909704, 0.32377539335975963, 0.30348995384426003,
                0.2852044872387397, 0.26720473104893033, 0.24699071913735865, 0.2232053270293963, 0.20091991460391792, 0.18052019092213378, 0.15943476224264294, 0.1387778991867189,
                0.11834960446336369, 0.09941108547907451, 0.08081439859929712, 0.06232179194238603, 0.04562331682777075, 0.03172837020577989, 0.024150443787635374, 0.019130169189090173,
                0.015389960110688182, 0.012974408414220045, 0.010636777740218968, 0.008766673201018106, 0.008065383998817783, 0.008065383998817783, 0.005549648130606819, 0.006351121504550237,
                0.006429042527016761, 0.006351121504550237, 0.005727753324816705, 0.00549399025741633, 0.00549399025741633, 0.00549399025741633, 0.00549399025741633, 0.005722558589985408,
                0.0055719112798828544, 0.00549399025741633, 0.00549399025741633, 0.0058835953697500215, 0.006506963549483285, 0.007208252751683876, 0.007208252751683876,
                0.007208252751683876, 0.006896568661816977, 0.00674072661688366, 0.006351121504550237, 0.008355856254790841};

    // based on best fit from Fig 13 in https://arxiv.org/pdf/1709.05002, using 0.605 as mean number of photons
    G4double wls_mean_num_photons = 0.605;

    // we control the WLS response as a function of incident energy using WLSABSLENGTH
    // Using Fig 13 best fit for VUV photons and Fig 11 for 250 onwards from https://arxiv.org/pdf/1709.05002
    // with some scaling for the absorption in the visible
    // let's figure out the scaling based on wavelength then convert to energy
    std::vector<G4double> TPB_WLSAbsLength_Wavelength = {276.7824815485308*nm, 284.47440666839077*nm, 292.1663317882507*nm, 300.1496177081054*nm, 306.1516501879961*nm, 312.44504346788153*nm,
                320.1369685877414*nm, 327.38390630397316*nm, 330.7036536008824*nm, 338.31788250741045*nm, 346.00980762727033*nm, 354.1679100271218*nm, 360.694391947003*nm, 368.22494800840434*nm,
                372.2397397895239*nm, 375.520222014443*nm, 378.9752009980987*nm, 382.16185569061213*nm, 384.4694332265701*nm, 386.9168639465255*nm, 389.0766639478518*nm, 392.09778871734034*nm,
                395.5001293094871*nm, 398.45475162631544*nm, 401.0187266662688*nm, 403.0699066982314*nm, 407.19557562615626*nm, 410.34227226609903*nm, 413.13933594604805*nm, 415.8198553059992*nm,
                418.5222909662401*nm, 421.8801599458889*nm, 425.02685658583164*nm, 427.82392026578066*nm, 430.97061690572343*nm, 433.76768058567245*nm, 436.5647442656216*nm, 439.3618079455706*nm,
                442.5085045855134*nm, 445.65520122545604*nm, 447.01643225177565*nm, 500.0*nm, 550.0*nm, 600.0*nm};

    // the absorption spectrum below ~380nm is fixed!!!
    std::vector<G4double> TPB_WLSAbsLength_FixedAbsorption = {104.57106779183636*nm, 104.57106779183636*nm, 101.86164369840671*nm, 95.71852757662245*nm, 80.96895332836402*nm, 60.211922218905926*nm,
                46.01246188212347*nm, 35.985035667274715*nm, 29.88492679795619*nm, 29.83741776132185*nm, 29.80183555478305*nm, 31.806087107658087*nm, 42.023350161001645*nm, 55.528067293309526*nm,
                67.12665530148168*nm, 90.18763865616405*nm, 113.91022322401744*nm, 182.59722816676643*nm};

    // but we allow different exponential behavior the absorption of visible light
    // now loop over the wavelength and grab/calculate correct absorption and then convert wavelength to energy
    std::vector<G4double> TPB_WLSAbsLength_Energy = {1.0*eV}; // padding lower bound
    std::vector<G4double> TPB_WLSAbsLength = {1e35*nm}; // padding lower bound (doesnt really matter, just need super big absorption length)
    // note -- looping over backwards to go from low energy (high wavelength) to high energy (low wavelength)
    for (size_t w = (TPB_WLSAbsLength_Wavelength.size() - 1); w > 0; w --){
        G4double this_wavelength = TPB_WLSAbsLength_Wavelength.at(w) / nm;
        G4double this_abs;
        if (w < TPB_WLSAbsLength_FixedAbsorption.size()){
            this_abs = TPB_WLSAbsLength_FixedAbsorption.at(w);
        } else {
            // ok now we need to calculate our absorption length
            this_abs = TPBAbsNorm_ * exp(TPBAbsTau_ * this_wavelength) * nm;
        }

        // multiply our absorption by the overall scaling
        this_abs *= TPBAbsScale_;

        // now convert wavelength to energy and save!
        double this_energy = ((197.326 * 2.0 * M_PI) / this_wavelength) * eV; // hc / wavelength (units are hardcoded -- energy in ev and wavelength in nm)
        TPB_WLSAbsLength_Energy.push_back(this_energy);
        TPB_WLSAbsLength.push_back(this_abs);
        std::cout << "for wavelength = " << this_wavelength << "nm, energy = " << this_energy/eV << "eV, abs = " << this_abs/nm << "nm" << std::endl;
    }

    // pad our spectrum once more...
    TPB_WLSAbsLength_Energy.push_back(12.0 * eV);
    TPB_WLSAbsLength.push_back(TPB_WLSAbsLength_FixedAbsorption.at(0));

    // I dont think this makes a difference, but adding index of refraction just in case
    std::vector<G4double> tpb_rin_energy = {1.0*eV, 14.0*eV};
    std::vector<G4double> tpb_rin = {1.62, 1.62};

    // let's try adding some mie scattering to our tpb! i have no idea what this is gonna do...we'll see
    std::vector<G4double> TPB_Mie_Scattering_Energy = {1.0 * eV, 4.0 * eV, 4.1 * eV, 12.0 * eV}; // 1 - 4 eV are vis, 4.1 - 12 eV are UV
    std::vector<G4double> TPB_Mie_Scattering_Length = {0.05 * mm, 0.05 * mm, 1.0 * m, 1.0 * m}; // mie scattering for vis light, no mie scattering for UV

    // now make our tpb foils!
    // making different ones for sides and top/bottom :)

    // side cylinder of TPB foil
    G4MaterialPropertiesTable* fTPBFoilSides_mt = new G4MaterialPropertiesTable();
    fTPBFoilSides_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission);
    fTPBFoilSides_mt->AddConstProperty("WLSTIMECONSTANT", 0.00001*ns); // setting to very small at the moment
    std::cout << "setting wls mean number of photons to " << WLSNPhotonsSideFoil_ << " for side tpb foils" << std::endl;
    fTPBFoilSides_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", WLSNPhotonsSideFoil_);
    fTPBFoilSides_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBFoilSides_mt->AddProperty("RINDEX", tpb_rin_energy, tpb_rin);
    // mie scattering!
    fTPBFoilSides_mt->AddProperty("MIEHG", TPB_Mie_Scattering_Energy, TPB_Mie_Scattering_Length);
    fTPBFoilSides_mt->AddConstProperty("MIEHG_FORWARD", Mie_GG_);
    fTPBFoilSides_mt->AddConstProperty("MIEHG_BACKWARD", Mie_GG_);
    fTPBFoilSides_mt->AddConstProperty("MIEHG_FORWARD_RATIO", Mie_Ratio_);
    fTPBFoilSides->SetMaterialPropertiesTable(fTPBFoilSides_mt);

    // top/bottom faces of tpb foil -- these have WLSNPhotonsFoil_!!!
    G4MaterialPropertiesTable* fTPBFoilTopBottom_mt = new G4MaterialPropertiesTable();
    fTPBFoilTopBottom_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission);
    fTPBFoilTopBottom_mt->AddConstProperty("WLSTIMECONSTANT", 0.00001*ns); // setting to very small at the moment
    std::cout << "setting wls mean number of photons to " << WLSNPhotonsEndCapFoil_ << " for top/bottom tpb foils" << std::endl;
    fTPBFoilTopBottom_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", WLSNPhotonsEndCapFoil_);
    fTPBFoilTopBottom_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBFoilTopBottom_mt->AddProperty("RINDEX", tpb_rin_energy, tpb_rin);
    // mie scattering!
    fTPBFoilTopBottom_mt->AddProperty("MIEHG", TPB_Mie_Scattering_Energy, TPB_Mie_Scattering_Length);
    fTPBFoilTopBottom_mt->AddConstProperty("MIEHG_FORWARD", Mie_GG_);
    fTPBFoilTopBottom_mt->AddConstProperty("MIEHG_BACKWARD", Mie_GG_);
    fTPBFoilTopBottom_mt->AddConstProperty("MIEHG_FORWARD_RATIO", Mie_Ratio_);
    fTPBFoilTopBottom->SetMaterialPropertiesTable(fTPBFoilTopBottom_mt);

    // tpb on pmts
    G4MaterialPropertiesTable* fTPBPMT_mt = new G4MaterialPropertiesTable();
    fTPBPMT_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission);
    fTPBPMT_mt->AddConstProperty("WLSTIMECONSTANT", 0.00001*ns); // setting to very small at the moment
    std::cout << "setting wls mean number of photons to " << WLSNPhotonsPMT_ << " for pmt tpb foils" << std::endl;
    fTPBPMT_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", WLSNPhotonsPMT_);
    fTPBPMT_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBPMT_mt->AddProperty("RINDEX", tpb_rin_energy, tpb_rin);
    // mie scattering!
    fTPBPMT_mt->AddProperty("MIEHG", TPB_Mie_Scattering_Energy, TPB_Mie_Scattering_Length);
    fTPBPMT_mt->AddConstProperty("MIEHG_FORWARD", Mie_GG_);
    fTPBPMT_mt->AddConstProperty("MIEHG_BACKWARD", Mie_GG_);
    fTPBPMT_mt->AddConstProperty("MIEHG_FORWARD_RATIO", Mie_Ratio_);
    fTPBPMT->SetMaterialPropertiesTable(fTPBPMT_mt);

    // Defines properties of the reflectors.
    std::vector<G4double> TefEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV,
                                       2.138*eV, 2.25*eV,  2.38*eV,  2.48*eV,
                                       2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,
                                       3.124*eV, 3.457*eV, 3.643*eV, 3.812*eV, 4.086*eV,
                                       4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                       7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };

    std::vector<G4double> TefRIndex = {1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64,
                                       1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64, 1.64,
                                        1.64, 1.64, 1.64, 1.64, 1.64};

    std::vector<G4double> ptfe_energy = {7.0*eV, 7.07*eV, 7.14*eV};
    std::vector<G4double> ptfe_AbsLength = {1e-12*mm, 1e-12*mm, 1e-12*mm};

    G4MaterialPropertiesTable* fPTFE_mt = new G4MaterialPropertiesTable();
    fPTFE_mt->AddProperty("RINDEX", TefEnergy, TefRIndex);
    fPTFE_mt->AddProperty("ABSLENGTH", ptfe_energy, ptfe_AbsLength);
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
        fMainVolume = new G4CCMMainVolume(nullptr, G4ThreeVector(), fExperimentalHall_log, false, 0, this,
                                          SourceRodIn_, SourceRodLocation_, CobaltSourceRun_, SodiumSourceRun_,
                                          EndCapFoilTPBThickness_, SideFoilTPBThickness_, PMTTPBThickness_);
    }

    return fExperimentalHall_phys;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// construct our pmts SD
void G4CCMDetectorConstruction::ConstructSDandField() {
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
            pmt_SD->SetReadout(readout_);
            G4SDManager::GetSDMpointer()->AddNewDetector(fPMT_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTCoatedWall(), fPMT_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTCoatedCaps(), fPMT_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTUncoatedCaps(), fPMT_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTUncoatedWall(), fPMT_SD.Get());
        }

    }

    if (LArSDStatus_){
        // Scint SD
        if(!fScint_SD.Get()) {
            G4cout << "Construction /LAr/scintSD" << G4endl;
            auto scint_SD = new G4CCMScintSD("/LAr/scintSD");
            scint_SD->SetPMTSDStatus(PMTSDStatus_);
            scint_SD->SetTimeCutStatus(TimeCut_);
            scint_SD->SetKillCherenkovStatus(KillCherenkov_);
            scint_SD->SetReadout(readout_);
            fScint_SD.Put(scint_SD);
            G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBCoatingWall(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBCoatingCaps(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBFoilSides(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBFoilTop(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBFoilBottom(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTCoatedWall(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTCoatedCaps(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTUncoatedWall(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTUncoatedCaps(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogReflectorFoil(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetShinyC406R0(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetShinyTop(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetShinyBottom(), fScint_SD.Get());
            // make sure to include source pellet + rod for SD if enabeled
            if (SourceRodIn_){
                SetSensitiveDetector(fMainVolume->GetLogSourceRod(), fScint_SD.Get());
            }
            if (CobaltSourceRun_ or SodiumSourceRun_){
                SetSensitiveDetector(fMainVolume->GetLogSourcePellet(), fScint_SD.Get());
            }
        }

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
    G4MTRunManager::GetRunManager()->ReinitializeGeometry();
}


void G4CCMDetectorConstruction::SetReadout(G4CCMReadout * readout) {
   readout_ = readout;
    G4CCMPMTSD* pmt = fPMT_SD.Get();
    if(pmt) {
        pmt->SetReadout(readout_);
    }
    G4CCMScintSD * scint = fScint_SD.Get();
    if(scint) {
        scint->SetReadout(readout_);
    }
}
