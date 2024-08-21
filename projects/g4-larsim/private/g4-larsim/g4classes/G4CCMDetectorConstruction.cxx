

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
G4CCMDetectorConstruction::G4CCMDetectorConstruction(G4double SingletTau, G4double TripletTau, G4double UVAbsLength, G4double Rayleigh128) {
    SingletTau_ = SingletTau;
    TripletTau_ = TripletTau;
    UVAbsLength_ = UVAbsLength;
    Rayleigh128_ = Rayleigh128;
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
  if(fTPBFoil != nullptr) delete fTPBFoil;
  if(fGlass_mt != nullptr) delete fGlass_mt;
  if(fPlastic_mt != nullptr) delete fPlastic_mt;
  if(fBlackPlastic_mt != nullptr) delete fBlackPlastic_mt;
  if(fLAr_mt != nullptr) delete fLAr_mt;
  if(fAlum_mt != nullptr) delete fAlum_mt;
  if(fSteel_mt != nullptr) delete fSteel_mt;
  if(fVacuum_mt != nullptr) delete fVacuum_mt;
  if(fPTFE_mt != nullptr) delete fPTFE_mt;
  if(fTPBFoil_mt != nullptr) delete fTPBFoil_mt;
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
    fTPBFoil = new G4Material("TPBFoil", density= 1.079*g/cm3, 2);
    fTPBFoil->AddElement(fC, 28);
    fTPBFoil->AddElement(fH, 22);

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

    std::vector<G4double> alum_abslen = { 1.0e-3*cm, 1.0e-3*cm};
    std::vector<G4double> alum_abseneg = { 1.5498*eV, 20.664*eV};

    G4MaterialPropertiesTable* fAlum_mt = new G4MaterialPropertiesTable();
    fAlum_mt->AddProperty("REFLECTIVITY", alum_energy, alum_reflect);
    fAlum_mt->AddProperty("ABSLENGTH", alum_abseneg, alum_abslen);
    fAlum->SetMaterialPropertiesTable(fAlum_mt);

    // now time to define foil + pmt TPB
    // at the moment, getting the same properties, but one day might assign different properties

    // Emission spectrum of TPB is controlled by "WLSCOMPONENT"
    // using energy and normalized spectrum pulled from lower right plot in Fig. 4 (axes are emission wavelength vs intensity) from https://arxiv.org/pdf/1210.3793
    // note -- noralized the intensity to integrate to one
    // also padding spectrum with 0.0 on either end
    std::vector<G4double> TPB_Emission_Energy = {14.0*eV, 6.191827936846959*eV, 5.985915350705116*eV, 5.797096561555699*eV, 5.619825649392834*eV, 5.453074671411217*eV, 5.29593417314788*eV, 5.147596594361829*eV, 5.007335004365769*eV, 4.874528346529596*eV, 4.748577742204919*eV, 4.628971952001313*eV, 4.515243321093366*eV, 4.406969051582301*eV, 4.30376595016281*eV, 4.205285896902584*eV, 4.1112094336510205*eV, 4.021254798859424*eV, 3.9351500718048067*eV, 3.852655461651022*eV, 3.7735485879663107*eV, 3.697624967454319*eV, 3.6246962490169308*eV, 3.5545849330007373*eV, 3.4871415918821045*eV, 3.422206435859456*eV, 3.3596454258854105*eV, 3.299330694912653*eV, 3.2411442800215085*eV, 3.188161671147352*eV, 3.1650988602873977*eV, 3.143008879989636*eV, 3.126053146367847*eV, 3.1140603757294354*eV, 3.1045536357451335*eV, 3.0927311249607796*eV, 3.0810005119229347*eV, 3.0716920690275646*eV, 3.0624410589818924*eV, 3.0509524010395817*eV, 3.0418267773138186*eV, 3.035031683887979*eV, 3.026002440867842*eV, 3.014775940025159*eV, 3.0014295188687776*eV, 2.9896455058927986*eV, 2.9860222588842675*eV, 2.9772916471557287*eV, 2.9642464379640887*eV, 2.9556287466209663*eV, 2.943843976767722*eV, 2.9300487609963777*eV, 2.9181210155948856*eV, 2.901263781362373*eV, 2.869912416373187*eV, 2.8273936530372312*eV, 2.7824018812493883*eV, 2.7439574872768895*eV, 2.7307391797993827*eV, 2.70544267579332*eV, 2.683719966807585*eV, 2.6628862604766415*eV, 2.652536961389666*eV, 2.634547901929533*eV, 2.6186539079884374*eV, 2.603447114507513*eV, 2.5917228450898007*eV, 2.571896476516571*eV, 2.555605444519763*eV, 2.5379231746909836*eV, 2.5236283984124293*eV, 2.512620790605399*eV, 2.498613217988123*eV, 2.4808292949218753*eV, 2.468041727171535*eV, 2.4545208082744616*eV, 2.435265951301829*eV, 2.417765505500488*eV, 2.399079069815766*eV, 2.3792869076076886*eV, 2.3598204192783827*eV, 2.3366072022471*eV, 2.3085649215936406*eV, 2.2799024796941434*eV, 2.2519482645733313*eV, 2.224675827928228*eV, 2.1980564645366782*eV, 2.1720678198737606*eV, 2.1466890592211167*eV, 2.1218984607956104*eV, 2.0976766090773373*eV, 2.075065069877949*eV, 1.0};

    std::vector<G4double> TPB_Emission_Normalized_Intensity = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0005622062600842365, 0.001503443199733779, 0.0026926709198067202, 0.00391638544559267, 0.004946015400875132, 0.006256319802282813, 0.007427329763164467, 0.008648983308141366, 0.00990864412545381, 0.011199957182825777, 0.01274446357640978, 0.01406109842582936, 0.015366503729714128, 0.016719714500980264, 0.018044765363484924, 0.020021405757084138, 0.021304732838133247, 0.02269521490746361, 0.024290389469406978, 0.025572157686232836, 0.026822322831527422, 0.0281135990125032, 0.02993660824221481, 0.03136630612234021, 0.03315213779863045, 0.034577574563899136, 0.035347364552504305, 0.035614350090406414, 0.03451754987676052, 0.03336382675093531, 0.032206445245944675, 0.03057497963396541, 0.02922781229918319, 0.028185856262140697, 0.02673140671075951, 0.025163638920597316, 0.024167427175013936, 0.022795658015919158, 0.021390003546106935, 0.02024183646397373, 0.019017679421433554, 0.017726046768627976, 0.0167113148774046, 0.015573301096354564, 0.01403743080355351, 0.0131925609585961, 0.01188573523056196, 0.010583875523880844, 0.009646674206949535, 0.008278977840944609, 0.007244057586158524, 0.0062790737094087275, 0.005132904976745019, 0.004124581383703269, 0.0031618124320831067, 0.002424637627796259, 0.0018900371598903282, 0.0013738525407463174, 0.0009144030554235511, 0.0005743442364584392, 0.0003292414245426639, 0.00021870894738393895, 6.332776967286204e-05, 0.0};
    
    // based on best fit from Fig 13 in https://arxiv.org/pdf/1709.05002, using 0.605 as mean number of photons
    G4double wls_mean_num_photons = 0.605;

    // we control the WLS response as a function of incident energy using WLSABSLENGTH
    // Using Fig 13 best fit for VUV photons and Fig 11 for 250 onwards from https://arxiv.org/pdf/1709.05002
    std::vector<G4double> TPB_WLSAbsLength_Energy = {12.39835823924519*eV, 4.520784914393557*eV, 4.403204903803938*eV, 4.314967417290872*eV, 4.236201285563993*eV, 4.157961750495743*eV, 4.085069083133975*eV, 4.0443177306841145*eV, 3.9817633828435217*eV, 3.9225321515018465*eV, 3.9205742918328146*eV, 3.8448631783865683*eV, 3.8025712231071136*eV, 3.7652970652514997*eV, 3.762862815865301*eV, 3.703191726710966*eV, 3.6412533430538785*eV, 3.5831817960708783*eV, 3.5290710713006566*eV, 3.48851021668889*eV, 3.4699110804911135*eV, 3.4361647119545338*eV, 3.3817661197965347*eV, 3.3459286055866837*eV, 3.3214274015068077*eV, 3.2755707229674424*eV, 3.2586507642702967*eV, 3.243867614987933*eV, 3.2225605676408287*eV, 3.2182531899038254*eV, 3.2039575153849125*eV, 3.1855614510225836*eV, 3.170859339074369*eV, 3.165441570938464*eV, 3.1390881863227924*eV, 3.1120639907749723*eV, 3.0930600046021595*eV, 3.079724230541368*eV, 3.0627397742099034*eV, 3.0456037796888635*eV, 3.032631087118836*eV, 3.019467387075414*eV, 2.999845732513397*eV, 2.983711394762438*eV, 2.971691329508609*eV, 2.96551484559134*eV, 2.9466766809674496*eV, 2.936310020028604*eV, 2.917122828423819*eV, 2.898127963726499*eV, 2.878918640363774*eV, 2.859390696395847*eV, 2.8394316919309692*eV, 2.820180762358608*eV, 2.8017086319312554*eV, 2.7875621787203357*eV, 2.7762497220260878*eV, 1.5497947799056488*eV, 1.0331965199370992*eV};

    std::vector<G4double> TPB_WLSAbsLength = {425.0*nm, 103.83329986337642*nm, 107.98073112246271*nm, 106.40677312247695*nm, 101.76697185101447*nm, 101.05172421787127*nm, 90.58455101512752*nm, 73.72675461199395*nm, 61.95847711489636*nm, 49.85890980317873*nm, 58.14673257532595*nm, 41.70021143667599*nm, 38.80279864198152*nm, 37.82386620617285*nm, 30.431746756020132*nm, 29.071521590291532*nm, 28.82888338356318*nm, 29.324649474231183*nm, 29.91384462882476*nm, 31.802265335516566*nm, 40.38391252461217*nm, 37.222095886705205*nm, 53.50121319296904*nm, 63.289058767089*nm, 88.9206377946115*nm, 108.02360627084354*nm, 153.5032805631189*nm, 199.85823080844185*nm, 255.85822744954626*nm, 200.48621550244906*nm, 353.05033779572574*nm, 480.5566616748119*nm, 629.3261989758672*nm, 791.1551509826332*nm, 1036.4431202910403*nm, 1573.3202940562826*nm, 2195.402325808518*nm, 2904.7461061188183*nm, 3531.1084066403014*nm, 4938.293668651343*nm, 5931.319058258765*nm, 7701.344462077818*nm, 11224.3487112599*nm, 15220.323839686669*nm, 18401.537167766*nm, 21435.851020054517*nm, 27132.275975142325*nm, 35858.47052646484*nm, 51842.61999393237*nm, 75873.60320297175*nm, 109030.37626309693*nm, 158887.77188864764*nm, 239537.8727919214*nm, 354962.3984602222*nm, 513931.57477543503*nm, 709624.1774047639*nm, 861003.9340108575*nm, 4.0740917253503145e+35*nm, 4.0740917253503145e+35*nm};

    
    G4MaterialPropertiesTable* fTPBFoil_mt = new G4MaterialPropertiesTable();
    fTPBFoil_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission_Normalized_Intensity);
    fTPBFoil_mt->AddConstProperty("WLSTIMECONSTANT", 0.00001*ns); // setting to very small at the moment
    fTPBFoil_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", wls_mean_num_photons);
    fTPBFoil_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBFoil->SetMaterialPropertiesTable(fTPBFoil_mt);
    
    G4MaterialPropertiesTable* fTPBPMT_mt = new G4MaterialPropertiesTable();
    fTPBPMT_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission_Normalized_Intensity);
    fTPBPMT_mt->AddConstProperty("WLSTIMECONSTANT", 0.00001*ns); // setting to very small at the moment
    fTPBPMT_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", wls_mean_num_photons);
    fTPBPMT_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBPMT->SetMaterialPropertiesTable(fTPBPMT_mt);

    // Defines properties of the PTFE reflectors.
    std::vector<G4double> TefEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV,
                                       2.138*eV, 2.25*eV,  2.38*eV,  2.48*eV,
                                       2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,
                                       3.124*eV, 3.457*eV, 3.643*eV, 3.812*eV, 4.086*eV,
                                       4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                       7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };

    std::vector<G4double> TefRIndex = {10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
                                       10., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
                                        10., 10., 10., 10., 10.};


    G4MaterialPropertiesTable* fPTFE_mt = new G4MaterialPropertiesTable();
    fPTFE_mt->AddProperty("RINDEX", TefEnergy, TefRIndex);
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
        fMainVolume = new G4CCMMainVolume(nullptr, G4ThreeVector(), fExperimentalHall_log, false, 0, this, SourceRodIn_, SourceRodLocation_, CobaltSourceRun_, SodiumSourceRun_);
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
            fScint_SD.Put(scint_SD);
            G4SDManager::GetSDMpointer()->AddNewDetector(fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBCoatingWall(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBCoatingCaps(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogTPBFoil(), fScint_SD.Get());

            SetSensitiveDetector(fMainVolume->GetLogPMTCoatedWall(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTCoatedCaps(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTUncoatedWall(), fScint_SD.Get());
            SetSensitiveDetector(fMainVolume->GetLogPMTUncoatedCaps(), fScint_SD.Get());


            SetSensitiveDetector(fMainVolume->GetLogReflectorFoil(), fScint_SD.Get());
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



