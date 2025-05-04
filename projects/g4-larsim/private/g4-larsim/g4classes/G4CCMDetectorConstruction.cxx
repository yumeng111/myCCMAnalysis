

#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4CCMMainVolume.h"
#include "g4-larsim/g4classes/G4CCMDetectorMessenger.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

#include <sstream>
#include <cmath>
#include <boost/preprocessor/stringize.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <functional>
#include <numeric>
#include <vector>

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
G4CCMDetectorConstruction::G4CCMDetectorConstruction(G4bool EnableUVAbsorption, G4double UVAbsA, G4double UVAbsB, G4double UVAbsD, G4double UVAbsScaling,
                                                     G4double WLSNPhotonsEndCapFoil, G4double WLSNPhotonsSideFoil, G4double WLSNPhotonsPMT,
                                                     G4double EndCapFoilTPBThickness, G4double SideFoilTPBThickness, G4double PMTTPBThickness,
                                                     G4double Rayleigh128, G4double TPBAbsTau, G4double TPBAbsNorm, G4double TPBAbsScale,
                                                     G4double Mie_GG, G4double Mie_Ratio, G4double Normalization, G4double PhotonSamplingFactor) {
    EnableUVAbsorption_ = EnableUVAbsorption;
    UVAbsA_ = UVAbsA;
    UVAbsB_ = UVAbsB;
    UVAbsD_ = UVAbsD;
    UVAbsScaling_ = UVAbsScaling;
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
    Normalization_ = Normalization;
    PhotonSamplingFactor_ = PhotonSamplingFactor;
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

double G4CCMDetectorConstruction::HarmonicOscillatorRefractiveIndex(double a0, double aUV, double gammaUV, double wlUV, double wl) {
    return a0 + aUV * ((std::pow(wlUV,-2.0) - std::pow(wl,-2.0)) / (std::pow((std::pow(wlUV,-2.0) - std::pow(wl,-2.0)),2.0) + (std::pow(gammaUV,2.0) * std::pow(wl,-2.0))));
}

double G4CCMDetectorConstruction::HarmonicOscillatorRefractiveIndexDerivative(double aUV, double lambda, double lambdaUV, double gammaUV) {
    double lambda2 = lambda * lambda;
    double lambda4 = lambda2 * lambda2;
    double lambdaUV2 = lambdaUV * lambdaUV;
    double lambdaUV4 = lambdaUV2 * lambdaUV2;
    double gammaUV2 = gammaUV * gammaUV;

    double numerator = 2.0 * aUV * lambda * lambdaUV4 * (2.0 * lambda2 * lambdaUV2 - lambdaUV4 + lambda4 * (-1.0 + gammaUV2 * lambdaUV2));

    double denominator = std::pow(lambda4 + lambdaUV4 + lambda2 * lambdaUV2 * (-2.0 + gammaUV2 * lambdaUV2), 2.0);

    return numerator / denominator;
}

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
    G4MaterialPropertiesTable* fLAr_mt = new G4MaterialPropertiesTable();

    // let's grab our scintillation profile from txt file
    // txt file is digitization of liquid argon (black line) in fig 1 from https://arxiv.org/pdf/1511.07718
    std::vector<G4double> lar_wavelength = {};
    std::vector<G4double> lar_intensity = {};

    std::string sourceDir = std::filesystem::path(__FILE__).parent_path().string();
    std::ifstream file(sourceDir + "/LArScintillationProfile.txt");

    // Check if the file is successfully opened
    if (!file.is_open()) {
        log_fatal("Cannot open file with LAr scintillation profile!");
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        double wavelength, intensity;
        char comma;

        ss >> wavelength >> comma >> intensity;
        lar_wavelength.push_back(wavelength);
        lar_intensity.push_back(intensity);
    }

    file.close();

    // before we are done, let's normalizse the intensity so it sums to 1
    double intensity_total = std::accumulate(lar_intensity.begin(), lar_intensity.end(), 0.0);
    for (size_t i = 0; i < lar_intensity.size(); i++){
        lar_intensity.at(i) /= intensity_total;
    }

    // finally, loop over backward and convert to energy
    std::vector<G4double> lar_energy = {};
    std::vector<G4double> lar_sorted_intensity = {};

    for (size_t w = lar_wavelength.size(); w > 0; w --){
        double this_wavelength = lar_wavelength.at(w-1);
        double this_intensity = lar_intensity.at(w-1);
        double this_energy = ((197.326 * 2.0 * M_PI) / this_wavelength) * eV; // hc / wavelength (units are hardcoded -- energy in ev and wavelength in nm)

        lar_energy.push_back(this_energy);
        lar_sorted_intensity.push_back(this_intensity);
    }

    fLAr_mt->AddProperty("SCINTILLATIONCOMPONENT1", lar_energy, lar_sorted_intensity);

    // for LAr index of refraction, using grace fit around section 2.3.3 from https://arxiv.org/pdf/2408.00817v1
    // can also calculate rayleigh scattering length at the same time
    std::vector<G4double> lar_light_properties_energy = {};
    std::vector<G4double> ho_rin_vals = {};
    std::vector<G4double> grace_rin_vals = {};
    std::vector<G4double> babicz_rin_vals = {};
    std::vector<G4double> group_velocity_vals = {};
    std::vector<G4double> rayl_scattering_length = {};

    // constants for index of refraction
    double a0 = 1.26;
    double aUV = 0.23;
    double aIR = 0.0023;
    double lamUV = 106.6;
    double lamIR = 908.3;
    double a0bab = 0.334;
    double aUVbab =  0.100;
    double aIRbab = 0.008;

    double a0ho = 1.10232;
    double aUVho = 0.00001058;
    //double gammaho = 0.002524;
    double gammaho = 0.002794965;

    double starting_wavelength = 100.0;
    double ending_wavelength = 850.0;
    size_t n_entries = 10000;

    // for rayl, we have the desired scattering length (in cm) at 128nm so first we need to solve for the scaling
    double rindex_128 = HarmonicOscillatorRefractiveIndex(a0ho, aUVho, gammaho, lamUV, 128.0);
    double rayl_scaling = (Rayleigh128_ / cm) * std::pow((std::pow(rindex_128, 2.0) - 1)*(std::pow(rindex_128, 2.0) + 2), 2.0) / std::pow(128.0*1e-7, 4.0);

    for(int i = (n_entries + 1); i >= 0; --i) {
        double this_wavelength = starting_wavelength + ((static_cast<double>(i-1) / static_cast<double>(n_entries)) * (ending_wavelength - starting_wavelength));
        double this_energy = ((197.326 * 2.0 * M_PI) / this_wavelength) * eV; // hc / wavelength (units are hardcoded -- energy in ev and wavelength in nm)
        double this_rindex = HarmonicOscillatorRefractiveIndex(a0ho, aUVho, gammaho, lamUV, this_wavelength);
        double this_rayl = (rayl_scaling * std::pow(this_wavelength*1e-7, 4.0) / std::pow((std::pow(this_rindex, 2.0) - 1)*(std::pow(this_rindex, 2.0) + 2), 2.0) ) * cm;

        double dn_dlambda = HarmonicOscillatorRefractiveIndexDerivative(aUVho, this_wavelength, lamUV, gammaho);
        double this_group_velocity = (c_light / this_rindex) * (1.0 + ((this_wavelength / this_rindex) * dn_dlambda));

        // save
        lar_light_properties_energy.push_back(this_energy);
        ho_rin_vals.push_back(this_rindex);
        group_velocity_vals.push_back(this_group_velocity);
        rayl_scattering_length.push_back(this_rayl);
    }

    fLAr_mt->AddProperty("RINDEX", lar_light_properties_energy, ho_rin_vals);
    fLAr_mt->AddProperty("GROUPVEL", lar_light_properties_energy, group_velocity_vals);
    fLAr_mt->AddProperty("RAYLEIGH", lar_light_properties_energy, rayl_scattering_length);

    // using this paper: https://link.springer.com/article/10.1140/epjc/s10052-012-2190-z#Bib1
    // set parameter for uv abs scaling function
    double a_param = UVAbsA_ / (1.0/nm); // per nm
    double b_param = UVAbsB_ / nm; // nm
    double d_param = UVAbsD_ / cm; // cm
    double scaling_param = UVAbsScaling_; // dimensionless

    double large_abs_length = 100000.0; // cm

    // now let's get some min and max bounds
    double min_wavelength = b_param + 0.1;
    double max_wavelength = std::min(800.0, b_param - log(1.0 - std::exp(-d_param/large_abs_length)) / a_param);

    double min_wavelength_function_T = 1.0 - std::exp( - a_param * (min_wavelength - b_param));
    double min_abs_length = (d_param / std::log(1.0 / min_wavelength_function_T));
    min_abs_length *= scaling_param;

    double max_wavelength_function_T = 1.0 - std::exp( - a_param * (max_wavelength - b_param));
    double max_abs_length = std::min((d_param / std::log(1.0 / max_wavelength_function_T)), large_abs_length /*cm*/);
    max_abs_length *= scaling_param;

    // now we need to fill in our absorption lengths
    std::vector<G4double> uv_abs_energy; uv_abs_energy.reserve(n_entries+1);
    std::vector<G4double> uv_abs_length; uv_abs_length.reserve(n_entries+1);
    for(int i = (n_entries + 1); i >= 0; --i) {
        double this_wavelength = starting_wavelength + ((static_cast<double>(i-1) / static_cast<double>(n_entries)) * (ending_wavelength - starting_wavelength));
        G4double this_energy = ((197.326 * 2.0 * M_PI) / this_wavelength) * eV; // hc / wavelength (units are hardcoded -- energy in ev and wavelength in nm)
        G4double this_abs;

        if(this_wavelength <= b_param) {
            this_abs = min_abs_length * cm;
        } else if (this_wavelength > b_param and this_wavelength < max_wavelength){
            double function_T = 1.0 - std::exp( - a_param * (this_wavelength - b_param));
            double abs_length = (d_param / std::log(1.0 / function_T));
            this_abs = scaling_param * abs_length * cm;
        } else {
            this_abs = max_abs_length * cm;
        }

        uv_abs_energy.push_back(this_energy);
        uv_abs_length.push_back(this_abs);
    }

    if(EnableUVAbsorption_)
        fLAr_mt->AddProperty("ABSLENGTH", uv_abs_energy, uv_abs_length);

    std::cout << "using normalization = " << Normalization_ << " for scintillation yields" << std::endl;
    G4double scint_yeild = Normalization_ * (1.0/(19.5*eV)) * PhotonSamplingFactor_; // scintillation yield: 50 per keV.
    fLAr_mt->AddConstProperty("SCINTILLATIONYIELD", scint_yeild);
    fLAr_mt->AddConstProperty("RESOLUTIONSCALE", std::sqrt(PhotonSamplingFactor_));
    fLAr_mt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 1e-5 * ns);
    fLAr_mt->AddConstProperty("SCINTILLATIONYIELD1",1.0); // for e/m scintillation
    fLAr->SetMaterialPropertiesTable(fLAr_mt);

    // Set the Birks Constant for the LAr scintillator
    fLAr->GetIonisation()->SetBirksConstant(0.145731*cm/MeV);
    fLAr->GetIonisation()->SetMeanExcitationEnergy(203.0 * eV);

    // Set PMT glass constants
    std::vector<G4double> glass_energy = {1.0*eV, 5.0*eV, 7.07*eV, 10.14*eV, 12.0*eV};
    std::vector<G4double> glass_AbsLength = {1e-12*mm, 1e-12*mm, 1e-12*mm, 1e-12*mm, 1e-12*mm};
    std::vector<G4double> glass_RIND = {1.49,1.49,1.49, 1.49, 1.49};

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
    fAlum_mt->AddProperty("ABSLENGTH", alum_abseneg, alum_abslen);
    fAlum->SetMaterialPropertiesTable(fAlum_mt);

    // now time to define foil + pmt TPB
    // emission spectrum digitized from https://arxiv.org/pdf/1709.05002
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
    }

    // pad our spectrum once more...
    TPB_WLSAbsLength_Energy.push_back(12.0 * eV);
    TPB_WLSAbsLength.push_back(TPB_WLSAbsLength_FixedAbsorption.at(0));

    // tpb rindex taken from discussion on pg 10 of https://arxiv.org/pdf/1709.05002
    std::vector<G4double> tpb_rin_energy = {1.0*eV, 2.0*eV, 3.0*eV, 4.0*eV, 5.0*eV, 6.0*eV, 7.0*eV,
                                            8.0*eV, 9.0*eV, 10.0*eV, 11.0*eV, 12.0*eV, 13.0*eV, 14.0*eV};
    std::vector<G4double> tpb_rin = {1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67,
                                     1.67, 1.67, 1.67, 1.67, 1.67, 1.67, 1.67};

    // let's try adding some mie scattering to our tpb! i have no idea what this is gonna do...we'll see
    std::vector<G4double> TPB_Scattering_Energy = {1.0 * eV, 2.0 * eV, 3.1 * eV, 3.11 * eV, 5.9*eV, 6.0 * eV, 12.0*eV}; // ~1200 nm, ~600nm, ~400nm, ~401nm, ~210nm, ~206nm, ~100nm
    //std::vector<G4double> TPB_Mie_Scattering_Length = {0.00075 * mm, 0.00075 * mm, 0.00075 * mm, 1.0 * m, 1.0 * m, 1.0 * m, 1.0 * m}; // mie scattering for vis light
    std::vector<G4double> TPB_Mie_Scattering_Length = {0.0003 * mm, 0.0003 * mm, 0.0003 * mm, 1.0 * m, 1.0 * m, 1.0 * m, 1.0 * m}; // mie scattering for vis light

    // now make our tpb foils!
    // making different ones for sides and top/bottom :)

    // side cylinder of TPB foil
    G4MaterialPropertiesTable* fTPBFoilSides_mt = new G4MaterialPropertiesTable();
    fTPBFoilSides_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission);
    fTPBFoilSides_mt->AddConstProperty("WLSTIMECONSTANT", 0.3*ns); // setting to very small at the moment
    std::cout << "setting wls mean number of photons to " << WLSNPhotonsSideFoil_ << " for side tpb foils" << std::endl;
    fTPBFoilSides_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", WLSNPhotonsSideFoil_);
    fTPBFoilSides_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBFoilSides_mt->AddProperty("RINDEX", tpb_rin_energy, tpb_rin);
    // mie scattering!
    fTPBFoilSides_mt->AddProperty("MIEHG", TPB_Scattering_Energy, TPB_Mie_Scattering_Length);
    fTPBFoilSides_mt->AddConstProperty("MIEHG_FORWARD", Mie_GG_);
    fTPBFoilSides_mt->AddConstProperty("MIEHG_BACKWARD", Mie_GG_);
    fTPBFoilSides_mt->AddConstProperty("MIEHG_FORWARD_RATIO", Mie_Ratio_);
    fTPBFoilSides->SetMaterialPropertiesTable(fTPBFoilSides_mt);

    // top/bottom faces of tpb foil -- these have WLSNPhotonsFoil_!!!
    G4MaterialPropertiesTable* fTPBFoilTopBottom_mt = new G4MaterialPropertiesTable();
    fTPBFoilTopBottom_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission);
    fTPBFoilTopBottom_mt->AddConstProperty("WLSTIMECONSTANT", 0.3*ns); // setting to very small at the moment
    std::cout << "setting wls mean number of photons to " << WLSNPhotonsEndCapFoil_ << " for top/bottom tpb foils" << std::endl;
    fTPBFoilTopBottom_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", WLSNPhotonsEndCapFoil_);
    fTPBFoilTopBottom_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBFoilTopBottom_mt->AddProperty("RINDEX", tpb_rin_energy, tpb_rin);
    // mie scattering!
    fTPBFoilTopBottom_mt->AddProperty("MIEHG", TPB_Scattering_Energy, TPB_Mie_Scattering_Length);
    fTPBFoilTopBottom_mt->AddConstProperty("MIEHG_FORWARD", Mie_GG_);
    fTPBFoilTopBottom_mt->AddConstProperty("MIEHG_BACKWARD", Mie_GG_);
    fTPBFoilTopBottom_mt->AddConstProperty("MIEHG_FORWARD_RATIO", Mie_Ratio_);
    fTPBFoilTopBottom->SetMaterialPropertiesTable(fTPBFoilTopBottom_mt);

    // tpb on pmts
    G4MaterialPropertiesTable* fTPBPMT_mt = new G4MaterialPropertiesTable();
    fTPBPMT_mt->AddProperty("WLSCOMPONENT", TPB_Emission_Energy, TPB_Emission);
    fTPBPMT_mt->AddConstProperty("WLSTIMECONSTANT", 0.3*ns); // setting to very small at the moment
    std::cout << "setting wls mean number of photons to " << WLSNPhotonsPMT_ << " for pmt tpb foils" << std::endl;
    fTPBPMT_mt->AddConstProperty("WLSMEANNUMBERPHOTONS", WLSNPhotonsPMT_);
    fTPBPMT_mt->AddProperty("WLSABSLENGTH", TPB_WLSAbsLength_Energy, TPB_WLSAbsLength);
    fTPBPMT_mt->AddProperty("RINDEX", tpb_rin_energy, tpb_rin);
    // mie scattering!
    fTPBPMT_mt->AddProperty("MIEHG", TPB_Scattering_Energy, TPB_Mie_Scattering_Length);
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

    // digitized from here : https://www.filmetrics.com/refractive-index-database/PET/Estar-Melinex-Mylar#:~:text=For%20a%20typical%20sample%20of,refractive%20index%20and%20extinction%20coefficients.
    std::vector<G4double> PET_rindex_energy = {1.7778300334109947*eV, 1.7872697079948843*eV, 1.7968107655233225*eV, 1.806453896957342*eV, 1.8162015044994946*eV, 1.826054531660671*eV,
        1.8360155399984264*eV, 1.8460853920026714*eV, 1.8562670291780863*eV, 1.86656123525517*eV, 1.876970547528567*eV, 1.887496314116809*eV, 1.8981413253203288*eV, 1.9089071638391237*eV,
        1.9197955149718025*eV, 1.9308097236821775*eV, 1.9419499434801575*eV, 1.95322057296573*eV, 1.9646224690150118*eV, 1.9761582634916426*eV, 1.9878307402957722*eV, 1.9996419268993484*eV,
        2.011594226376238*eV, 2.023690779616533*eV, 2.0359334373439077*eV, 2.0483257361004954*eV, 2.060869463686366*eV, 2.073568397369212*eV, 2.0864239856362023*eV, 2.0994415308743246*eV,
        2.1126224402804925*eV, 2.1259699025813332*eV, 2.1394870946714972*eV, 2.1531770817267284*eV, 2.167043784983864*eV, 2.181089856068025*eV, 2.195320804465512*eV, 2.2097379661424887*eV,
        2.2243466680731285*eV, 2.2391498133794525*eV, 2.2541515217246206*eV, 2.269355922451397*eV, 2.2847669340950563*eV, 2.3003894578729946*eV, 2.3162270964299503*eV, 2.3322844368517197*eV,
        2.348566771634016*eV, 2.365078047854428*eV, 2.3818234826322553*eV, 2.3988082124844183*eV, 2.4160377682141903*eV, 2.4335172351702834*eV, 2.451251089704145*eV, 2.469246068240786*eV,
        2.4875079809296183*eV, 2.5060429423227224*eV, 2.5248560606325587*eV, 2.5439551267223606*eV, 2.5633465729270943*eV, 2.583035914715776*eV, 2.6030297876112307*eV, 2.623336601263155*eV,
        2.6439643428117443*eV, 2.6649193499467736*eV, 2.6862106749340904*eV, 2.7078452585472292*eV, 2.729833175754107*eV, 2.7521816892141318*eV, 2.7749010022782667*eV, 2.797999465601588*eV,
        2.8204104545663076*eV, 2.8431836870882465*eV, 2.8663310619680407*eV, 2.888730728292982*eV, 2.9103373563879384*eV, 2.9311084674490284*eV, 2.9510009747212393*eV, 2.9699732331440534*eV,
        2.989192041978536*eV, 3.008661720346588*eV, 3.0283885251345364*eV, 3.04837864755239*eV, 3.0686336209552505*eV, 3.083487539137929*eV, 3.0976469047471906*eV};

    std::vector<G4double> PET_rindex_vals = {1.6302321152772228, 1.6304785070225736, 1.6307741771169946, 1.6310424703508208, 1.6313436158173609, 1.6316173844233062, 1.6319294806340838,
        1.632208724612148, 1.6325427223114013, 1.6328493431500601, 1.6331778654771945, 1.6334844863158533, 1.6338294347593445, 1.6341798585749545, 1.6345083809020888, 1.6349026076946502,
        1.6352201792775467, 1.635614406070108, 1.6359867313741936, 1.6363590566782793, 1.6367587588429595, 1.6371584610076397, 1.637552687800201, 1.6379797668254756, 1.6383904197343937,
        1.6388394002481441, 1.6392664792734188, 1.639731885903526, 1.6401480141845628, 1.640657223791621, 1.6411609580265605, 1.6416646922615, 1.6421684264964393, 1.642661209987141,
        1.6431758949663182, 1.6436686784570198, 1.644249067901624, 1.6447911297413957, 1.6453824699302375, 1.6459738101190795, 1.6465761010521591, 1.6471948181015956, 1.6478190105231512,
        1.6484815305495388, 1.6491440505759265, 1.6498120459744332, 1.650518368977772, 1.651224691981111, 1.6519474411008066, 1.652692091708978, 1.6534750699219818, 1.6542854249955798,
        1.6550793539528212, 1.6559061351427762, 1.6567657685654447, 1.6576637295929453, 1.658556215248327, 1.659503454624898, 1.660499972350539, 1.66149649007618, 1.6624820570575831,
        1.6635059516438186, 1.6645900753233622, 1.6656851497471437, 1.666834977892114, 1.6679957567813222, 1.6692277155080764, 1.6704800113312719, 1.6717956650318435, 1.6731426065730945,
        1.6744566958816323, 1.6757817359344076, 1.6772151883551374, 1.6786245491385439, 1.679971003980051, 1.6812893519113479, 1.6825725662048816, 1.6838298485276857, 1.6851172453971435,
        1.6864196995399285, 1.6877748541393578, 1.6892128237420858, 1.6906282074348231, 1.691521200126798, 1.6925655766028962};

    std::vector<G4double> ptfe_energy = {7.0*eV, 7.07*eV, 7.14*eV};
    std::vector<G4double> ptfe_AbsLength = {1e-12*mm, 1e-12*mm, 1e-12*mm};

    G4MaterialPropertiesTable* fPTFE_mt = new G4MaterialPropertiesTable();
    fPTFE_mt->AddProperty("RINDEX", PET_rindex_energy, PET_rindex_vals);
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
                                          TrainingSource_, DecayX_, DecayY_, DecayZ_,
                                          EndCapFoilTPBThickness_, SideFoilTPBThickness_, PMTTPBThickness_);
    }

    return fExperimentalHall_phys;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// construct our pmts SD
void G4CCMDetectorConstruction::ConstructSDandField() {
    if(!fMainVolume)
        return;

    if(RecordHits_) {
        // PMT SD
        G4CCMPMTSD* pmt = fPMT_SD.Get();
        if(!pmt) {
            // Created here so it exists as pmts are being placed
            G4cout << "Construction /LAr/pmtSD" << G4endl;
            auto pmt_SD = new G4CCMPMTSD("/LAr/pmtSD");
            fPMT_SD.Put(pmt_SD);

            pmt_SD->InitPMTs();
            pmt_SD->SetPmtPositions(fMainVolume->GetPMTPositions());
            pmt_SD->SetPhotonTracking(DetailedPhotonTracking_);
            pmt_SD->SetReadout(readout_);
            G4SDManager::GetSDMpointer()->AddNewDetector(fPMT_SD.Get());
            for(G4LogicalVolume * log : fMainVolume->GetPMTLogicalVolumes()) {
                if(log == nullptr)
                    continue;
                SetSensitiveDetector(log, fPMT_SD.Get());
            }
        }
    }

    bool tree_tracker = TrackParticles_ or TrackEnergyLosses_
                     or DetailedPhotonTracking_ or TimeCut_
                     or KillNeutrinos_ or KillCherenkov_ or KillScintillation_ or KillPhotons_
                     or VetoSDSaveEnergyLossesTree_ or InteriorSDSaveEnergyLossesTree_
                     or SaveAllEnergyLossesTree_;

    bool track_particles = TrackParticles_ or TrackEnergyLosses_
                     or SaveAllEnergyLossesTree_ or VetoSDSaveEnergyLossesTree_
                     or InteriorSDSaveEnergyLossesTree_;

    bool track_energy_losses = TrackEnergyLosses_ or SaveAllEnergyLossesTree_
                     or VetoSDSaveEnergyLossesTree_ or InteriorSDSaveEnergyLossesTree_;

    if(tree_tracker) {
        // Tree tracker
        if(!fTreeTracker_SD.Get()) {
            G4cout << "Construction /LAr/treeTracker" << G4endl;
            auto tree_tracker = new G4CCMTreeTracker("/LAr/treeTracker");
            tree_tracker->SetTrackParticles(track_particles);
            tree_tracker->SetTrackEnergyLosses(track_energy_losses);
            tree_tracker->SetDetailedPhotonTracking(DetailedPhotonTracking_);
            tree_tracker->SetTimeCut(TimeCut_);
            tree_tracker->SetKillNeutrinos(KillNeutrinos_);
            tree_tracker->SetKillCherenkov(KillCherenkov_);
            tree_tracker->SetKillScintillation(KillScintillation_);
            tree_tracker->SetKillPhotons(KillPhotons_);
            tree_tracker->SetReadout(readout_);
            tree_tracker->SetG4RangeCut(G4RangeCut_);
            tree_tracker->SetG4EDepMin(G4EDepMin_);
            tree_tracker->SetG4ETrackingMin(G4ETrackingMin_);
            fTreeTracker_SD.Put(tree_tracker);
            G4SDManager::GetSDMpointer()->AddNewDetector(fTreeTracker_SD.Get());
            for(G4LogicalVolume * log : fMainVolume->GetAllLogicalVolumes()) {
                if(log == nullptr)
                    continue;
                SetSensitiveDetector(log, fTreeTracker_SD.Get());
            }
        }
    }
    if(VetoSDSaveEnergyLossesTree_ or VetoSDSaveEnergyLossesVector_) {
        // Veto SD
        if(!fVeto_SD.Get()) {
            G4cout << "Construction /LAr/vetoSD" << G4endl;
            auto veto_SD = new G4CCMEDepSD("/LAr/vetoSD", G4CCMReadout::VolumeType::Veto);
            veto_SD->SetSaveEnergyLossesTree(VetoSDSaveEnergyLossesTree_);
            veto_SD->SetSaveEnergyLossesVector(VetoSDSaveEnergyLossesVector_);
            veto_SD->SetPruneTree(VetoSDPruneTree_);
            veto_SD->SetReadout(readout_);
            if(fTreeTracker_SD.Get())
                veto_SD->SetTreeTracker(fTreeTracker_SD.Get());
            fVeto_SD.Put(veto_SD);
            G4SDManager::GetSDMpointer()->AddNewDetector(fVeto_SD.Get());
            for(G4LogicalVolume * log : fMainVolume->GetVetoLArLogicalVolumes()) {
                if(log == nullptr)
                    continue;
                SetSensitiveDetector(log, fVeto_SD.Get());
            }
        }
    }
    if(InteriorSDSaveEnergyLossesTree_ or InteriorSDSaveEnergyLossesVector_) {
        // Interior SD
        if(!fInterior_SD.Get()) {
            G4cout << "Construction /LAr/interiorSD" << G4endl;
            auto interior_SD = new G4CCMEDepSD("/LAr/interiorSD", G4CCMReadout::VolumeType::Inner);
            interior_SD->SetSaveEnergyLossesTree(InteriorSDSaveEnergyLossesTree_);
            interior_SD->SetSaveEnergyLossesVector(InteriorSDSaveEnergyLossesVector_);
            interior_SD->SetPruneTree(InteriorSDPruneTree_);
            interior_SD->SetTreeTracker(fTreeTracker_SD.Get());
            interior_SD->SetReadout(readout_);
            if(fTreeTracker_SD.Get())
                interior_SD->SetTreeTracker(fTreeTracker_SD.Get());
            fInterior_SD.Put(interior_SD);
            G4SDManager::GetSDMpointer()->AddNewDetector(fInterior_SD.Get());
            for(G4LogicalVolume * log : fMainVolume->GetInteriorLArLogicalVolumes()) {
                if(log == nullptr)
                    continue;
                SetSensitiveDetector(log, fInterior_SD.Get());
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
}
