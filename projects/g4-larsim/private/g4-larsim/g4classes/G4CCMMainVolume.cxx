
#include "g4-larsim/g4classes/G4CCMMainVolume.h"

#include <map>
#include <math.h> 

#include "dataclasses/I3Position.h"
#include "dataclasses/I3Orientation.h"
#include "dataclasses/I3Map.h"
#include "icetray/CCMPMTKey.h"
#include "dataclasses/geometry/CCMGeometry.h"

#include "globals.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4OpticalSurface.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

// let's define some things relevant for getting the geometry of our pmts
std::vector<std::string> position_id = {"C101R0",
                                        "C401R0",
                                        "C105R0",
                                        "C210R0",
                                        "C201R0",
                                        "C103R0",
                                        "C315R0",
                                        "C301R0",
                                        "C102R0",
                                        "C314R0",
                                        "C419R0",
                                        "C313R0",
                                        "C209R0",
                                        "C417R0",
                                        "C104R0",
                                        "C312R0",
                                        "C207R0",
                                        "C208R0",
                                        "C415R0",
                                        "C311R0",
                                        "C302R0",
                                        "C305R0",
                                        "C303R0",
                                        "C204R0",
                                        "C306R0",
                                        "C405R0",
                                        "C304R0",
                                        "C403R0",
                                        "C202R0",
                                        "C203R0",
                                        "C407R0",
                                        "C307R0",
                                        "C205R0",
                                        "C409R0",
                                        "C308R0",
                                        "C309R0",
                                        "C411R0",
                                        "C206R0",
                                        "C310R0",
                                        "C413R0",
                                        "C2R3",
                                        "C2R4",
                                        "C2R2",
                                        "C2R5",
                                        "C2R1",
                                        "C1R3",
                                        "C1R4",
                                        "C1R5",
                                        "C1R1",
                                        "C1R2",
                                        "C3R3",
                                        "C3R5",
                                        "C3R1",
                                        "C3R4",
                                        "C3R2",
                                        "C5R2",
                                        "C5R5",
                                        "C5R1",
                                        "C5R4",
                                        "C5R3",
                                        "C24R5",
                                        "C24R2",
                                        "C24R4",
                                        "C24R1",
                                        "C24R3",
                                        "C23R4",
                                        "C23R3",
                                        "C23R1",
                                        "C23R2",
                                        "C23R5",
                                        "C4R5",
                                        "C4R3",
                                        "C4R4",
                                        "C4R2",
                                        "C4R1",
                                        "C8R1",
                                        "C8R5",
                                        "C8R4",
                                        "C8R3",
                                        "C8R2",
                                        "C10R1",
                                        "C9R4",
                                        "C9R5",
                                        "C9R2",
                                        "C10R3",
                                        "C10R2",
                                        "C10R5",
                                        "C9R3",
                                        "C9R1",
                                        "C10R4",
                                        "C11R1",
                                        "C11R3",
                                        "C11R2",
                                        "C11R4",
                                        "C11R5",
                                        "C13R5",
                                        "C12R1",
                                        "C12R3",
                                        "C12R2",
                                        "C13R1",
                                        "C13R2",
                                        "C13R4",
                                        "C12R5",
                                        "C12R4",
                                        "C13R3",
                                        "C14R1",
                                        "C14R2",
                                        "C14R5",
                                        "C14R4",
                                        "C14R3",
                                        "C15R3",
                                        "C15R2",
                                        "C15R1",
                                        "C15R4",
                                        "C15R5",
                                        "C16R3",
                                        "C16R2",
                                         "C16R5",
                                        "C16R4",
                                        "C16R1",
                                        "C19R1",
                                        "C19R5",
                                        "C19R4",
                                        "C19R2",
                                        "C19R3",
                                        "C20R1",
                                        "C20R2",
                                        "C20R5",
                                        "C20R4",
                                        "C20R3",
                                        "C21R3",
                                        "C21R4",
                                        "C21R5",
                                        "C21R2",
                                        "C21R1",
                                        "C22R2",
                                        "C22R4",
                                        "C22R1",
                                        "C22R3",
                                        "C22R5",
                                        "C18R2",
                                        "C18R1",
                                        "C18R4",
                                        "C17R3",
                                        "C18R3",
                                        "C17R2",
                                        "C17R5",
                                        "C18R5",
                                        "C17R4",
                                        "C17R1",
                                        "C6R4",
                                        "C6R5",
                                        "C6R2",
                                        "C7R2",
                                        "C7R3",
                                        "C6R3",
                                        "C7R4",
                                        "C6R1",
                                        "C7R5",
                                        "C7R1",
                                        "C410R6",
                                        "C308R6",
                                        "C206R6",
                                        "C309R6",
                                        "C412R6",
                                        "C310R6",
                                        "C207R6",
                                        "C416R6",
                                        "C313R6",
                                        "C209R6",
                                        "C208R6",
                                        "C312R6",
                                        "C311R6",
                                        "C105R6",
                                        "C104R6",
                                        "C101R6",
                                        "C314R6",
                                        "C418R6",
                                        "C414R6",
                                        "C315R6",
                                        "C210R6",
                                        "C420R6",
                                        "C201R6",
                                        "C301R6",
                                        "C202R6",
                                        "C404R6",
                                        "C303R6",
                                        "C302R6",
                                        "C402R6",
                                        "C203R6",
                                        "C304R6",
                                        "C102R6",
                                        "C103R6",
                                        "C307R6",
                                        "C204R6",
                                        "C406R6",
                                        "C205R6",
                                        "C408R6",
                                        "C306R6",
                                        "C305R6"};

// Z positions of the rows of pmts
std::map<size_t, double> pmt_region_z_positions = {
    {-2, 75.00 * I3Units::cm},
    {-1, 65.00 * I3Units::cm},
    {0,  58.00 * I3Units::cm},
    {1,  46.22 * I3Units::cm},
    {2,  23.11 * I3Units::cm},
    {3,   0.00 * I3Units::cm},
    {4, -23.11 * I3Units::cm},
    {5, -46.22 * I3Units::cm},
    {6, -58.00 * I3Units::cm},
    {7, -65.00 * I3Units::cm},
    {8, -75.00 * I3Units::cm},
};

// The radii of rings of pmts
std::map<size_t, double> ring_radii = {
    {0,  96.0 * I3Units::cm}, // Outer wall
    {1,  20.0 * I3Units::cm}, // Inner ring on top/bottom
    {2,  42.0 * I3Units::cm}, //
    {3,  64.2 * I3Units::cm}, //
    {4,  85.5 * I3Units::cm}, // Outer ring on top/bottom
    {5, 101.6 * I3Units::cm},  // Veto VT and VB rings
    {6, 111.6 * I3Units::cm}  // Veto VCT and VCB rings
};

// The number of possible pmt positions in each ring (not all positions will be filled)
std::map<size_t, size_t> ring_pmt_pos_count = {
    {0, 24 * I3Units::cm}, // Outer wall
    {1,  5 * I3Units::cm}, // Inner ring on top/bottom
    {2, 10 * I3Units::cm}, //
    {3, 15 * I3Units::cm}, //
    {4, 20 * I3Units::cm}, // Outer ring on top/bottom
    {5, 24 * I3Units::cm},  // Veto VT and VB rings
    {6, 24 * I3Units::cm}  // Veto VCT and VCB rings
};

std::set<std::tuple<size_t, size_t, size_t>> cap_uncoated_pmts = {
    // Ring, col, row
    {1,  1, 0},
    {1,  3, 0},
    {2,  3, 0},
    {2,  8, 0},
    {3,  1, 0},
    {3,  8, 0},
    {4,  5, 0},
    {4, 15, 0},
    {1,  2, 6},
    {1,  4, 6},
    {2,  5, 6},
    {2, 10, 6},
    {3,  3, 6},
    {3, 10, 6},
    {4,  8, 6},
    {4, 18, 6}
};

double default_start_pmt_number = 1;
int default_num_pmts = 24;
int default_middle_pmt_row = 3;

double angle_from_cylinder_pmt_number(
        int pmt_number,
        double start_pmt_number = default_start_pmt_number,
        int num_pmts = default_num_pmts) {
    return double(pmt_number - start_pmt_number) / double(num_pmts) * 2.0 * M_PI;
}

std::pair<double, double> xy_position_from_cylinder_pmt_number(
        int pmt_number,
        double radius,
        double start_pmt_number = default_start_pmt_number,
        int num_pmts = default_num_pmts,
        double angular_offset = 0.0) {
    double angle = angle_from_cylinder_pmt_number(pmt_number, start_pmt_number, num_pmts);
    angle += angular_offset;
    return {radius * cos(angle), radius * sin(angle)};
}

I3Position get_pmt_cap_position(int pmt_row, int ring_number, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
    double z = pmt_region_z_positions[pmt_row];
    std::pair<double, double> xy = xy_position_from_cylinder_pmt_number(
            pmt_number, 
            ring_radii[ring_number], 
            starting_pmt_number, 
            ring_pmt_pos_count[ring_number],
            angular_offset);
    I3Position pos(xy.first, xy.second, z);
    return pos;
}

I3Position get_pmt_wall_position(int pmt_row, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
    return get_pmt_cap_position(pmt_row, 0, pmt_number, starting_pmt_number, angular_offset);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMMainVolume::G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                                 G4LogicalVolume* pMotherLogical, G4bool pMany,
                                 G4int pCopyNo, G4CCMDetectorConstruction* c)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(pRot, tlate,
                  // Temp logical volume must be created here
                  new G4LogicalVolume(new G4Box("temp", 1, 1, 1),
                                    G4Material::GetMaterial("Vacuum"), "temp"),
                  "housing", pMotherLogical, pMany, pCopyNo)
  , fConstructor(c) {

    // now let's build our detector
    // Outer cryogen
    fCryoVessel = new G4Tubs("Cryogen", 0*cm, 138*cm, 131*cm, 0*deg,360*deg);
    fCryoVessel_log  = new G4LogicalVolume(fCryoVessel, G4Material::GetMaterial("Steel"), "Cryogen");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fCryoVessel_log, "Cryogen", pMotherLogical, false, 0, true);

    // Vacuum jacket
    fVacuum = new G4Tubs("Vacuum", 0*cm, 135*cm, 126*cm, 0*deg, 360*deg);
    fVacuum_log = new G4LogicalVolume(fVacuum, G4Material::GetMaterial("Vacuum"), "Vacuum");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fVacuum_log, "Vacuum", fCryoVessel_log, false, 0, true);

    // Inner cryogen
    fInnerJacket = new G4Tubs("InnerJacket", 0*cm, 125*cm, 120*cm, 0*deg, 360*deg);
    fInnerJacket_log = new G4LogicalVolume(fInnerJacket, G4Material::GetMaterial("Steel"), "InnerJacket");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fInnerJacket_log, "InnerJacket", fVacuum_log, false, 0, true);

    // Argon outside the fiducial volume
    fArgonOuter = new G4Tubs("OuterLiquidArgon", 0*cm, 120*cm, 115*cm, 0*deg, 360*deg);
    fArgonOuter_log = new G4LogicalVolume(fArgonOuter, G4Material::GetMaterial("LAr"), "OuterLiquidArgon");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fArgonOuter_log, "OuterLiquidArgon", fInnerJacket_log, false, 0, true);

    // Aluminum frame holding PMTs and instrumentation
    fInnerFrame = new G4Tubs("InnerFrame",0*cm, 106*cm, 75*cm, 0*deg, 360*deg);
    fInnerFrame_log= new G4LogicalVolume(fInnerFrame, G4Material::GetMaterial("Al"), "InnerFrame");
    new G4PVPlacement(0, G4ThreeVector(0,0,0), fInnerFrame_log, "InnerFrame", fArgonOuter_log, false, 0, true);

    // now let's place the TPB foils
    // note -- using the same thickness for TPB foils across the entire detector, if this isn't true, fix at some point
    G4double basethick = 0.00019*cm;
    G4double totalH = 123.2/2.0;

    fTPBFoil = new G4Tubs("TPBFoil", 0*cm, (96.0*cm+basethick), (totalH*cm+basethick), 0*deg, 360*deg);
    fTPBFoil_log = new G4LogicalVolume(fTPBFoil, G4Material::GetMaterial("TPBFoil"), "TPBFoil");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fTPBFoil_log, "TPBFoil", fInnerFrame_log, false, 0, true);

    // now fiducial LAr!
    fFiducialAr = new G4Tubs("FiducialArgon", 0*cm, 96*cm, totalH*cm, 0*deg, 360*deg);
    fFiducialAr_log = new G4LogicalVolume(fFiducialAr, G4Material::GetMaterial("LAr"), "FiducialArgon");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fFiducialAr_log, "FiducialArgon", fTPBFoil_log, false, 0, true);

    // now let's build PMTs using J4SolidMaker
    fPMT = J4PMTSolidMaker::Get8inchPMTSolid();  
    fPMT_log = new G4LogicalVolume(fPMT, G4Material::GetMaterial("Glass"), "pmt_log");
    
    // now let's build also PMT photocathods using J4SolidMaker (fix at some point)
    fPhotocath = J4PMTSolidMaker::Get8inchPMTSolid();  
    fPhotocath_log = new G4LogicalVolume(fPhotocath, G4Material::GetMaterial("Alum"), "photocath_log");
    G4double height_pmt = 10.16; // idk, fix
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., -height_pmt / 2.), fPhotocath_log, "photocath", fPMT_log, false, 0);


    // now that we've defined the pmt logical volume, we can get pmt locations using CCMGeometryGenerator logic
    G4int k = 0;
    for (size_t i = 0; i < position_id.size(); i++){
        std::string position_string = position_id[i]; 
        size_t r_pos = position_string.find("r");
        size_t n_row_chars = r_pos - 1;
        bool on_caps = n_row_chars > 2;
        int row;
        I3Position position;
        CCMOMGeo::OMType omtype = CCMOMGeo::OMType::CCM8inCoated;
        int pmt_number = 0;
        if(on_caps) {
            int col = std::atoi(position_string.substr(2, n_row_chars - 1).c_str());
            int ring = std::atoi(position_string.substr(1, 1).c_str());
            row = std::atoi(position_string.substr(r_pos + 1, std::string::npos).c_str());
            position = get_pmt_cap_position(row, ring, col);
            for(std::pair<const int, int> const & p : ring_pmt_pos_count)
                if(p.first < ring and p.first > 0)
                    pmt_number += p.second;
            pmt_number += col;
            std::tuple<int, int, int> cap_coated_key = {ring, col, row};
            if(cap_uncoated_pmts.count(cap_coated_key) > 0)
                omtype = CCMOMGeo::OMType::CCM8inUncoated;
        } else {
            int col = std::atoi(position_string.substr(1, r_pos - 1).c_str());
            row = std::atoi(position_string.substr(r_pos + 1, std::string::npos).c_str());
            position = get_pmt_wall_position(row, col);
            pmt_number = col;
            if( (row ==  1 and col % 4 == 1) or
                (row ==  2 and col % 4 == 3) or
                (row ==  4 and col % 4 == 0) or
                (row ==  5 and col % 4 == 2)) {
                omtype = CCMOMGeo::OMType::CCM8inUncoated;
            }
        }
        // so we have our pmt info
        // let's make the string key to save it to
        G4String pmt_name = std::to_string(row) + "_" + std::to_string(pmt_number);
        new G4PVPlacement(0, G4ThreeVector(position.GetX()/I3Units::cm, position.GetY()/I3Units::cm, position.GetZ()/I3Units::cm), fPMT_log, pmt_name, fFiducialAr_log, false, k);
        ++k;
        fPMTPositions.push_back(G4ThreeVector(position.GetX()/I3Units::cm, position.GetY()/I3Units::cm, position.GetZ()/I3Units::cm));
    }

    VisAttributes();
    SurfaceProperties();
    SetLogicalVolume(fFiducialAr_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::VisAttributes()
{
    auto argon_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
    fFiducialAr_log->SetVisAttributes(argon_va);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::SurfaceProperties()
{
    // now it's time to define optical surface properties
    // let's start with the TPB foils
    
    std::vector<G4double> TPBEnergy = { 0.602*eV/*(2066nm)*/, 0.689*eV/*(1799nm)*/, 1.030*eV/*(1204nm)*/, 1.926*eV/*(644nm)*/, 2.138*eV/* (580nm)*/,
                                        2.250*eV/*(551nm)*/,  2.380*eV/*(521nm)*/,  2.480*eV/*(500nm)*/,  2.583*eV/*(480nm)*/, 2.800*eV/*(443nm)*/,
                                        2.880*eV/*(431nm)*/,  2.980*eV/*(416nm)*/,  3.124*eV/*(397nm)*/,  3.457*eV/*(359nm)*/, 3.643*eV/*(341nm)*/,
                                        3.812*eV/*(325nm)*/,  4.086*eV/*(304nm)*/,  4.511*eV/*(275nm)*/,  5.166*eV/*(240nm)*/, 5.821*eV/*(213nm)*/,
                                        6.526*eV/*(190nm)*/,  8.266*eV/*(150nm)*/,  9.686*eV/*(128nm)*/,  11.27*eV/*(110nm)*/, 12.60*eV/*(98nm)*/  };
    
    std::vector<G4double> TPBfoilOSTransmit = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                                               1., 1., 1., 1., 1., 1., 1., 1.}; // set to 1 and have all absorption in bulk

    std::vector<G4double> TPBfoilOSReflect = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                              0.0, 0.0, 0.0, 0., 0., 0., 0., 0.};

    std::vector<G4double> TPBfoilOSEff = {0., 0., 0., 0., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0, 1.0, 1., 1., 1., 1., 1., 1.,
                                          1., 1., 1., 1., 1.};

    // define our TPB foil material properties table
    G4MaterialPropertiesTable *TPBfoilMPT = new G4MaterialPropertiesTable();
    TPBfoilMPT->AddProperty("REFLECTIVITY", TPBEnergy, TPBfoilOSReflect);
    TPBfoilMPT->AddProperty("TRANSMITTANCE", TPBEnergy, TPBfoilOSTransmit);
    TPBfoilMPT->AddProperty("EFFICIENCY", TPBEnergy, TPBfoilOSEff);
    
    // now define our optical surface
    G4OpticalSurface *TPBfoilOS = new G4OpticalSurface("TPBFoilsSurface");
    TPBfoilOS->SetModel(unified); 
    TPBfoilOS->SetType(dielectric_dielectric);
    TPBfoilOS->SetFinish(ground);
    TPBfoilOS->SetSigmaAlpha(0.05);
    TPBfoilOS->SetMaterialPropertiesTable(TPBfoilMPT);

    // create logical skin surfaces
    new G4LogicalSkinSurface("TPBFoils_Surface", fTPBFoil_log, TPBfoilOS);

    //**Photocathode surface properties
    std::vector<G4double> ephoton = { 7.0 * eV, 7.14 * eV };
    std::vector<G4double> photocath_EFF     = { 1., 1. };
    std::vector<G4double> photocath_ReR     = { 1.92, 1.92 };
    std::vector<G4double> photocath_ImR     = { 1.69, 1.69 };
    auto photocath_mt = new G4MaterialPropertiesTable();
    photocath_mt->AddProperty("EFFICIENCY", ephoton, photocath_EFF);
    photocath_mt->AddProperty("REALRINDEX", ephoton, photocath_ReR);
    photocath_mt->AddProperty("IMAGINARYRINDEX", ephoton, photocath_ImR);
    auto photocath_opsurf = new G4OpticalSurface(
    "photocath_opsurf", glisur, polished, dielectric_metal);
    photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

    //**Create logical skin surfaces
    new G4LogicalSkinSurface("photocath_surf", fPhotocath_log, photocath_opsurf);

}

