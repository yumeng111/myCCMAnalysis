
#include "g4-larsim/g4classes/G4CCMMainVolume.h"

#include <map>
#include <set>
#include <math.h>
#include <tuple>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <globals.hh>
#include <G4Colour.hh>
#include <G4Sphere.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4VisAttributes.hh>
#include <G4SystemOfUnits.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4MaterialPropertiesTable.hh>

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
    {-2, 75.00},
    {-1, 65.00},
    {0,  58.00},
    {1,  46.22},
    {2,  23.11},
    {3,   0.00},
    {4, -23.11},
    {5, -46.22},
    {6, -58.00},
    {7, -65.00},
    {8, -75.00},
};

// The radii of rings of pmts
std::map<size_t, double> ring_radii = {
    {0,  96.0}, // Outer wall
    {1,  20.0}, // Inner ring on top/bottom
    {2,  42.0}, //
    {3,  64.2}, //
    {4,  85.5}, // Outer ring on top/bottom
    {5, 101.6},  // Veto VT and VB rings
    {6, 111.6}  // Veto VCT and VCB rings
};

// The number of possible pmt positions in each ring (not all positions will be filled)
std::map<size_t, size_t> ring_pmt_pos_count = {
    {0, 24}, // Outer wall
    {1,  5}, // Inner ring on top/bottom
    {2, 10}, //
    {3, 15}, //
    {4, 20}, // Outer ring on top/bottom
    {5, 24},  // Veto VT and VB rings
    {6, 24}  // Veto VCT and VCB rings
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

std::vector<double> get_pmt_cap_position(int pmt_row, int ring_number, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
    double z = pmt_region_z_positions[pmt_row];
    std::pair<double, double> xy = xy_position_from_cylinder_pmt_number(
            pmt_number,
            ring_radii[ring_number],
            starting_pmt_number,
            ring_pmt_pos_count[ring_number],
            angular_offset);
    std::vector<double> pos;
    pos.push_back(xy.first);
    pos.push_back(xy.second);
    pos.push_back(z);
    return pos;
}

std::vector<double> get_pmt_wall_position(int pmt_row, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0) {
    return get_pmt_cap_position(pmt_row, 0, pmt_number, starting_pmt_number, angular_offset);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMMainVolume::G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                                 G4LogicalVolume* pMotherLogical, G4bool pMany,
                                 G4int pCopyNo, G4CCMDetectorConstruction* c,
                                 G4bool SodiumSourceOn, G4double SodiumSourceLocation)
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


    // 3cm of steel until we get to vacuum

    // Vacuum jacket
    fVacuum = new G4Tubs("Vacuum", 0*cm, 135*cm, 126*cm, 0*deg, 360*deg);
    fVacuum_log = new G4LogicalVolume(fVacuum, G4Material::GetMaterial("Vacuum"), "Vacuum");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fVacuum_log, "Vacuum", fCryoVessel_log, false, 0, true);

    // 10cm of vacuum until we get to inner steel jacket

    // Inner cryogen
    fInnerJacket = new G4Tubs("InnerJacket", 0*cm, 125*cm, 120*cm, 0*deg, 360*deg);
    fInnerJacket_log = new G4LogicalVolume(fInnerJacket, G4Material::GetMaterial("Steel"), "InnerJacket");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fInnerJacket_log, "InnerJacket", fVacuum_log, false, 0, true);

    // 5cm of steel until we get to liquid argon

    // Argon outside the fiducial volume
    fArgonOuter = new G4Tubs("OuterLiquidArgon", 0*cm, 120*cm, 115*cm, 0*deg, 360*deg);
    fArgonOuter_log = new G4LogicalVolume(fArgonOuter, G4Material::GetMaterial("LAr"), "OuterLiquidArgon");
    new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fArgonOuter_log, "OuterLiquidArgon", fInnerJacket_log, false, 0, true);

    // Radius of inner frame is 1037mm, inner radius is 1036mm
    // Total height of inner frame is 1239.6mm (half height is 619.8mm)

    G4double frame_height = 1239.6*mm;
    G4double frame_half_height = frame_height / 2.0;
    G4double frame_radius = 1037.0*mm;

    // Aluminum frame holding PMTs and instrumentation
    fInnerFrame = new G4Tubs("InnerFrame", 0*cm, frame_radius, frame_half_height, 0*deg, 360*deg);
    fInnerFrame_log= new G4LogicalVolume(fInnerFrame, G4Material::GetMaterial("Alum"), "InnerFrame");
    new G4PVPlacement(0, G4ThreeVector(0,0,0), fInnerFrame_log, "InnerFrame", fArgonOuter_log, false, 0, true);

    G4double frame_thickness = 1.0 * mm;
    G4double ptfe_thickness = 0.5 * mm;
    G4double tpb_thickness = 0.0019 * mm;

    G4double ptfe_half_height = frame_half_height - frame_thickness;
    G4double ptfe_radius = frame_radius - frame_thickness;

    // Teflon reflector foils
    fReflectorFoil = new G4Tubs("PTFEFoil", 0*cm, ptfe_radius, ptfe_half_height, 0*deg, 360*deg);
    fReflectorFoil_log = new G4LogicalVolume(fReflectorFoil, G4Material::GetMaterial("PTFE"), "PTFEFoil");
    fReflectorFoil_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0), fReflectorFoil_log, "PTFEFoil", fInnerFrame_log, false, 0, true);

    G4double tpb_half_height = ptfe_half_height - ptfe_thickness;
    G4double tpb_radius = ptfe_radius - ptfe_thickness;

    fTPBFoil = new G4Tubs("TPBFoil", 0*cm, tpb_radius, tpb_half_height, 0*deg, 360*deg);
    fTPBFoil_log = new G4LogicalVolume(fTPBFoil, G4Material::GetMaterial("TPBFoil"), "TPBFoil");
    fTPBFoil_phys = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fTPBFoil_log, "TPBFoil", fReflectorFoil_log, false, 0, true);

    G4double fiducial_lar_half_height = tpb_half_height - tpb_thickness;
    G4double fiducial_lar_radius = tpb_radius - tpb_thickness;

    // now fiducial LAr!
    fFiducialAr = new G4Tubs("FiducialArgon", 0*cm, fiducial_lar_radius, fiducial_lar_half_height, 0*deg, 360*deg);
    fFiducialAr_log = new G4LogicalVolume(fFiducialAr, G4Material::GetMaterial("LAr"), "FiducialArgon");
    fFiducialAr_phys = new G4PVPlacement(0, G4ThreeVector(0*cm, 0*cm, 0*cm), fFiducialAr_log, "FiducialArgon", fTPBFoil_log, false, 0, true);

    G4double pmt_protrusion_distance = 61.89 * mm;

    // now let's build PMTs using J4SolidMaker
    // note -- we are construction coated PMTs and uncoated PMTs seperately
    //fPMTCoated = J4PMTSolidMaker::Get8inchPMTSolid();
    fPMTCoatedWall = J4PMTSolidMaker::Get8inchPMTWallSolid(fiducial_lar_radius, pmt_protrusion_distance);
    fPMTCoatedCaps = J4PMTSolidMaker::Get8inchPMTCapsSolid(pmt_protrusion_distance);

    fPMTCoatedWall_log = new G4LogicalVolume(fPMTCoatedWall, G4Material::GetMaterial("Glass"), "PMTCoatedWallLog");
    fPMTCoatedCaps_log = new G4LogicalVolume(fPMTCoatedCaps, G4Material::GetMaterial("Glass"), "PMTCoatedCapsLog");

    fPMTUncoatedWall = J4PMTSolidMaker::Get8inchPMTWallSolid(fiducial_lar_radius, pmt_protrusion_distance);
    fPMTUncoatedCaps = J4PMTSolidMaker::Get8inchPMTCapsSolid(pmt_protrusion_distance);

    fPMTUncoatedWall_log = new G4LogicalVolume(fPMTUncoatedWall, G4Material::GetMaterial("Glass"), "PMTUncoatedWallLog");
    fPMTUncoatedCaps_log = new G4LogicalVolume(fPMTUncoatedCaps, G4Material::GetMaterial("Glass"), "PMTUncoatedCapsLog");

    // now get TPB coating -- NOT a daughter volume of coated pmt logical volume
    fTPBCoatingWall = J4PMTSolidMaker::GetTPBCoatingWallSolid(fiducial_lar_radius, pmt_protrusion_distance);
    fTPBCoatingCaps = J4PMTSolidMaker::GetTPBCoatingCapsSolid(pmt_protrusion_distance);
    fTPBCoatingWall_log = new G4LogicalVolume(fTPBCoatingWall, G4Material::GetMaterial("TPBPMT"), "TPBCoatingWallLog");
    fTPBCoatingCaps_log = new G4LogicalVolume(fTPBCoatingCaps, G4Material::GetMaterial("TPBPMT"), "TPBCoatingCapsLog");

    // now let's also build PMT photocathode!
    //fPhotocathCoated = J4PMTSolidMaker::GetPhotcathodeSolid();
    //fPhotocathCoated_log = new G4LogicalVolume(fPhotocathCoated, G4Material::GetMaterial("Vacuum"), "photocathCoated_log");
    //new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), fPhotocathCoated_log, "photocathCoated", fPMTCoated_log, false, 0);

    //fPhotocathUncoated = J4PMTSolidMaker::GetPhotcathodeSolid();
    //fPhotocathUncoated_log = new G4LogicalVolume(fPhotocathUncoated, G4Material::GetMaterial("Vacuum"), "photocathUncoated_log");
    //new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), fPhotocathUncoated_log, "photocathUncoated", fPMTUncoated_log, false, 0);

    // now that we've defined the pmt logical volume, we can get pmt locations using CCMGeometryGenerator logic
    G4int k = 0;
    for (size_t i = 0; i < position_id.size(); i++) {
        std::string position_string = position_id[i];
        size_t r_pos = position_string.find("R");
        size_t n_row_chars = r_pos - 1;
        bool on_caps = n_row_chars > 2;
        int row;
        std::vector<double> position;
        int pmt_number = 0;
        bool coated = true;
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
                coated = false;
        } else {
            int col = std::atoi(position_string.substr(1, r_pos - 1).c_str());
            row = std::atoi(position_string.substr(r_pos + 1, std::string::npos).c_str());
            position = get_pmt_wall_position(row, col);
            pmt_number = col;
            if( (row ==  1 and col % 4 == 1) or
                (row ==  2 and col % 4 == 3) or
                (row ==  4 and col % 4 == 0) or
                (row ==  5 and col % 4 == 2)) {
                coated = false;
            }
        }
        // so we have our pmt info
        // let's make the string key to save it to
        G4String pmt_name = std::to_string(row) + "_" + std::to_string(pmt_number);
        G4ThreeVector pmt_pos;

        // let's make appropriate rotation matrix
        G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
        double distance_to_translate = 10.0; // distance in cm to pull pmts back by

        if (row == 0) {
            // top face pmts
            // rotate so they are pointing downwards
            G4double rotationAngle = M_PI;
            rotationMatrix->rotateX(rotationAngle); // Rotate around the x axis
            position[2] += distance_to_translate; // pulling pmts back a tad
            pmt_pos.setX(position[0]*cm);
            pmt_pos.setY(position[1]*cm);
            pmt_pos.setZ(position[2]*cm);

        }
        if (row == 6){
            // bottom face pmts
            // note -- don't actually need to do a rotation for bottom face pmts
            // but do need to adjust z position slightly
            position[2] -= distance_to_translate;
            pmt_pos.setX(position[0]*cm);
            pmt_pos.setY(position[1]*cm);
            pmt_pos.setZ(position[2]*cm);
        }
        if (row > 0 and row < 6){
            // side pmts
            // we want to calculate facing direction essentially
            G4double facing_radius = std::pow(position[0], 2) + std::pow(position[1], 2);
            G4double facing_r_x = position[0] / facing_radius;
            G4double facing_r_y = position[1] / facing_radius;
            G4double facing_dir_norm_factor = std::sqrt(std::pow(facing_r_x, 2) + std::pow(facing_r_y, 2));
            G4double facing_dir_x = facing_r_x / facing_dir_norm_factor;
            G4double facing_dir_y = facing_r_y / facing_dir_norm_factor;
            G4double facing_dir_z = 0.0;

            G4double phi = std::atan2(facing_dir_y, facing_dir_x);
            G4double theta = std::acos(facing_dir_z);

            rotationMatrix->rotateZ(-phi); // rotate rotate
            rotationMatrix->rotateY(theta); // rotate rotate

            // also need to pull the pmts back a touch
            // we're going to change the magnitude of facing direction vector to equal distance_to_translate
            G4double translation_x = facing_dir_x * distance_to_translate / 2 + position[0];
            G4double translation_y = facing_dir_y * distance_to_translate / 2 + position[1];
            pmt_pos.setX(translation_x*cm);
            pmt_pos.setY(translation_y*cm);
            pmt_pos.setZ(position[2]*cm);

        }
        if(coated) {
            // place coated pmts
            G4String descriptive_name = "CoatedPMT_" + pmt_name;
            G4String tpb_descriptive_name = "TPBCoating_" + pmt_name;
            if(row > 0 and row < 6) {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTCoatedWall_log, descriptive_name, fFiducialAr_log, false, k);
                fTPBPMT_phys = new G4PVPlacement(rotationMatrix, pmt_pos, fTPBCoatingWall_log, tpb_descriptive_name, fFiducialAr_log, false, k);
            } else {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTCoatedCaps_log, descriptive_name, fFiducialAr_log, false, k);
                fTPBPMT_phys = new G4PVPlacement(rotationMatrix, pmt_pos, fTPBCoatingCaps_log, tpb_descriptive_name, fFiducialAr_log, false, k);
            }
        } else {
            // place uncoated pmts
            G4String descriptive_name = "UncoatedPMT_" + pmt_name;
            if(row > 0 and row < 6) {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTUncoatedWall_log, descriptive_name, fFiducialAr_log, false, k);
            } else {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTUncoatedCaps_log, descriptive_name, fFiducialAr_log, false, k);
            }
        }

        // now save positions and increment counter
        fPMTPositions.push_back(G4ThreeVector(position[0], position[1], position[2]));
        ++k;
    }

    // now let's make our sodium source pellet + rod
    // i'm just doing the thin rod for now -- TODO : add all rod specifics (main part, nuts, thin extension)
    G4double rod_inner_radius = 0.4*cm;
    G4double rod_outer_radius = 0.5*cm;
    G4double rod_height = 150.0*cm;
    fSourceRod = new G4Tubs("SourceRod", rod_inner_radius, rod_outer_radius, rod_height/2, 0, 360*deg);
    fSourceRod_log = new G4LogicalVolume(fSourceRod,  G4Material::GetMaterial("Steel"), "fSourceRodLog");

    // now make sodium pellet -- really guessing on these measurements
    G4double pellet_radius = 4.0*mm;
    G4double pellet_height = 3.0*mm;
    fSodiumSourcePellet = new G4Tubs("SodiumSourcePellet", 0, pellet_radius, pellet_height/2.0, 0, 360*deg);
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* fNa = nistManager->FindOrBuildMaterial("G4_Na");
    fSodiumSourcePellet_log  = new G4LogicalVolume(fSodiumSourcePellet, fNa, "SodiumSourcePelletLog");

    if (SodiumSourceOn){
        G4ThreeVector rodPosition(0.0*cm, 0.0*cm, SodiumSourceLocation + rod_height/2);
        new G4PVPlacement(nullptr, rodPosition, fSourceRod_log, "SourceRod", fFiducialAr_log, false, 0);

        // and now put the sodium pellet at the end of the rod (inset 1/4cm)
        G4double inset = 0.25 * cm;
        G4ThreeVector pelletPosition(0.0*cm, 0.0*cm, SodiumSourceLocation + pellet_height/2.0 + inset);
        new G4PVPlacement(nullptr, pelletPosition, fSodiumSourcePellet_log, "SodiumSourcePellet", fFiducialAr_log, false, 0);
    }

    VisAttributes();
    SurfaceProperties();
    SetLogicalVolume(fFiducialAr_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::VisAttributes()
{
    auto argon_va = new G4VisAttributes(G4Colour(0., 0., 1.0)); //blue
    fFiducialAr_log->SetVisAttributes(argon_va);

    // let's make other layers of detector invisible
    fCryoVessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fVacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fInnerJacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fArgonOuter_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fInnerFrame_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fTPBFoil_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTCoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTUncoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBCoating_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPhotocathCoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPhotocathUncoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    //auto pmt_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8)); //grey idk
    //pmt_va->SetForceSolid(true);
    //fPMTCoated_log->SetVisAttributes(pmt_va);
    //fPMTUncoated_log->SetVisAttributes(pmt_va);
    //
    //auto tpb_coating_va = new G4VisAttributes(G4Colour(0., 1., 0.)); //green
    //tpb_coating_va->SetForceSolid(true);
    //fTPBCoating_log->SetVisAttributes(tpb_coating_va);

    //auto photocath_va = new G4VisAttributes(G4Colour(0., 0., 1.)); // blue
    //photocath_va->SetForceSolid(true);
    //fPhotocathCoated_log->SetVisAttributes(photocath_va);
    //fPhotocathUncoated_log->SetVisAttributes(photocath_va);

    // make our source rod green
    auto source_rod_va = new G4VisAttributes(G4Colour(0., 1., 0.)); //green
    source_rod_va->SetForceSolid(true);
    fSourceRod_log->SetVisAttributes(source_rod_va);
    //fSourceRod_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    // and make sodium pellet blue
    auto sodium_pellet_va = new G4VisAttributes(G4Colour(0., 0., 1.)); // blue
    sodium_pellet_va->SetForceSolid(true);
    fSodiumSourcePellet_log->SetVisAttributes(sodium_pellet_va);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::SurfaceProperties()
{
    // now it's time to define optical surface properties

    // TPB foil optical surface
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

    G4OpticalSurface *TPBFoilOpticalSurface = new G4OpticalSurface("TPBFoilOpticalSurface");

    TPBFoilOpticalSurface->SetModel(unified);
    TPBFoilOpticalSurface->SetType(dielectric_dielectric);
    TPBFoilOpticalSurface->SetFinish(ground);
    TPBFoilOpticalSurface->SetSigmaAlpha(0.05);

    G4MaterialPropertiesTable *TPBFoil_mt = new G4MaterialPropertiesTable();
    TPBFoil_mt->AddProperty("REFLECTIVITY", TPBEnergy, TPBfoilOSReflect);
    TPBFoil_mt->AddProperty("TRANSMITTANCE", TPBEnergy, TPBfoilOSTransmit);
    TPBFoil_mt->AddProperty("EFFICIENCY", TPBEnergy, TPBfoilOSEff);
    TPBFoilOpticalSurface->SetMaterialPropertiesTable(TPBFoil_mt);

    // note -- usign logical skin surface, might want to use border surface but since using same properties for all TPB on walls, doesnt seem necessary
    new G4LogicalSkinSurface("TPBFoils_Surface", fTPBFoil_log, TPBFoilOpticalSurface);
    new G4LogicalSkinSurface("TPBCoatingWall_Surface", fTPBCoatingWall_log, TPBFoilOpticalSurface);
    new G4LogicalSkinSurface("TPBCoatingCaps_Surface", fTPBCoatingCaps_log, TPBFoilOpticalSurface);

    //// create logical border surfaces for TPB on walls of detector and on PMTs
	//new G4LogicalBorderSurface("TPBFoils_SurfaceForward", fFiducialAr_phys, fTPBFoil_phys, TPBFoilOpticalSurface);
	//new G4LogicalBorderSurface("TPBFoils_SurfaceBackward", fTPBFoil_phys, fFiducialAr_phys, TPBFoilOpticalSurface);
	//
    //// same for TPB coating -- note pretty sure coating + foils have same optical surface properties
    //new G4LogicalBorderSurface("TPBPMT_SurfaceForward", fFiducialAr_phys, fTPBPMT_phys, TPBFoilOpticalSurface);
	//new G4LogicalBorderSurface("TPBPMT_SurfaceBackward", fTPBPMT_phys, fFiducialAr_phys, TPBFoilOpticalSurface);

    // finally define optical surface for PTFE reflector foils
    G4OpticalSurface *PTFEFoilOpticalSurface = new G4OpticalSurface("PTFEFoilOpticalSurface");

    PTFEFoilOpticalSurface->SetModel(glisur); //Optical model
    PTFEFoilOpticalSurface->SetType(dielectric_metal);
    PTFEFoilOpticalSurface->SetFinish(ground);
    PTFEFoilOpticalSurface->SetSigmaAlpha(1.0);

    G4MaterialPropertiesTable *reflfoilMPT = new G4MaterialPropertiesTable();
    std::vector<G4double> PTFEFoilOpticalSurfaceEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV,  2.38*eV,
                                                          2.48*eV,  2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,  3.124*eV, 3.457*eV,
                                                          3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                                          7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
    G4double uvRf = 0.10; //change uvreflection
    G4double vsRf = 0.95; //change visible reflectivity
    //enter the uv and visible reflectivities defined above.
    std::vector<G4double> PTFEFoilOpticalSurfaceReflect = { vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, vsRf,
                                                            vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, vsRf,
                                                            vsRf, vsRf, vsRf, vsRf, uvRf, uvRf, uvRf,
                                                            uvRf, uvRf, uvRf, uvRf};
    std::vector<G4double> PTFEFoilOpticalSurfaceTransmit = {.02, .02, .02, .02, .02, .02, .02, .02, .02, .02,
                                                            .02, .02, .02, .02, .02, .02, .02, .02, .02, .02,
                                                            .02, .02, .02, .02, .02}; //minimal transmittance through foil
    std::vector<G4double> PTFEFoilOpticalSurfaceEff = { 0., 0., 0., 0., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                        0.0, 0.0, 0.0, 0.0, 0., 0., 0., 0., 0., 0.,
                                                        0., 0., 0., 0., 0.};
    reflfoilMPT->AddProperty("TRANSMITTANCE", PTFEFoilOpticalSurfaceEnergy, PTFEFoilOpticalSurfaceTransmit);
    reflfoilMPT->AddProperty("REFLECTIVITY", PTFEFoilOpticalSurfaceEnergy, PTFEFoilOpticalSurfaceReflect);
    reflfoilMPT->AddProperty("EFFICIENCY", PTFEFoilOpticalSurfaceEnergy, PTFEFoilOpticalSurfaceEff);
    PTFEFoilOpticalSurface->SetMaterialPropertiesTable(reflfoilMPT);

    // again using skin, keeping border surface example
    new G4LogicalSkinSurface("PTFE_Surface", fReflectorFoil_log, PTFEFoilOpticalSurface);

	//// now add logical borders for PTFE foil --> TPB foil and LAr
    //new G4LogicalBorderSurface("TPBFoilPTFE_SurfaceForward", fTPBFoil_phys, fReflectorFoil_phys, PTFEFoilOpticalSurface);
 	//new G4LogicalBorderSurface("TPBFoilPTFE_SurfaceBackward", fReflectorFoil_phys, fTPBFoil_phys, PTFEFoilOpticalSurface);

}

