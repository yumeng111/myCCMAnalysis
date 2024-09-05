
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
                                        "C406R0",
                                        "C205R6",
                                        "C408R6",
                                        "C306R6",
                                        "C305R6"};

// Z positions of the rows of pmts
std::map<size_t, double> pmt_region_z_positions = {
    {-2, 75.00},
    {-1, 65.00},
    {0,  58.00},
    {1,  46.00},
    {2,  23.00},
    {3,   0.00},
    {4, -23.00},
    {5, -46.00},
    {6, -58.00},
    {7, -65.00},
    {8, -75.00},
};

// The radii of rings of pmts
std::map<size_t, double> ring_radii = {
    {0,  96.0}, // Outer wall
    {1,  24.0}, // Inner ring on top/bottom
    {2,  45.0}, //
    {3,  66.0}, //
    {4,  87.0}, // Outer ring on top/bottom
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

std::vector<double> get_pmt_cap_position(int pmt_row, int ring_number, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0, double z=0.0) {
    if(z == 0.0) {
        z = pmt_region_z_positions[pmt_row];
    }
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

std::vector<double> get_pmt_wall_position(int pmt_row, int pmt_number, double starting_pmt_number = 1, double angular_offset = 0.0, double radius=0.0) {
    int ring_number = 0;
    double z = pmt_region_z_positions[pmt_row];
    if(radius <= 0.0) {
        radius = ring_radii[ring_number];
    }
    std::pair<double, double> xy = xy_position_from_cylinder_pmt_number(
            pmt_number,
            radius,
            starting_pmt_number,
            ring_pmt_pos_count[ring_number],
            angular_offset);
    std::vector<double> pos;
    pos.push_back(xy.first);
    pos.push_back(xy.second);
    pos.push_back(z);
    return pos;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMMainVolume::G4CCMMainVolume(G4RotationMatrix* pRot, const G4ThreeVector& tlate,
                                 G4LogicalVolume* pMotherLogical, G4bool pMany,
                                 G4int pCopyNo, G4CCMDetectorConstruction* c,
                                 G4bool SourceRodIn, G4double SourceRodLocation, G4bool CobaltSourceRun, G4bool SodiumSourceRun)
  // Pass info to the G4PVPlacement constructor
  : G4PVPlacement(pRot, tlate,
            new G4LogicalVolume(
                new G4Tubs("Cryogen", 0*cm, 138*cm, 131*cm, 0*deg,360*deg),
                G4Material::GetMaterial("Steel"), "Cryogen"),
            "Cryogen", pMotherLogical, false, 0, true)
  , fConstructor(c) {

    // now let's build our detector
    // Outer cryogen
    fCryoVessel_log = GetLogicalVolume();

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
    G4double frame_thickness = 1.0 * mm;
    G4double ptfe_thickness = 0.5 * mm;
    G4double tpb_thickness = 0.00278035 * mm; // from Vincent Basque's thesis https://pure.manchester.ac.uk/ws/portalfiles/portal/205622566/FULL_TEXT.PDF

    G4double frame_height = 1239.6*mm + (frame_thickness + ptfe_thickness + tpb_thickness) * 2.0;
    G4double frame_half_height = frame_height / 2.0;
    G4double frame_radius = 1037.0*mm - 3.15*mm + (frame_thickness + ptfe_thickness + tpb_thickness); // reducing radius by ~1/2cm to account for flexing of PTFE sheets

    // Aluminum frame holding PMTs and instrumentation
    fInnerFrame = new G4Tubs("InnerFrame", 0*cm, frame_radius, frame_half_height, 0*deg, 360*deg);
    fInnerFrame_log= new G4LogicalVolume(fInnerFrame, G4Material::GetMaterial("Alum"), "InnerFrame");
    new G4PVPlacement(0, G4ThreeVector(0,0,0), fInnerFrame_log, "InnerFrame", fArgonOuter_log, false, 0, true);

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
    fPMTCoatedWall = J4PMTSolidMaker::Get8inchPMTWallSolid(fiducial_lar_radius, pmt_protrusion_distance);
    fPMTCoatedCaps = J4PMTSolidMaker::Get8inchPMTCapsSolid(pmt_protrusion_distance);
    fPMTCoatedWall_log = new G4LogicalVolume(fPMTCoatedWall, G4Material::GetMaterial("Glass"), "PMTCoatedWallLog");
    fPMTCoatedCaps_log = new G4LogicalVolume(fPMTCoatedCaps, G4Material::GetMaterial("Glass"), "PMTCoatedCapsLog");

    fPMTUncoatedWall = J4PMTSolidMaker::Get8inchPMTWallSolid(fiducial_lar_radius, pmt_protrusion_distance);
    fPMTUncoatedCaps = J4PMTSolidMaker::Get8inchPMTCapsSolid(pmt_protrusion_distance);
    fPMTUncoatedWall_log = new G4LogicalVolume(fPMTUncoatedWall, G4Material::GetMaterial("Glass"), "PMTUncoatedWallLog");
    fPMTUncoatedCaps_log = new G4LogicalVolume(fPMTUncoatedCaps, G4Material::GetMaterial("Glass"), "PMTUncoatedCapsLog");

    // now let's build our bridle + frills that go around the PMTs
    G4double bridle_width = 12.0 * mm;
    G4double bridle_radius = 209.55 / 2.0 * mm;
    G4double frill_width = 0.813 * mm;
    G4double frill_radius = 228.6 / 2.0 * mm;
    fBridleWall = J4PMTSolidMaker::GetBridleWall(bridle_radius, bridle_width, tpb_radius, pmt_protrusion_distance);
    fBridleCaps = J4PMTSolidMaker::GetBridleCaps(bridle_radius, bridle_width, pmt_protrusion_distance);
    fFrillWall = J4PMTSolidMaker::GetFrillWall(bridle_radius, frill_radius,  frill_width, tpb_radius);
    fFrillCaps = J4PMTSolidMaker::GetFrillCaps(bridle_radius, frill_radius, frill_width);
    fBridleWall_log = new G4LogicalVolume(fBridleWall, G4Material::GetMaterial("BlackPlastic"), "BridleWallLog");
    fBridleCaps_log = new G4LogicalVolume(fBridleCaps, G4Material::GetMaterial("BlackPlastic"), "BridleCapsLog");
    fFrillWall_log = new G4LogicalVolume(fFrillWall, G4Material::GetMaterial("Plastic"), "FrillWallLog");
    fFrillCaps_log = new G4LogicalVolume(fFrillCaps, G4Material::GetMaterial("Plastic"), "FrillCapsLog");

    // now get TPB coating
    G4double tpb_protrusion = pmt_protrusion_distance - bridle_width;
    fTPBCoatingWall = J4PMTSolidMaker::GetTPBCoatingWallSolid(fiducial_lar_radius, pmt_protrusion_distance);
    fTPBCoatingCaps = J4PMTSolidMaker::GetTPBCoatingCapsSolid(pmt_protrusion_distance);
    fTPBCoatingWall_log = new G4LogicalVolume(fTPBCoatingWall, G4Material::GetMaterial("TPBPMT"), "TPBCoatingWallLog");
    fTPBCoatingCaps_log = new G4LogicalVolume(fTPBCoatingCaps, G4Material::GetMaterial("TPBPMT"), "TPBCoatingCapsLog");

    // let's build the shiny guy who is at C406R0
    // we are placing it in our loop over PMT position IDs!
    G4double shiny_radius = 205.0 / 2.0 * mm; // making the assumption that it is the same radius as the bridle
    G4double shiny_half_height = 0.125 * mm;
    fShinyC406R0 = new G4Tubs("ShinyC406R0", 0*cm, shiny_radius, shiny_half_height, 0*deg, 360*deg);
    fShinyC406R0_log= new G4LogicalVolume(fShinyC406R0, G4Material::GetMaterial("Alum"), "ShinyC406R0");

    // speaking of shiny, let's build our shiny reflective circles that go on top and bottom of detector
    // documentation (https://docdb.lns.mit.edu/cgi-bin/captainmills/ShowDocument?docid=567) says it's 5.25in radius
    G4double top_shiny_radius = 13.335 * cm;
    fShinyTop = new G4Tubs("ShinyTop", 0*cm, top_shiny_radius, shiny_half_height, 0*deg, 360*deg);
    fShinyTop_log= new G4LogicalVolume(fShinyTop, G4Material::GetMaterial("Alum"), "ShinyTop");
    new G4PVPlacement(0, G4ThreeVector(0, 0, fiducial_lar_half_height - shiny_half_height), fShinyTop_log, "ShinyTop", fFiducialAr_log, false, 0);
    fShinyBottom = new G4Tubs("ShinyBottom", 0*cm, top_shiny_radius, shiny_half_height, 0*deg, 360*deg);
    fShinyBottom_log= new G4LogicalVolume(fShinyBottom, G4Material::GetMaterial("Alum"), "ShinyBottom");
    new G4PVPlacement(0, G4ThreeVector(0, 0, - fiducial_lar_half_height + shiny_half_height), fShinyBottom_log, "ShinyBottom", fFiducialAr_log, false, 0);

    double pmt_radius_cm = (fiducial_lar_radius + (J4PMTSolidMaker::Get8inchPMTRadius() - pmt_protrusion_distance)) / cm;
    double pmt_height_cm = (fiducial_lar_half_height + (J4PMTSolidMaker::Get8inchPMTRadius() - pmt_protrusion_distance)) / cm;
    double bridle_height_cm = (fiducial_lar_half_height - (bridle_width / 2.0)) / cm;
    double frill_height_cm = (fiducial_lar_half_height - (frill_width / 2.0)) / cm;
    double bridle_radius_cm = fiducial_lar_radius / cm;
    double frill_radius_cm =  fiducial_lar_radius / cm;

    // now that we've defined the pmt logical volume, we can get pmt locations using CCMGeometryGenerator logic
    G4int k = 0;
    for (size_t i = 0; i < position_id.size(); i++) {
        std::string position_string = position_id[i];
        size_t r_pos = position_string.find("R");
        size_t n_row_chars = r_pos - 1;
        bool on_caps = n_row_chars > 2;
        int row;
        std::vector<double> position;
        std::vector<double> bridle_position;
        std::vector<double> frill_position;
        int pmt_number = 0;
        bool coated = true;
        if(on_caps) {
            int col = std::atoi(position_string.substr(2, n_row_chars - 1).c_str());
            int ring = std::atoi(position_string.substr(1, 1).c_str());
            row = std::atoi(position_string.substr(r_pos + 1, std::string::npos).c_str());
            if(row == 0) {
                position = get_pmt_cap_position(row, ring, col, 1, 0.0, pmt_height_cm);
                bridle_position = get_pmt_cap_position(row, ring, col, 1, 0.0, bridle_height_cm);
                frill_position = get_pmt_cap_position(row, ring, col, 1, 0.0, frill_height_cm);
            } else if(row == 6) {
                position = get_pmt_cap_position(row, ring, col, 1, 0.0, -pmt_height_cm);
                bridle_position = get_pmt_cap_position(row, ring, col, 1, 0.0, -bridle_height_cm);
                frill_position = get_pmt_cap_position(row, ring, col, 1, 0.0, -frill_height_cm);
            } else {
                throw std::runtime_error("Bad cap row");
            }
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
            position = get_pmt_wall_position(row, col, 1, 0.0, pmt_radius_cm);
            bridle_position = get_pmt_wall_position(row, col, 1, 0.0, bridle_radius_cm);
            frill_position = get_pmt_wall_position(row, col, 1, 0.0, frill_radius_cm);
            pmt_number = col;
            if( (row ==  1 and col % 4 == 1) or
                (row ==  2 and col % 4 == 3) or
                (row ==  4 and col % 4 == 0) or
                (row ==  5 and col % 4 == 2)) {
                coated = false;
            }
        }
        if (position_string == "C406R0"){
            // now place our shiny circle!
            new G4PVPlacement(0, G4ThreeVector(position[0]*cm, position[1]*cm, fiducial_lar_half_height - shiny_half_height), fShinyC406R0_log, "ShinyC406R0", fFiducialAr_log, false, 0);
            continue;
        }
        // so we have our pmt info
        // let's make the string key to save it to
        G4String pmt_name = std::to_string(row) + "_" + std::to_string(pmt_number);
        G4ThreeVector pmt_pos;
        pmt_pos.setX(position[0]*cm);
        pmt_pos.setY(position[1]*cm);
        pmt_pos.setZ(position[2]*cm);

        G4ThreeVector bridle_pos;
        bridle_pos.setX(bridle_position[0]*cm);
        bridle_pos.setY(bridle_position[1]*cm);
        bridle_pos.setZ(bridle_position[2]*cm);

        G4ThreeVector frill_pos;
        frill_pos.setX(frill_position[0]*cm);
        frill_pos.setY(frill_position[1]*cm);
        frill_pos.setZ(frill_position[2]*cm);

        // let's make appropriate rotation matrix
        G4RotationMatrix* rotationMatrix = new G4RotationMatrix();

        if (row == 0) {
            // top face pmts
            // rotate so they are pointing downwards
            G4double rotationAngle = M_PI;
            rotationMatrix->rotateX(rotationAngle); // Rotate around the x axis

        }
        // note -- don't actually need to do a rotation for bottom face pmts -- skipping for now
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
        }
        if(coated) {
            // place coated pmts
            G4String descriptive_name = "CoatedPMT_" + pmt_name;
            G4String tpb_descriptive_name = "TPBCoating_" + pmt_name;
            G4String bridle_descriptive_name = "Bridle_" + pmt_name;
            G4String frill_descriptive_name = "Frill_" + pmt_name;

            if(row > 0 and row < 6) {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTCoatedWall_log, descriptive_name, fFiducialAr_log, false, k);
                fTPBPMT_phys = new G4PVPlacement(rotationMatrix, pmt_pos, fTPBCoatingWall_log, tpb_descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, bridle_pos, fBridleWall_log, bridle_descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, frill_pos, fFrillWall_log, frill_descriptive_name, fFiducialAr_log, false, k);
            } else {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTCoatedCaps_log, descriptive_name, fFiducialAr_log, false, k);
                fTPBPMT_phys = new G4PVPlacement(rotationMatrix, pmt_pos, fTPBCoatingCaps_log, tpb_descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, bridle_pos, fBridleCaps_log, bridle_descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, frill_pos, fFrillCaps_log, frill_descriptive_name, fFiducialAr_log, false, k);
            }
        } else {
            // place uncoated pmts
            G4String descriptive_name = "UncoatedPMT_" + pmt_name;
            G4String bridle_descriptive_name = "Bridle_" + pmt_name;
            G4String frill_descriptive_name = "Frill_" + pmt_name;

            if(row > 0 and row < 6) {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTUncoatedWall_log, descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, bridle_pos, fBridleWall_log, bridle_descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, frill_pos, fFrillWall_log, frill_descriptive_name, fFiducialAr_log, false, k);
            } else {
                new G4PVPlacement(rotationMatrix, pmt_pos, fPMTUncoatedCaps_log, descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, bridle_pos, fBridleCaps_log, bridle_descriptive_name, fFiducialAr_log, false, k);
                new G4PVPlacement(rotationMatrix, frill_pos, fFrillCaps_log, frill_descriptive_name, fFiducialAr_log, false, k);
            }
        }

        // now save positions and increment counter
        fPMTPositions.push_back(G4ThreeVector(position[0], position[1], position[2]));
        ++k;
    }

    if (SourceRodIn){
        // now let's make our source rod
        // i'm just doing the bayonet for now since it is >1m tall, it should only matter for sodium runs above -50cm
        G4double rod_inner_radius = 0.0*cm;
        G4double rod_outer_radius = 6.31 * mm;
        G4double rod_height = fiducial_lar_half_height - SourceRodLocation - (shiny_half_height * 2.0);
        fSourceRod = new G4Tubs("SourceRod", rod_inner_radius, rod_outer_radius, rod_height/2, 0, 360*deg);
        fSourceRod_log = new G4LogicalVolume(fSourceRod,  G4Material::GetMaterial("Steel"), "fSourceRodLog");
        G4ThreeVector rodPosition(0.0*cm, 0.0*cm, SourceRodLocation + rod_height/2);
        new G4PVPlacement(nullptr, rodPosition, fSourceRod_log, "SourceRod", fFiducialAr_log, false, 0);

        if (CobaltSourceRun or SodiumSourceRun){
            // now make source pellet -- really guessing on these measurements
            G4double pellet_radius = 4.0*mm;
            G4double pellet_height = 3.0*mm;
            fSourcePellet = new G4Tubs("SourcePellet", 0, pellet_radius, pellet_height/2.0, 0, 360*deg);
            G4NistManager* nistManager = G4NistManager::Instance();
            G4Material* fSource = nistManager->FindOrBuildMaterial("G4_Na");
            if (CobaltSourceRun){
                fSource = nistManager->FindOrBuildMaterial("G4_Co");
            }
            fSourcePellet_log  = new G4LogicalVolume(fSourcePellet, fSource, "SourcePelletLog");

            // and now put the source pellet at the end of the rod (inset 1/4cm)
            G4double inset = 0.25 * cm;
            G4ThreeVector pelletPosition(0.0*cm, 0.0*cm, SourceRodLocation + pellet_height/2.0 + inset);
            new G4PVPlacement(nullptr, pelletPosition, fSourcePellet_log, "SourcePellet", fFiducialAr_log, false, 0);
        }
    }

    VisAttributes(SourceRodIn);
    SurfaceProperties();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::VisAttributes(G4bool SourceRodIn)
{
    // let's make other layers of detector invisible
    auto lilac = new G4VisAttributes(G4Colour(153/255., 153/255., 153/255.));
    //lilac->SetForceSolid(true);
    auto salmon = new G4VisAttributes(G4Colour(255/255., 128/255., 128/255.));
    //salmon->SetForceSolid(true);
    auto dark_blue = new G4VisAttributes(G4Colour(0., 0., 128/255.));
    //dark_blue->SetForceSolid(true);
    auto green = new G4VisAttributes(G4Colour(51/255., 153/255., 102/255.));
    //green->SetForceSolid(true);
    auto orange = new G4VisAttributes(G4Colour(255/255., 153/255., 0.));
    //orange->SetForceSolid(true);
    auto teal = new G4VisAttributes(G4Colour(0., 128/255., 128/255.));
    //teal->SetForceSolid(true);
    auto pink = new G4VisAttributes(G4Colour(255/255., 0., 255/255.));
    //pink->SetForceSolid(true);

    fCryoVessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fVacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fInnerJacket_log->SetVisAttributes(dark_blue);
    fArgonOuter_log->SetVisAttributes(green);
    fInnerFrame_log->SetVisAttributes(orange);
    fTPBFoil_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    fFiducialAr_log->SetVisAttributes(pink);

    //fCryoVessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fVacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fInnerJacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fArgonOuter_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fInnerFrame_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBFoil_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTCoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTUncoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBCoating_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPhotocathCoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPhotocathUncoated_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto pmt_va = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8)); //grey idk
    pmt_va->SetForceSolid(true);
    fPMTCoatedWall_log->SetVisAttributes(pmt_va);
    fPMTCoatedCaps_log->SetVisAttributes(pmt_va);
    fPMTUncoatedWall_log->SetVisAttributes(pmt_va);
    fPMTUncoatedCaps_log->SetVisAttributes(pmt_va);
    //fPMTCoatedWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTCoatedCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTUncoatedWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fPMTUncoatedCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto tpb_coating_va = new G4VisAttributes(G4Colour(0., 1., 0.)); //green
    tpb_coating_va->SetForceSolid(true);
    fTPBCoatingWall_log->SetVisAttributes(tpb_coating_va);
    fTPBCoatingCaps_log->SetVisAttributes(tpb_coating_va);
    //fTPBCoatingWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fTPBCoatingCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    salmon->SetForceSolid(true);
    teal->SetForceSolid(true);
    fBridleWall_log->SetVisAttributes(salmon);
    fBridleCaps_log->SetVisAttributes(salmon);
    fFrillWall_log->SetVisAttributes(teal);
    fFrillCaps_log->SetVisAttributes(teal);
    //fBridleWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fBridleCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fFrillWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    //fFrillCaps_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    // add shiny guy!
    auto red = new G4VisAttributes(G4Colour(1., 0., 0.));
    red->SetForceSolid(true);
    fShinyC406R0_log->SetVisAttributes(red);
    fShinyTop_log->SetVisAttributes(red);
    fShinyBottom_log->SetVisAttributes(red);

    // make our source rod green
    if (SourceRodIn)
        fSourceRod_log->SetVisAttributes(tpb_coating_va);
    //fSourceRod_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    //// and make sodium pellet blue
    //auto sodium_pellet_va = new G4VisAttributes(G4Colour(0., 0., 1.)); // blue
    //sodium_pellet_va->SetForceSolid(true);
    //fSodiumSourcePellet_log->SetVisAttributes(sodium_pellet_va);
    //fSourcePellet_log->SetVisAttributes(G4VisAttributes::GetInvisible());


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMMainVolume::SurfaceProperties()
{
    // now it's time to define optical surface properties
    G4OpticalSurface *ReflectorOpticalSurface = new G4OpticalSurface("ReflectorOpticalSurface");

    ReflectorOpticalSurface->SetModel(unified);
    ReflectorOpticalSurface->SetType(dielectric_metal);
    ReflectorOpticalSurface->SetFinish(polished);

    // let's add reflection properties as a function of energy
    std::vector<G4double> MylarReflectionEnergy = {1.239835823924519*eV, 2.1758440843642446*eV, 2.2155183016259334*eV, 2.2540941739092317*eV, 2.2965084198409396*eV,
            2.3398220945804553*eV, 2.3853127919447585*eV, 2.432339257160777*eV, 2.4816594267685304*eV, 2.5316156616832233*eV, 2.584054978035275*eV, 2.639307340578113*eV,
            2.6963833000212216*eV, 2.7566402380421113*eV, 2.820482241642379*eV, 2.8840814650497304*eV, 2.9538999741438468*eV, 3.0256914075766623*eV, 3.104173784836434*eV,
            3.1423952308218324*eV, 3.1815179074033235*eV, 3.2230571642902683*eV, 3.2634554851366535*eV, 3.306493319706959*eV, 3.352069231820206*eV, 3.445795378257322*eV,
            3.541505417163459*eV, 3.647364059603557*eV, 3.7597652297046125*eV, 3.876624049427723*eV, 4.000980170596548*eV, 4.135062961692673*eV, 4.276792392741077*eV,
            4.430351962448966*eV, 8.265572159496793*eV, 12.39835823924519*eV};

    std::vector<G4double> MylarReflection = {99.3226740868548, 99.3226740868548, 99.7628366612554, 99.7214862392674, 98.91393146613707, 98.91530448108364, 98.96965285731727, 98.76966534546983,
            98.1497114120479, 98.72409290812726, 98.32154569778484, 98.03407282162186, 97.82027684556721, 98.53168216346992, 97.48482826997746, 96.10633382322926, 96.51396948976577,
            96.65116783853793, 95.91151864322745, 95.30515359961677, 91.25574897695003, 81.90244998989932, 41.8931179140657, 24.101847317428096, 12.26957594631375, 12.06518399527252,
            12.58223308670786, 12.888602883272128, 14.106365509586524, 14.817742868938694, 15.529120228290893, 14.721537496610026, 13.103799846600765, 12.245584179933772, 0.0, 0.0};

    G4MaterialPropertiesTable *ReflectiveFoilMPT = new G4MaterialPropertiesTable();
    ReflectiveFoilMPT->AddProperty("REFLECTIVITY", MylarReflectionEnergy, MylarReflection);
    ReflectorOpticalSurface->SetMaterialPropertiesTable(ReflectiveFoilMPT);

    // TPB foil optical surface
    std::vector<G4double> TPBEnergy = { 0.602*eV/*(2066nm)*/, 0.689*eV/*(1799nm)*/, 1.030*eV/*(1204nm)*/, 1.926*eV/*(644nm)*/, 2.138*eV/* (580nm)*/,
                                        2.250*eV/*(551nm)*/,  2.380*eV/*(521nm)*/,  2.480*eV/*(500nm)*/,  2.583*eV/*(480nm)*/, 2.800*eV/*(443nm)*/,
                                        2.880*eV/*(431nm)*/,  2.980*eV/*(416nm)*/,  3.124*eV/*(397nm)*/,  3.457*eV/*(359nm)*/, 3.643*eV/*(341nm)*/,
                                        3.812*eV/*(325nm)*/,  4.086*eV/*(304nm)*/,  4.511*eV/*(275nm)*/,  5.166*eV/*(240nm)*/, 5.821*eV/*(213nm)*/,
                                        6.526*eV/*(190nm)*/,  8.266*eV/*(150nm)*/,  9.686*eV/*(128nm)*/,  11.27*eV/*(110nm)*/, 12.60*eV/*(98nm)*/  };

    std::vector<G4double> TPBfoilOSTransmit = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                                               1., 1., 1., 1., 1., 1., 1., 1.}; // set to 1 and have all absorption in bulk

    std::vector<G4double> TPBfoilOSEff = {0., 0., 0., 0., 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                          1.0, 1.0, 1.0, 1.0, 1., 1., 1., 1., 1., 1.,
                                          1., 1., 1., 1., 1.};

    // define optical surface for PTFE reflector foils
    // note -- having PTFE describe the TPB reflections as well -- using Lambertian distribution (diffuse)
    //G4OpticalSurface *PTFEFoilOpticalSurface = new G4OpticalSurface("PTFEFoilOpticalSurface");

    //PTFEFoilOpticalSurface->SetModel(unified);
    ////PTFEFoilOpticalSurface->SetType(dielectric_dielectric);
    //PTFEFoilOpticalSurface->SetType(dielectric_metal);
    ////PTFEFoilOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections
    //PTFEFoilOpticalSurface->SetFinish(polished); // 100% specular reflections

    //G4MaterialPropertiesTable *reflfoilMPT = new G4MaterialPropertiesTable();
    //std::vector<G4double> PTFEFoilOpticalSurfaceEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV,  2.38*eV,
    //                                                      2.48*eV,  2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,  3.124*eV, 3.457*eV,
    //                                                      3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
    //                                                      7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };
    //G4double uvRf = 0.10; //change uvreflection
    //G4double vsRf = 0.95; //change visible reflectivity
    ////enter the uv and visible reflectivities defined above.
    //std::vector<G4double> PTFEFoilOpticalSurfaceReflect = { vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, vsRf,
    //                                                        vsRf, vsRf, vsRf, vsRf, vsRf, vsRf, vsRf,
    //                                                        vsRf, vsRf, vsRf, vsRf, uvRf, uvRf, uvRf,
    //                                                        uvRf, uvRf, uvRf, uvRf};
    //reflfoilMPT->AddProperty("REFLECTIVITY", PTFEFoilOpticalSurfaceEnergy, PTFEFoilOpticalSurfaceReflect);
    //PTFEFoilOpticalSurface->SetMaterialPropertiesTable(reflfoilMPT);

    // now add logical skin surface
    new G4LogicalSkinSurface("PTFE_Surface", fReflectorFoil_log, ReflectorOpticalSurface);

    // and surface properties for the frill + bridle
    // Definition of MPT for Plastic frills
    std::vector<G4double> plastic_Energy = { 1.0*eV,1.2*eV,2.5*eV,3.0*eV,3.4*eV,6.5*eV,10.0*eV,12.6*eV };
    std::vector<G4double> plastic_reflect = {0.10, 0.10, 0.25, 0.30, 0.10, 0.05, 0.01, 0.01};

    G4OpticalSurface *PlasticOpticalSurface = new G4OpticalSurface("PlasticOpticalSurface");

    PlasticOpticalSurface->SetModel(unified);
    PlasticOpticalSurface->SetType(dielectric_dielectric);
    PlasticOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections

    G4MaterialPropertiesTable *Plastic_MT = new G4MaterialPropertiesTable();
    Plastic_MT->AddProperty("REFLECTIVITY", plastic_Energy, plastic_reflect);
    Plastic_MT->AddProperty("TRANSMITTANCE", TPBEnergy, TPBfoilOSTransmit);
    Plastic_MT->AddProperty("EFFICIENCY", TPBEnergy, TPBfoilOSEff);
    PlasticOpticalSurface->SetMaterialPropertiesTable(Plastic_MT);

    new G4LogicalSkinSurface("FrillWall_Surface", fFrillWall_log, PlasticOpticalSurface);
    new G4LogicalSkinSurface("FrillCaps_Surface", fFrillCaps_log, PlasticOpticalSurface);
    new G4LogicalSkinSurface("BridleWall_Surface", fBridleWall_log, PlasticOpticalSurface);
    new G4LogicalSkinSurface("BridleCaps_Surface", fBridleCaps_log, PlasticOpticalSurface);

    // let's also give the pmt glass some reflection properties
    G4OpticalSurface *CoatedPMTGlassOpticalSurface = new G4OpticalSurface("CoatedPMTGlassOpticalSurface");
    G4OpticalSurface *UncoatedPMTGlassOpticalSurface = new G4OpticalSurface("UncoatedPMTGlassOpticalSurface");

    // define uncoated pmts --> ground
    UncoatedPMTGlassOpticalSurface->SetModel(unified);
    UncoatedPMTGlassOpticalSurface->SetType(dielectric_dielectric);
    UncoatedPMTGlassOpticalSurface->SetFinish(groundfrontpainted); // 100% Lambertian (diffuse) reflections

    // define coated pmts --> groundfrontpainted (all diffuse since modelling TPB + ground glass reflections at once)
    CoatedPMTGlassOpticalSurface->SetModel(unified);
    CoatedPMTGlassOpticalSurface->SetType(dielectric_dielectric);
    CoatedPMTGlassOpticalSurface->SetFinish(groundfrontpainted);

    // define reflectivity for PMT glass (to be used on both coated and uncoated pmts)
    G4MaterialPropertiesTable *PMTGlassMPT = new G4MaterialPropertiesTable();
    std::vector<G4double> PMTGlassEnergy = {0.602*eV, 0.689*eV, 1.03*eV,  1.926*eV, 2.138*eV, 2.25*eV,  2.38*eV,
                                            2.48*eV,  2.583*eV, 2.845*eV, 2.857*eV, 2.95*eV,  3.124*eV, 3.457*eV,
                                            3.643*eV, 3.812*eV, 4.086*eV, 4.511*eV, 4.953*eV, 5.474*eV, 6.262*eV,
                                            7.000*eV, 8.300*eV, 10.00*eV, 12.60*eV };

    G4double uvTransmittance = 0.0;
    G4double vsTransmittance = 0.90;

    std::vector<G4double> PMTGlassTransmittance = { vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance,
                                                    vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance,
                                                    vsTransmittance, vsTransmittance, vsTransmittance, vsTransmittance, uvTransmittance, uvTransmittance, uvTransmittance,
                                                    uvTransmittance, uvTransmittance, uvTransmittance, uvTransmittance};

    G4double uvReflection = 0.06;
    G4double vsReflection = 0.10;

    std::vector<G4double> PMTGlassReflection = { vsReflection, vsReflection, vsReflection, vsReflection, vsReflection, vsReflection, vsReflection,
                                                    vsReflection, vsReflection, vsReflection, vsReflection, vsReflection, vsReflection, vsReflection,
                                                    vsReflection, vsReflection, vsReflection, vsReflection, uvReflection, uvReflection, uvReflection,
                                                    uvReflection, uvReflection, uvReflection, uvReflection};

    PMTGlassMPT->AddProperty("REFLECTIVITY", PMTGlassEnergy, PMTGlassReflection);
    PMTGlassMPT->AddProperty("TRANSMITTANCE", PMTGlassEnergy, PMTGlassTransmittance);

    CoatedPMTGlassOpticalSurface->SetMaterialPropertiesTable(PMTGlassMPT);
    UncoatedPMTGlassOpticalSurface->SetMaterialPropertiesTable(PMTGlassMPT);

    new G4LogicalSkinSurface("CoatedPMTGlassWall_Surface", fPMTCoatedWall_log, CoatedPMTGlassOpticalSurface);
    new G4LogicalSkinSurface("CoatedPMTGlassCaps_Surface", fPMTCoatedCaps_log, CoatedPMTGlassOpticalSurface);
    new G4LogicalSkinSurface("UncoatedPMTGlassWall_Surface", fPMTUncoatedWall_log, UncoatedPMTGlassOpticalSurface);
    new G4LogicalSkinSurface("UncoatedPMTGlassCaps_Surface", fPMTUncoatedCaps_log, UncoatedPMTGlassOpticalSurface);

    // finally, add surface properties for the shiny guys
    // we have the shiny reflectors at (0, 0) on top + bottom and shiny guy on top off center
    new G4LogicalSkinSurface("ShinyC406R0_Surface", fShinyC406R0_log, ReflectorOpticalSurface);
    new G4LogicalSkinSurface("ShinyTop_Surface", fShinyTop_log, ReflectorOpticalSurface);
    new G4LogicalSkinSurface("ShinyBottom_Surface", fShinyBottom_log, ReflectorOpticalSurface);
}

