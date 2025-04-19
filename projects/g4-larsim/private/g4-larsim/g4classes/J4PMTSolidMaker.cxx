// $Id: J4PMTSolidMaker.cc,v 1.2 2007/03/13 17:29:14 hoshina Exp $
/**
* @file J4PMTSolidMaker.cc
* @brief namespace that provides geometory functions to generate G4Solids
* @date 2021/04/11
* @note many numbers are hard coded, if you want to change it on the fly move these parameters to J4PartsParameterList.
* (Update Record)
*	2021/04/11  K.Hoshina	collected functions from each detector classes
*/

#include <set>

#include "g4-larsim/g4classes/J4PMTSolidMaker.h"

#include <G4Box.hh>
#include <G4Cons.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Polycone.hh>
#include <G4MultiUnion.hh>
#include <G4UserLimits.hh>
#include <G4SystemOfUnits.hh>
#include <G4OpticalSurface.hh>
#include <G4SubtractionSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4SolidStore.hh>

// let's define some things relevant for getting the geometry of our pmts
std::vector<std::string> pmt_ids = {"C101R0",
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

G4double ComputeDiskHeight(G4double cylinder_radius, G4double disk_radius, G4double disk_width) {
    G4double r0 = cylinder_radius - disk_width;
    G4double r1 = cylinder_radius;
    G4double r = disk_radius;
    G4double h = r1 - sqrt(r0*r0 - r*r);
    return h;
}

G4double WallSolidsOverlap(G4ThreeVector pmt_pos_0, G4ThreeVector pmt_pos_1, G4double inset_distance, G4double object_radius) {
    G4double z_dist = std::abs(pmt_pos_0.z() - pmt_pos_1.z());
    G4double x0 = pmt_pos_0.x(); G4double y0 = pmt_pos_0.y();
    G4double x1 = pmt_pos_1.x(); G4double y1 = pmt_pos_1.y();
    G4double dx = x1 - x0; G4double dy = y1 - y0;
    G4double xy_dist = sqrt(dx*dx + dy*dy);
    if(z_dist > 1.0 * cm) { // Not the same z-position
        if(xy_dist < 1.0 * cm) { // The same xy-position
            double distance = (pmt_pos_0 - pmt_pos_1).mag();
            return distance < 2*object_radius;
        } else {
            return z_dist < 2*object_radius; // They are not the same z-position but close enough
        }
    }

    if(xy_dist > 3.0 * object_radius) { // Too far apart in xy
        return false;
    }

    pmt_pos_0.setZ(0);
    pmt_pos_1.setZ(0);
    G4ThreeVector d0 = pmt_pos_0.unit();
    G4ThreeVector d1 = pmt_pos_1.unit();
    G4double R = (pmt_pos_0.mag() + pmt_pos_1.mag())/2.0 - inset_distance;
    G4double r = object_radius;
    double cos_theta = d0.dot(d1);
    G4double D = sqrt(2)*R*sqrt(1.0 - cos_theta);
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    G4double d = D - D*r/(R*sin_theta);
    return d <= 0;
}

G4double CapSolidsOverlap(G4ThreeVector pmt_pos_0, G4ThreeVector pmt_pos_1, G4double object_radius) {
    bool same_side = std::abs(pmt_pos_0.z() - pmt_pos_1.z()) < 1.0 * cm;
    if(not same_side)
        return false;
    double distance = (pmt_pos_0 - pmt_pos_1).mag();
    return distance < 2*object_radius;
}

//////////////////////////////////////////////////////////////////////////////////////
std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTPolyconeSolid(double scale_factor) {
    G4ThreeVector centerOfPolycone = G4ThreeVector(0, 0, 59.5 * scale_factor * mm);
    std::vector<G4double> tempZ, tempInner, tempOuter;

    G4int segment = 100;

    G4double sr0 = 57.0894 * scale_factor; // radius of spindle torus sphere
    G4double centerOfsr0 = 43.9106 * scale_factor; // distance from center of torus sphere to z-axis
    G4double rplanet = sr0 * 2. / segment; // z length of each segment

    G4double rInner[segment + 1], rOuter[segment + 1], zPlane[segment + 1];

    for (G4int j = 0; j <= segment; ++j) {
        tempZ.push_back((sr0 - j * rplanet) * mm);
        tempInner.push_back(0.);
        tempOuter.push_back((centerOfsr0 + sqrt(sr0 * sr0 - (sr0 - j * rplanet) * (sr0 - j * rplanet))) * mm);
    }

    for (G4int i = 0; i <= segment; i++) {
        rInner[i] = tempInner[i];
        rOuter[i] = tempOuter[i];
        zPlane[i] = tempZ[i];
    }

    G4Polycone* polycone1 = new G4Polycone("polycone1", 0, 2 * M_PI, segment + 1, zPlane, rInner, rOuter);
    temp_solids.emplace_back(polycone1);
    G4VSolid * solid = new G4DisplacedSolid("solid", polycone1, 0, centerOfPolycone);
    temp_solids.emplace_back(solid);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTSphereSolid(double scale_factor) {
    // G4ThreeVector centerOfSphere = G4ThreeVector(0, 0, 0);
    G4double rSmin = 0.0;
    G4double rSmax = pmt_radius * scale_factor;
    G4double dSphi = 2 * M_PI;
    G4double dStheta = 37.6392 * degree;

    G4VSolid * sphere = new G4Sphere("sphere", rSmin, rSmax, 0, dSphi, 0, dStheta);
    temp_solids.emplace_back(sphere);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTConeSolid(double scale_factor) {
    G4ThreeVector centerOfCons = G4ThreeVector(0, 0, (70.4284 + 33.8363 / 2 - 89.) * scale_factor * mm);
    G4double rmin = 0;
    G4double rmax = (84.5 / 2) * scale_factor * mm;
    G4double rmin2 = 0;
    G4double rmax2 = (160. / 2) * scale_factor * mm;
    G4double dz = (33.8363 / 2) * scale_factor * mm;
    G4double sphi = 0;
    G4double dphi = 2 * M_PI;

    G4VSolid * cons = new G4Cons("cons", rmin, rmax, rmin2, rmax2, dz, sphi, dphi);
    temp_solids.emplace_back(cons);
    G4VSolid * solid = new G4DisplacedSolid("solid", cons, 0, centerOfCons);
    temp_solids.emplace_back(solid);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTTubsSolid(double scale_factor) {
    G4ThreeVector centerOfTubs = G4ThreeVector(0, 0, -(89. - 70.4284 / 2) * scale_factor * mm);
    G4double rmin = 0;
    G4double rmax = (84.5 / 2) * scale_factor * mm;
    G4double dz = (70.4284 / 2) * scale_factor * mm;
    G4double sphi = 0;
    G4double dphi = 2 * M_PI;

    G4VSolid * tubs = new G4Tubs("tubs", rmin, rmax, dz, sphi, dphi);
    temp_solids.emplace_back(tubs);
    G4VSolid * solid = new G4DisplacedSolid("solid", tubs, 0, centerOfTubs);
    temp_solids.emplace_back(solid);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTSolid(double scale_factor) {
    std::shared_ptr<G4VSolid> polycone = CreatePMTPolyconeSolid(scale_factor);
    std::shared_ptr<G4VSolid> sphere = CreatePMTSphereSolid(scale_factor);
    std::shared_ptr<G4VSolid> cons = CreatePMTConeSolid(scale_factor);
    std::shared_ptr<G4VSolid> tubs = CreatePMTTubsSolid(scale_factor);

    G4Transform3D transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, 0));

    G4MultiUnion * solid = new G4MultiUnion("solid");
    solid->AddNode(*sphere, transform);
    solid->AddNode(*polycone, transform);
    solid->AddNode(*tubs, transform);
    solid->AddNode(*cons, transform);
    solid->Voxelize();
    temp_solids.emplace_back(solid);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTGlassFilledSolid() {
    if(pmt_glass_filled_solid)
        throw std::runtime_error("PMT glass solid already exists");
    pmt_glass_filled_solid = CreatePMTSolid(1.0);
    return pmt_glass_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTTPBFilledSolid() {
    if(pmt_tpb_filled_solid)
        throw std::runtime_error("PMT TPB solid already exists");
    double scale_factor = (pmt_radius + pmt_tpb_thickness) / pmt_radius;
    pmt_tpb_filled_solid = CreatePMTSolid(scale_factor);
    return pmt_tpb_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTVacuumFilledSolid() {
    if(pmt_vacuum_filled_solid)
        throw std::runtime_error("PMT vacuum solid already exists");
    double scale_factor = (pmt_radius - pmt_glass_thickness) / pmt_radius;
    pmt_vacuum_filled_solid = CreatePMTSolid(scale_factor);
    return pmt_vacuum_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateCylinderCurvatureSolid() {
    if(cylinder_curvature_solid)
        throw std::runtime_error("Cylinder curvature solid already exists");
    G4Tubs * tub = new G4Tubs("DetectorCylinderCurvature", 0*cm, cylinder_radius, pmt_radius * 2, 0*deg, 360*deg);
    temp_solids.emplace_back(tub);

    cylinder_curvature_rotation = new G4RotationMatrix();
    cylinder_curvature_rotation->rotateY(M_PI / 2.0); // Rotate 90degrees around the x axis
    cylinder_curvature_position = G4ThreeVector(0, 0, pmt_radius - pmt_protrusion_distance + cylinder_radius);

    cylinder_curvature_solid = std::shared_ptr<G4VSolid>(tub);
    return cylinder_curvature_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateCylinderFaceSolid() {
    if(cylinder_face_solid)
        throw std::runtime_error("Cylinder face solid already exists");
    G4double radius = frill_radius * 2.0;
    G4double half_height = pmt_protrusion_distance / 2.0;
    G4Tubs * tub = new G4Tubs("DetectorCylinderFace", 0*cm, frill_radius * 2.0, half_height, 0*deg, 360*deg);
    temp_solids.emplace_back(tub);

    cylinder_face_rotation = new G4RotationMatrix();
    cylinder_face_position = G4ThreeVector(0, 0, half_height + pmt_radius - pmt_protrusion_distance);

    cylinder_face_solid = std::shared_ptr<G4VSolid>(tub);
    return cylinder_face_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateWallSolid(std::shared_ptr<G4VSolid> pmt_solid) {
    G4VSolid * solid = new G4IntersectionSolid("WallSolid" + pmt_solid->GetName(), pmt_solid.get(), cylinder_curvature_solid.get(), cylinder_curvature_rotation, cylinder_curvature_position);
    temp_solids.emplace_back(solid);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateCapSolid(std::shared_ptr<G4VSolid> pmt_solid) {
    G4VSolid * solid = new G4IntersectionSolid("CapSolid" + pmt_solid->GetName(), pmt_solid.get(), cylinder_face_solid.get(), cylinder_face_rotation, cylinder_face_position);
    temp_solids.emplace_back(solid);
    return temp_solids.back();
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTWallGlassFilledSolid() {
    pmt_wall_glass_filled_solid = CreateWallSolid(pmt_glass_filled_solid);
    return pmt_wall_glass_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTWallTPBFilledSolid() {
    pmt_wall_tpb_filled_solid = CreateWallSolid(pmt_tpb_filled_solid);
    return pmt_wall_tpb_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTWallVacuumFilledSolid() {
    pmt_wall_vacuum_filled_solid = CreateWallSolid(pmt_vacuum_filled_solid);
    return pmt_wall_vacuum_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTCapGlassFilledSolid() {
    pmt_cap_glass_filled_solid = CreateCapSolid(pmt_glass_filled_solid);
    return pmt_cap_glass_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTCapTPBFilledSolid() {
    pmt_cap_tpb_filled_solid = CreateCapSolid(pmt_tpb_filled_solid);
    return pmt_cap_tpb_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTCapVacuumFilledSolid() {
    pmt_cap_vacuum_filled_solid = CreateCapSolid(pmt_vacuum_filled_solid);
    return pmt_cap_vacuum_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleFilledSolid() {
    if(bridle_filled_solid)
        throw std::runtime_error("Frill center solid already exists");
    G4double max_disk_height = ComputeDiskHeight(cylinder_radius, bridle_radius, bridle_width) + bridle_width;

    G4VSolid * fBridlePunchOut = new G4Tubs("BridlePunchOut", 0*cm, bridle_radius, max_disk_height*4, 0*deg, 360*deg);
    temp_solids.emplace_back(fBridlePunchOut);
    bridle_filled_solid = temp_solids.back();
    return bridle_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleCylinderSolid() {
    G4VSolid * fBridleSurface = new G4Tubs("BridleSurface", cylinder_radius - bridle_width, cylinder_radius, bridle_radius * 2, 0*deg, 360*deg);
    temp_solids.emplace_back(fBridleSurface);
    bridle_cylinder_solid = temp_solids.back();
    return bridle_cylinder_solid;
}


std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillFilledSolid() {
    G4double max_disk_height = ComputeDiskHeight(cylinder_radius, frill_radius, frill_width) + frill_width;

    G4VSolid * fFrillPunchOut = new G4Tubs("FrillPunchOut", 0, frill_radius, max_disk_height*4, 0*deg, 360*deg);
    temp_solids.emplace_back(fFrillPunchOut);
    frill_filled_solid = temp_solids.back();
    return frill_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillCylinderSolid() {
    G4VSolid * fFrillSurface = new G4Tubs("FrillSurface", cylinder_radius - frill_width, cylinder_radius, frill_radius * 2, 0*deg, 360*deg);
    temp_solids.emplace_back(fFrillSurface);
    frill_cylinder_solid = temp_solids.back();
    return frill_cylinder_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleWallFilledSolid() {
    // Creates a thin disk of material that has a circular shape in its X/Y projection and curvature in the Y/Z plane
    // Curvature bends into the positive Z direction
    // Origin of solid is in line with its X/Y center and at the lowest tangent point of its exterior curvature


    // bridle_radius --> the radius of the disk --> r2
    // bridle_width --> the thickness of the disk
    // cylinder_radius --> the radius of the drum that we are punching the disk out of --> r1


    // Solve for the height of the disk
    //G4double max_disk_height = ComputeDiskHeight(cylinder_radius, bridle_radius, bridle_width);

    //G4VSolid * fBridlePunchOut = new G4Tubs("BridlePunchOut", 0*cm, bridle_radius, max_disk_height*4, 0*deg, 360*deg);
    //temp_solids.emplace_back(fBridlePunchOut);
    G4VSolid * fBridlePunchOut = bridle_filled_solid.get();

    // Create the surface that we are punching the disk out of
    G4VSolid * fBridleSurface = bridle_cylinder_solid.get();
    //G4VSolid * fBridleSurface = new G4Tubs("BridleSurface", cylinder_radius - bridle_width, cylinder_radius, bridle_radius * 2, 0*deg, 360*deg);
    //temp_solids.emplace_back(fBridleSurface);

    // Need to have the surface curvature in the Y/Z plane
    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(M_PI / 2.0); // Rotate 90degrees around the Y axis

    // We want the outer edge of the surface cylinder to touch the origin of the punchout
    G4double z_center = cylinder_radius;
    G4ThreeVector centerOfTub(0, 0, z_center);
    G4VSolid * disk = new G4IntersectionSolid("BridlePunchOut", fBridlePunchOut, fBridleSurface, rotationMatrix, centerOfTub);
    temp_solids.emplace_back(disk);

    // Now we need to offset relative to the PMT center
    G4double displacement = pmt_radius - pmt_protrusion_distance;
    G4ThreeVector disk_offset(0, 0, displacement);

    G4DisplacedSolid * solid = new G4DisplacedSolid("Bridle", disk, 0, disk_offset);
    temp_solids.emplace_back(solid);
    bridle_wall_filled_solid = temp_solids.back();
    delete rotationMatrix;
    return bridle_wall_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleCapFilledSolid() {
    G4VSolid * fBridleDisk = new G4Tubs("BridleDisk", 0*cm, bridle_radius, bridle_width/2.0, 0*deg, 360*deg);
    temp_solids.emplace_back(fBridleDisk);

    G4double displacement = pmt_radius - pmt_protrusion_distance + bridle_width / 2.0;
    G4ThreeVector disk_offset(0, 0, displacement);

    G4DisplacedSolid * solid = new G4DisplacedSolid("Bridle", fBridleDisk, 0, disk_offset);
    temp_solids.emplace_back(solid);
    bridle_cap_filled_solid = temp_solids.back();
    return bridle_cap_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillWallFilledSolid() {
    //G4VSolid * fFrillPunchOut = new G4Tubs("FrillPunchOut", 0, frill_radius, max_disk_height*4, 0*deg, 360*deg);
    //temp_solids.emplace_back(fFrillPunchOut);
    G4VSolid * fFrillPunchOut = frill_filled_solid.get();
    //G4VSolid * fFrillPunchOut = frill_center_solid.get(); // Use the same solid as the center of the frill
    //G4VSolid * fFrillSurface = new G4Tubs("FrillSurface", cylinder_radius - frill_width, cylinder_radius, frill_radius * 2, 0*deg, 360*deg);
    //temp_solids.emplace_back(fFrillSurface);
    G4VSolid * fFrillSurface = frill_cylinder_solid.get();

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(M_PI / 2.0); // Rotate 90degrees around the Y axis
                                         //
    G4double z_center = cylinder_radius;
    G4ThreeVector centerOfTub(0, 0, z_center);
    G4VSolid * disk = new G4IntersectionSolid("FrillPunchOut", fFrillPunchOut, fFrillSurface, rotationMatrix, centerOfTub);
    temp_solids.emplace_back(disk);

    // Now we need to offset relative to the PMT center
    G4double displacement = pmt_radius - pmt_protrusion_distance;
    G4ThreeVector disk_offset(0, 0, displacement);

    G4DisplacedSolid * solid = new G4DisplacedSolid("Frill", disk, 0, disk_offset);
    temp_solids.emplace_back(solid);
    frill_wall_filled_solid = temp_solids.back();
    delete rotationMatrix;
    return frill_wall_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillCapFilledSolid() {
    G4VSolid * fFrillDisk = new G4Tubs("FrillDisk", 0*cm, frill_radius, frill_width/2.0, 0*deg, 360*deg);
    temp_solids.emplace_back(fFrillDisk);

    G4double displacement = pmt_radius - pmt_protrusion_distance + frill_width / 2.0;
    G4ThreeVector disk_offset(0, 0, displacement);

    G4DisplacedSolid * solid = new G4DisplacedSolid("Frill", fFrillDisk, 0, disk_offset);
    temp_solids.emplace_back(solid);
    frill_cap_filled_solid = temp_solids.back();
    return frill_cap_filled_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTGlassSolid() {
    if(pmt_glass_solid)
        throw std::runtime_error("PMT glass solid already exists");
    temp_solids.emplace_back(
        new G4SubtractionSolid("PMTGlassHollow",
            pmt_glass_filled_solid.get(),
            pmt_vacuum_filled_solid.get()
        )
    );
    pmt_glass_solid = temp_solids.back();
    return pmt_glass_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTTPBSolid() {
    if(pmt_tpb_solid)
        throw std::runtime_error("PMT TPB solid already exists");
    temp_solids.emplace_back(
        new G4SubtractionSolid("PMTTPBHollow",
            pmt_tpb_filled_solid.get(),
            pmt_glass_filled_solid.get()
        )
    );
    pmt_tpb_solid = temp_solids.back();
    return pmt_tpb_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTWallGlassSolid() {
    if(pmt_wall_glass_solid)
        throw std::runtime_error("PMT wall glass solid already exists");
    pmt_wall_glass_solid = CreateWallSolid(pmt_glass_solid);
    return pmt_wall_glass_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTCapGlassSolid() {
    if(pmt_cap_glass_solid)
        throw std::runtime_error("PMT cap glass solid already exists");
    pmt_cap_glass_solid = CreateCapSolid(pmt_glass_solid);
    return pmt_cap_glass_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTWallTPBSolid() {
    if(pmt_wall_tpb_solid)
        throw std::runtime_error("PMT wall TPB solid already exists");
    pmt_wall_tpb_solid = CreateWallSolid(pmt_tpb_solid);
    return pmt_wall_tpb_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreatePMTCapTPBSolid() {
    if(pmt_cap_tpb_solid)
        throw std::runtime_error("PMT cap TPB solid already exists");
    pmt_cap_tpb_solid = CreateCapSolid(pmt_tpb_solid);
    return pmt_cap_tpb_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleWallCoatedSolid() {
    if(bridle_wall_coated_solid)
        throw std::runtime_error("Bridle wall coated solid already exists");
    temp_solids.emplace_back(
        new G4SubtractionSolid("BridleWallCoated",
            bridle_wall_filled_solid.get(),
            pmt_tpb_filled_solid.get()
        )
    );
    bridle_wall_coated_solid = temp_solids.back();
    return bridle_wall_coated_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleWallUncoatedSolid() {
    if(bridle_wall_uncoated_solid)
        throw std::runtime_error("Bridle wall uncoated solid already exists");
    temp_solids.emplace_back(
        new G4SubtractionSolid("BridleWallUncoated",
            bridle_wall_filled_solid.get(),
            pmt_glass_filled_solid.get()
        )
    );
    bridle_wall_uncoated_solid = temp_solids.back();
    return bridle_wall_uncoated_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleCapCoatedSolid() {
    if(bridle_cap_coated_solid)
        throw std::runtime_error("Bridle cap coated solid already exists");
    temp_solids.emplace_back(
        new G4SubtractionSolid("BridleCapCoated",
            bridle_cap_filled_solid.get(),
            pmt_tpb_filled_solid.get()
        )
    );
    bridle_cap_coated_solid = temp_solids.back();
    return bridle_cap_coated_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateBridleCapUncoatedSolid() {
    if(bridle_cap_uncoated_solid)
        throw std::runtime_error("Bridle cap uncoated solid already exists");
    temp_solids.emplace_back(
        new G4SubtractionSolid("BridleCapUncoated",
            bridle_cap_filled_solid.get(),
            pmt_glass_filled_solid.get()
        )
    );
    bridle_cap_uncoated_solid = temp_solids.back();
    return bridle_cap_uncoated_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillCenterSolid() {
    frill_center_solid = bridle_filled_solid;
    return frill_center_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillWallSolid() {
    if(frill_wall_solid)
        throw std::runtime_error("Frill wall solid already exists");
    //G4VSolid * fBridlePunchOut = new G4Tubs("BridlePunchOut", 0*cm, bridle_radius, frill_width*2, 0*deg, 360*deg);
    //temp_solids.emplace_back(fBridlePunchOut);
    G4double displacement = pmt_radius - pmt_protrusion_distance + frill_width/2.0;
    G4ThreeVector disk_offset(0, 0, displacement);
    temp_solids.emplace_back(
        new G4SubtractionSolid("FrillWall",
            frill_wall_filled_solid.get(),
            frill_center_solid.get(),
            0,
            disk_offset
        )
    );
    frill_wall_solid = temp_solids.back();
    return frill_wall_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateFrillCapSolid() {
    if(frill_cap_solid)
        throw std::runtime_error("Frill cap solid already exists");
    //G4VSolid * fBridlePunchOut = new G4Tubs("BridlePunchOut", 0*cm, bridle_radius, frill_width*2, 0*deg, 360*deg);
    //temp_solids.emplace_back(fBridlePunchOut);
    G4double displacement = pmt_radius - pmt_protrusion_distance + frill_width;
    G4ThreeVector disk_offset(0, 0, displacement);
    temp_solids.emplace_back(
        new G4SubtractionSolid("FrillCap",
            frill_cap_filled_solid.get(),
            frill_center_solid.get(),
            0,
            disk_offset
        )
    );
    frill_cap_solid = temp_solids.back();
    return frill_cap_solid;
}

void J4PMTSolidMaker::ParsePMTID(std::string pmt_id, bool & pmt_on_cap, int & pmt_row, int & pmt_col, int & pmt_ring, int & pmt_number, bool & coated) {
    size_t r_pos = pmt_id.find("R");
    size_t n_row_chars = r_pos - 1;
    pmt_on_cap = n_row_chars > 2;
    pmt_row = 0;
    pmt_col = 0;
    pmt_ring = 0;
    pmt_number = 0;
    coated = true;
    if(pmt_on_cap) {
        pmt_col = std::atoi(pmt_id.substr(2, n_row_chars - 1).c_str());
        pmt_ring = std::atoi(pmt_id.substr(1, 1).c_str());
        pmt_row = std::atoi(pmt_id.substr(r_pos + 1, std::string::npos).c_str());
        for(std::pair<size_t const, size_t> const & p : ring_pmt_pos_count)
            if(p.first < pmt_ring and p.first > 0)
                pmt_number += p.second;
        pmt_number += pmt_col;
        std::tuple<int, int, int> cap_coated_key = {pmt_ring, pmt_col, pmt_row};
        if(cap_uncoated_pmts.count(cap_coated_key) > 0)
            coated = false;
    } else {
        pmt_col = std::atoi(pmt_id.substr(1, r_pos - 1).c_str());
        pmt_row = std::atoi(pmt_id.substr(r_pos + 1, std::string::npos).c_str());
        pmt_number = pmt_col;
        if( (pmt_row ==  1 and pmt_col % 4 == 1) or
            (pmt_row ==  2 and pmt_col % 4 == 3) or
            (pmt_row ==  4 and pmt_col % 4 == 0) or
            (pmt_row ==  5 and pmt_col % 4 == 2)) {
            coated = false;
        }
    }
}

G4RotationMatrix* J4PMTSolidMaker::GetPMTRotationMatrix(std::string pmt_id) {
    bool pmt_on_cap, coated;
    int pmt_row, pmt_col, pmt_ring, pmt_number;
    ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    if(pmt_row == 0) {
        // top face pmts
        // rotate so they are pointing downwards
        G4double rotationAngle = M_PI;
        rotationMatrix->rotateX(rotationAngle); // Rotate around the x axis
    }
    // note -- don't actually need to do a rotation for bottom face pmts -- skipping for now
    if(pmt_row > 0 and pmt_row < 6) {
        G4ThreeVector position = GetPMTPosition(pmt_id);
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
    return rotationMatrix;
}

G4ThreeVector J4PMTSolidMaker::GetPMTPosition(std::string pmt_id) {
    bool pmt_on_cap, coated;
    int pmt_row, pmt_col, pmt_ring, pmt_number;
    ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

    double pmt_radius_cm = (cylinder_radius + (pmt_radius - pmt_protrusion_distance)) / cm;
    double pmt_height_cm = (cylinder_half_height + (pmt_radius - pmt_protrusion_distance)) / cm;

    std::vector<double> position;
    if(pmt_on_cap) {
        if(pmt_row == 0) {
            position = get_pmt_cap_position(pmt_row, pmt_ring, pmt_col, 1, 0.0, pmt_height_cm);
        } else if(pmt_row == 6) {
            position = get_pmt_cap_position(pmt_row, pmt_ring, pmt_col, 1, 0.0, -pmt_height_cm);
        } else {
            throw std::runtime_error("Bad cap row");
        }
    } else {
        position = get_pmt_wall_position(pmt_row, pmt_col, 1, 0.0, pmt_radius_cm);
    }
    G4ThreeVector pmt_pos;
    pmt_pos.setX(position[0]*cm);
    pmt_pos.setY(position[1]*cm);
    pmt_pos.setZ(position[2]*cm);
    return pmt_pos;
}

void J4PMTSolidMaker::CreateBaseSolids() {
    CreateCylinderFaceSolid();
    CreateCylinderCurvatureSolid();
    CreateBridleFilledSolid();
    CreateBridleCylinderSolid();
    CreateFrillFilledSolid();
    CreateFrillCylinderSolid();
    std::shared_ptr<G4VSolid> pmt_vacuum_filled_solid = CreatePMTVacuumFilledSolid();
    std::shared_ptr<G4VSolid> pmt_glass_filled_solid = CreatePMTGlassFilledSolid();
    std::shared_ptr<G4VSolid> pmt_tpb_filled_solid = CreatePMTTPBFilledSolid();
    std::shared_ptr<G4VSolid> pmt_wall_glass_filled_solid = CreatePMTWallGlassFilledSolid();
    std::shared_ptr<G4VSolid> pmt_wall_tpb_filled_solid = CreatePMTWallTPBFilledSolid();
    std::shared_ptr<G4VSolid> pmt_wall_vacuum_filled_solid = CreatePMTWallVacuumFilledSolid();
    std::shared_ptr<G4VSolid> pmt_cap_glass_filled_solid = CreatePMTCapGlassFilledSolid();
    std::shared_ptr<G4VSolid> pmt_cap_tpb_filled_solid = CreatePMTCapTPBFilledSolid();
    std::shared_ptr<G4VSolid> pmt_cap_vacuum_filled_solid = CreatePMTCapVacuumFilledSolid();
    std::shared_ptr<G4VSolid> frill_center_solid = CreateFrillCenterSolid();
    std::shared_ptr<G4VSolid> bridle_wall_filled_solid = CreateBridleWallFilledSolid();
    std::shared_ptr<G4VSolid> bridle_cap_filled_solid = CreateBridleCapFilledSolid();
    std::shared_ptr<G4VSolid> frill_wall_filled_solid = CreateFrillWallFilledSolid();
    std::shared_ptr<G4VSolid> frill_cap_filled_solid = CreateFrillCapFilledSolid();

    std::shared_ptr<G4VSolid> pmt_glass_solid = CreatePMTGlassSolid();
    std::shared_ptr<G4VSolid> pmt_tpb_solid = CreatePMTTPBSolid();
    std::shared_ptr<G4VSolid> pmt_wall_glass_solid = CreatePMTWallGlassSolid();
    std::shared_ptr<G4VSolid> pmt_cap_glass_solid = CreatePMTCapGlassSolid();
    std::shared_ptr<G4VSolid> pmt_wall_tpb_solid = CreatePMTWallTPBSolid();
    std::shared_ptr<G4VSolid> pmt_cap_tpb_solid = CreatePMTCapTPBSolid();

    std::shared_ptr<G4VSolid> bridle_wall_coated_solid = CreateBridleWallCoatedSolid();
    std::shared_ptr<G4VSolid> bridle_wall_uncoated_solid = CreateBridleWallUncoatedSolid();
    std::shared_ptr<G4VSolid> bridle_cap_coated_solid = CreateBridleCapCoatedSolid();
    std::shared_ptr<G4VSolid> bridle_cap_uncoated_solid = CreateBridleCapUncoatedSolid();
    std::shared_ptr<G4VSolid> frill_wall_solid = CreateFrillWallSolid();
    std::shared_ptr<G4VSolid> frill_cap_solid = CreateFrillCapSolid();

    for(std::string const & pmt_id : pmt_ids) {
        bool pmt_on_cap, coated;
        int pmt_row, pmt_col, pmt_ring, pmt_number;
        ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

        std::shared_ptr<G4VSolid> vacuum_solid = (pmt_on_cap) ? pmt_cap_vacuum_filled_solid : pmt_wall_vacuum_filled_solid;
        std::shared_ptr<G4VSolid> glass_solid = (pmt_on_cap) ? pmt_cap_glass_filled_solid : pmt_wall_glass_filled_solid;
        std::shared_ptr<G4VSolid> tpb_solid = (coated) ? ((pmt_on_cap) ? pmt_cap_tpb_filled_solid : pmt_wall_tpb_filled_solid) : nullptr;
        std::shared_ptr<G4VSolid> bridle_solid = (pmt_on_cap) ? bridle_cap_filled_solid : bridle_wall_filled_solid;
        std::shared_ptr<G4VSolid> frill_solid = (pmt_on_cap) ? frill_cap_filled_solid : frill_wall_filled_solid;

        std::shared_ptr<G4RotationMatrix> rotationMatrix(GetPMTRotationMatrix(pmt_id));
        G4ThreeVector position = GetPMTPosition(pmt_id);

        tpb_solids.insert({pmt_id, tpb_solid});
        tpb_positions.insert({pmt_id, position});
        tpb_rotations.insert({pmt_id, rotationMatrix});
        pmt_solids.insert({pmt_id, glass_solid});
        pmt_positions.insert({pmt_id, position});
        pmt_rotations.insert({pmt_id, rotationMatrix});
        vacuum_solids.insert({pmt_id, vacuum_solid});
        vacuum_positions.insert({pmt_id, position});
        vacuum_rotations.insert({pmt_id, rotationMatrix});

        std::shared_ptr<G4VSolid> bridle_punchout_solid = (coated) ? pmt_tpb_filled_solid : pmt_glass_filled_solid;
        std::shared_ptr<G4VSolid> frill_punchout_solid = frill_center_solid;

        base_solids.insert({pmt_id, {{bridle_solid, bridle_punchout_solid, frill_solid, frill_punchout_solid}, rotationMatrix, position}});
        // transformed_solids.insert({
        //     pmt_id, {
        //         new G4DisplacedSolid("PMTVacuum", pmt_vacuum_filled_solid, rotationMatrix, position),
        //         new G4DisplacedSolid("PMTGlass", pmt_glass_filled_solid, rotationMatrix, position),
        //         new G4DisplacedSolid("PMTTPB", tpb_solid, rotationMatrix, position),
        //         new G4DisplacedSolid("Bridle", bridle_solid, rotationMatrix, position),
        //         new G4DisplacedSolid("Frill", frill_solid, rotationMatrix, position)
        //     },
        //     rotationMatrix,
        //     position
        // });
    }

    // 1. Create a dummy mother volume that is large enough to contain your daughters.
    G4double motherHalfSize = 10.0*m;  // Adjust size as needed.
    G4Box * motherSolid = new G4Box("MotherBox", motherHalfSize, motherHalfSize, motherHalfSize);
    G4Material* dummyMaterial = new G4Material("DummyMaterial", 1.0*g/cm3, 1);
    dummyMaterial->AddElement(G4NistManager::Instance()->FindOrBuildElement("C"), 1);
    G4LogicalVolume * motherLV = new G4LogicalVolume(motherSolid, dummyMaterial, "MotherLV");
    G4PVPlacement * motherPV = new G4PVPlacement(0, G4ThreeVector(), motherLV, "MotherPV", 0, false, 0, false);

    std::function<bool(G4VSolid *, G4VSolid *)> test_overlap = [&](G4VSolid * solidA, G4VSolid * solidB) {
        G4ThreeVector pos(0, 0, 0);
        G4RotationMatrix * rot = new G4RotationMatrix(); // Identity
        G4LogicalVolume * lvA = new G4LogicalVolume(solidA, dummyMaterial, "LV_A");
        G4LogicalVolume * lvB = new G4LogicalVolume(solidB, dummyMaterial, "LV_B");
        G4PVPlacement * pvA = new G4PVPlacement(rot, pos, lvA, "PV_A", motherLV, false, 0, false);
        G4PVPlacement * pvB = new G4PVPlacement(rot, pos, lvB, "PV_B", motherLV, false, 0, false);
        bool overlaps = motherPV->CheckOverlaps(int(1e6));

        G4PhysicalVolumeStore::GetInstance()->DeRegister(pvA);
        G4PhysicalVolumeStore::GetInstance()->DeRegister(pvB);
        G4LogicalVolumeStore::GetInstance()->DeRegister(lvA);
        G4LogicalVolumeStore::GetInstance()->DeRegister(lvB);
        G4SolidStore::GetInstance()->DeRegister(solidA);
        G4SolidStore::GetInstance()->DeRegister(solidB);

        delete rot;
        delete lvA;
        delete lvB;
        delete pvA;
        delete pvB;

        return overlaps;
    };

    unsigned int next_frill_group_id = 0;
    unsigned int next_bridle_group_id = 0;
    std::map<unsigned int, std::vector<std::string>> frill_overlap_groups;
    std::map<unsigned int, std::vector<std::string>> bridle_overlap_groups;
    std::map<std::string, unsigned int> frill_overlaps;
    std::map<std::string, unsigned int> bridle_overlaps;

    for(std::string const & pmt_id_0 : pmt_ids) {
        std::tuple<std::array<std::shared_ptr<G4VSolid>, 5>, std::shared_ptr<G4RotationMatrix>, G4ThreeVector> const & base_solid_0 = base_solids[pmt_id_0];
        G4VSolid * bridle_solid_0 = std::get<0>(base_solid_0)[0].get();
        G4VSolid * frill_solid_0 = std::get<0>(base_solid_0)[2].get();
        G4RotationMatrix * rotationMatrix_0 = std::get<1>(base_solid_0).get();
        G4ThreeVector position_0 = std::get<2>(base_solid_0);
        bool pmt_on_cap_0, coated_0;
        int pmt_row_0, pmt_col_0, pmt_ring_0, pmt_number_0;
        ParsePMTID(pmt_id_0, pmt_on_cap_0, pmt_row_0, pmt_col_0, pmt_ring_0, pmt_number_0, coated_0);

        G4DisplacedSolid * frill_0 = new G4DisplacedSolid("Frill 0", frill_solid_0, rotationMatrix_0, position_0);
        G4DisplacedSolid * bridle_0 = new G4DisplacedSolid("Bridle 0", bridle_solid_0, rotationMatrix_0, position_0);

        for(std::string const & pmt_id_1 : pmt_ids) {
            std::tuple<std::array<std::shared_ptr<G4VSolid>, 5>, std::shared_ptr<G4RotationMatrix>, G4ThreeVector> const & base_solid_1 = base_solids[pmt_id_1];
            G4VSolid * bridle_solid_1 = std::get<0>(base_solid_1)[0].get();
            G4VSolid * frill_solid_1 = std::get<0>(base_solid_1)[2].get();
            G4RotationMatrix * rotationMatrix_1 = std::get<1>(base_solid_1).get();
            G4ThreeVector position_1 = std::get<2>(base_solid_1);
            bool pmt_on_cap_1, coated_1;
            int pmt_row_1, pmt_col_1, pmt_ring_1, pmt_number_1;
            ParsePMTID(pmt_id_1, pmt_on_cap_1, pmt_row_1, pmt_col_1, pmt_ring_1, pmt_number_1, coated_1);

            if(pmt_id_0 == pmt_id_1)
                continue;
            if(pmt_on_cap_0 != pmt_on_cap_1)
                continue;

            G4DisplacedSolid * frill_1 = new G4DisplacedSolid("Frill 1", frill_solid_1, rotationMatrix_1, position_1);
            G4DisplacedSolid * bridle_1 = new G4DisplacedSolid("Bridle 1", bridle_solid_1, rotationMatrix_1, position_1);

            bool frill_overlap = false;
            bool bridle_overlap = false;

            if(pmt_on_cap_0) {
                frill_overlap =
                    CapSolidsOverlap(position_0, position_1, frill_radius) or
                    test_overlap(frill_0, frill_1);
                bridle_overlap = frill_overlap and (
                    CapSolidsOverlap(position_0, position_1, bridle_radius) or
                    test_overlap(bridle_0, bridle_1));
            } else {
                G4double frill_inset_distance = pmt_radius - pmt_protrusion_distance + frill_width;
                frill_overlap =
                    WallSolidsOverlap(position_0, position_1, frill_inset_distance, frill_radius) or
                    test_overlap(frill_0, frill_1);
                G4double bridle_inset_distance = pmt_radius - pmt_protrusion_distance + bridle_width;
                bridle_overlap = frill_overlap and (
                    WallSolidsOverlap(position_0, position_1, bridle_inset_distance, bridle_radius) or
                    test_overlap(bridle_0, bridle_1));
            }

            delete frill_1;
            delete bridle_1;

            std::vector<std::tuple<
                bool,
                unsigned int &,
                std::map<std::string, unsigned int> &,
                std::map<unsigned int, std::vector<std::string>> &
            >> overlaps = {
                {frill_overlap, next_frill_group_id, frill_overlaps, frill_overlap_groups},
                {bridle_overlap, next_bridle_group_id, bridle_overlaps, bridle_overlap_groups}
            };
            for(std::tuple<bool, unsigned int &, std::map<std::string, unsigned int> &, std::map<unsigned int, std::vector<std::string>> &> & overlap : overlaps) {
                bool object_overlap = std::get<0>(overlap);
                unsigned int & next_group_id = std::get<1>(overlap);
                std::map<std::string, unsigned int> & object_overlaps = std::get<2>(overlap);
                std::map<unsigned int, std::vector<std::string>> & overlap_groups = std::get<3>(overlap);
                if(object_overlap) {
                    bool pmt_0_counted = object_overlaps.count(pmt_id_0) > 0;
                    bool pmt_1_counted = object_overlaps.count(pmt_id_1) > 0;
                    if(pmt_0_counted != pmt_1_counted) {
                        unsigned int group_id = pmt_0_counted ? object_overlaps[pmt_id_0] : object_overlaps[pmt_id_1];
                        if(pmt_0_counted) {
                            object_overlaps[pmt_id_1] = group_id;
                            overlap_groups[group_id].push_back(pmt_id_1);
                        } else {
                            object_overlaps[pmt_id_0] = group_id;
                            overlap_groups[group_id].push_back(pmt_id_0);
                        }
                    } else if(!pmt_0_counted and !pmt_1_counted) {
                        unsigned int group_id = next_group_id;
                        next_group_id += 1;
                        object_overlaps[pmt_id_0] = group_id;
                        object_overlaps[pmt_id_1] = group_id;
                        overlap_groups[group_id] = {pmt_id_0, pmt_id_1};
                    } else {
                        unsigned int group_id_0 = object_overlaps[pmt_id_0];
                        unsigned int group_id_1 = object_overlaps[pmt_id_1];
                        std::vector<std::string> & group_0 = overlap_groups[group_id_0];
                        std::vector<std::string> & group_1 = overlap_groups[group_id_1];
                        if(group_id_0 != group_id_1) {
                            for(std::string const & pmt_id : group_1) {
                                group_0.push_back(pmt_id);
                                object_overlaps[pmt_id] = group_id_0;
                            }
                            overlap_groups.erase(group_id_1);
                        }
                    }
                } else {
                    if(object_overlaps.count(pmt_id_0) == 0) {
                        unsigned int group_id = next_group_id;
                        next_group_id += 1;
                        object_overlaps[pmt_id_0] = group_id;
                        overlap_groups[group_id] = {pmt_id_0};
                    }
                    if(object_overlaps.count(pmt_id_1) == 0) {
                        unsigned int group_id = next_group_id;
                        next_group_id += 1;
                        object_overlaps[pmt_id_1] = group_id;
                        overlap_groups[group_id] = {pmt_id_1};
                    }
                }
            }
        }
        delete frill_0;
        delete bridle_0;
    }

    // frill_overlap_groups
    for(std::pair<unsigned int const, std::vector<std::string>> const & overlap_group : frill_overlap_groups) {
        std::vector<std::string> const & group = overlap_group.second;
        std::string group_id = std::to_string(overlap_group.first);
        if(group.size() == 1) {
            std::string pmt_id = group[0];
            bool pmt_on_cap, coated;
            int pmt_row, pmt_col, pmt_ring, pmt_number;
            ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

            std::string name = "Frill_" + pmt_id;
            std::shared_ptr<G4VSolid> frill;
            if(pmt_on_cap)
                frill = frill_cap_solid;
            else
                frill = frill_wall_solid;
            //std::shared_ptr<G4VSolid> frill_group_punchout(new G4SubtractionSolid("Frill Group Punchout " + group_id,
            //        std::get<0>(base_solids[pmt_id])[4].get(),
            //        std::get<0>(base_solids[pmt_id])[3].get()));
            frill_solids[name] = frill;
            frill_rotations[name] = std::get<1>(base_solids[pmt_id]);
            frill_positions[name] = std::get<2>(base_solids[pmt_id]);
            continue;
        }
        G4MultiUnion * frill_group = new G4MultiUnion("Frill Group " + group_id);
        temp_solids.emplace_back(frill_group);
        G4MultiUnion * punchout_group = new G4MultiUnion("Punchout Group " + group_id);
        temp_solids.emplace_back(punchout_group);
        std::string name = "Frill";
        for(std::string const & pmt_id : group) {
            name += "_" + pmt_id;
            std::tuple<std::array<std::shared_ptr<G4VSolid>, 5>, std::shared_ptr<G4RotationMatrix>, G4ThreeVector> const & base_solid = base_solids[pmt_id];
            std::shared_ptr<G4VSolid> frill = std::get<0>(base_solid)[2];
            std::shared_ptr<G4VSolid> bridle = std::get<0>(base_solid)[3];
            G4RotationMatrix * rotationMatrix = std::get<1>(base_solid).get();
            G4ThreeVector position = std::get<2>(base_solid);
            //G4Transform3D transform = G4Transform3D(*rotationMatrix, position);
            G4Transform3D transform;
            G4DisplacedSolid * displaced_frill = new G4DisplacedSolid("Frill" + group_id, frill.get(), rotationMatrix, position);
            temp_solids.emplace_back(displaced_frill);
            G4DisplacedSolid * displaced_bridle = new G4DisplacedSolid("Bridle" + group_id, bridle.get(), rotationMatrix, position);
            temp_solids.emplace_back(displaced_bridle);
            frill_group->AddNode(*displaced_frill, transform);
            punchout_group->AddNode(*displaced_bridle, transform);
        }
        frill_group->Voxelize();
        punchout_group->Voxelize();
        std::shared_ptr<G4VSolid> frill_group_punchout(new G4SubtractionSolid("Frill Group Punchout " + group_id, frill_group, punchout_group));
        temp_solids.emplace_back(frill_group_punchout);
        frill_solids[name] = frill_group_punchout;
        frill_rotations[name] = std::make_shared<G4RotationMatrix>();
        frill_positions[name] = G4ThreeVector(0, 0, 0);
    }
    for(std::pair<unsigned int const, std::vector<std::string>> const & overlap_group : bridle_overlap_groups) {
        std::vector<std::string> const & group = overlap_group.second;
        std::string group_id = std::to_string(overlap_group.first);
        if(group.size() == 1) {
            std::string pmt_id = group[0];
            bool pmt_on_cap, coated;
            int pmt_row, pmt_col, pmt_ring, pmt_number;
            ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

            std::string name = "Bridle_" + pmt_id;
            std::shared_ptr<G4VSolid> bridle;
            if(pmt_on_cap)
                if(coated)
                    bridle = bridle_cap_coated_solid;
                else
                    bridle = bridle_cap_uncoated_solid;
            else
                if(coated)
                    bridle = bridle_wall_coated_solid;
                else
                    bridle = bridle_wall_uncoated_solid;
            //std::shared_ptr<G4VSolid> pmt = std::get<0>(base_solids[pmt_id])[1];
            //std::shared_ptr<G4VSolid> tpb = std::get<0>(base_solids[pmt_id])[2];
            //std::shared_ptr<G4VSolid> bridle = std::get<0>(base_solids[pmt_id])[3];
            //std::shared_ptr<G4VSolid> bridle_group_punchout(new G4SubtractionSolid("Bridle Group Punchout " + group_id,
            //        bridle.get(),
            //        tpb ? tpb.get() : pmt.get()));
            bridle_solids[name] = bridle;
            bridle_rotations[name] = std::get<1>(base_solids[pmt_id]);
            bridle_positions[name] = std::get<2>(base_solids[pmt_id]);
            continue;
        }
        G4MultiUnion * bridle_group = new G4MultiUnion("Bridle Group " + group_id);
        temp_solids.emplace_back(bridle_group);
        G4MultiUnion * punchout_group = new G4MultiUnion("Punchout Group " + group_id);
        temp_solids.emplace_back(punchout_group);
        std::string name = "Bridle";
        for(std::string const & pmt_id : group) {
            name += "_" + pmt_id;
            std::tuple<std::array<std::shared_ptr<G4VSolid>, 5>, std::shared_ptr<G4RotationMatrix>, G4ThreeVector> const & base_solid = base_solids[pmt_id];
            std::shared_ptr<G4VSolid> bridle = std::get<0>(base_solid)[0];
            std::shared_ptr<G4VSolid> punchout = std::get<0>(base_solid)[1];
            G4RotationMatrix * rotationMatrix = std::get<1>(base_solid).get();
            G4ThreeVector position = std::get<2>(base_solid);
            //G4Transform3D transform = G4Transform3D(*rotationMatrix, position);
            G4Transform3D transform;
            G4DisplacedSolid * displaced_bridle = new G4DisplacedSolid("Bridle" + group_id, bridle.get(), rotationMatrix, position);
            temp_solids.emplace_back(displaced_bridle);
            G4DisplacedSolid * displaced_punchout = new G4DisplacedSolid("Punchout" + group_id, punchout.get(), rotationMatrix, position);
            temp_solids.emplace_back(displaced_punchout);
            bridle_group->AddNode(*displaced_bridle, transform);
            punchout_group->AddNode(*displaced_punchout, transform);

            //bridle_group->AddNode(*bridle, transform);
            //punchout_group->AddNode(
            //    tpb ? *tpb : *pmt,
            //    transform);
        }
        bridle_group->Voxelize();
        punchout_group->Voxelize();
        std::shared_ptr<G4VSolid> bridle_group_punchout(new G4SubtractionSolid("Bridle Group Punchout " + group_id, bridle_group, punchout_group));
        temp_solids.emplace_back(bridle_group_punchout);
        bridle_solids[name] = bridle_group_punchout;
        bridle_rotations[name] = std::make_shared<G4RotationMatrix>();
        bridle_positions[name] = G4ThreeVector(0, 0, 0);
    }

    G4PhysicalVolumeStore::GetInstance()->DeRegister(motherPV);
    G4LogicalVolumeStore::GetInstance()->DeRegister(motherLV);
    G4SolidStore::GetInstance()->DeRegister(motherSolid);
    delete motherPV;
    delete motherLV;
    delete motherSolid;
}

    std::shared_ptr<G4VSolid> CreateShinySolid();
    std::shared_ptr<G4VSolid> CreateShinyTopSolid();
    std::shared_ptr<G4VSolid> CreateShinyBottomSolid();

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateShinySolid() {
    temp_solids.emplace_back(new G4Tubs("ShinyC406R0", 0*cm, shiny_radius, shiny_half_height, 0*deg, 360*deg));
    shiny_solid = temp_solids.back();

    G4ThreeVector displacement(0, 0, pmt_protrusion_distance - pmt_radius + shiny_half_height);

    std::string pmt_id = "C406R0";
    bool pmt_on_cap, coated;
    int pmt_row, pmt_col, pmt_ring, pmt_number;
    ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

    std::shared_ptr<G4RotationMatrix> rotationMatrix(GetPMTRotationMatrix(pmt_id));
    G4ThreeVector position = GetPMTPosition(pmt_id) + displacement;

    shiny_solids[pmt_id] = shiny_solid;
    shiny_rotations[pmt_id] = rotationMatrix;
    shiny_positions[pmt_id] = position;

    return shiny_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateShinyTopSolid() {
    temp_solids.emplace_back(new G4Tubs("ShinyTopOriginal", 0*cm, shiny_top_radius, shiny_half_height, 0*deg, 360*deg));
    G4VSolid * shiny_disk = temp_solids.back().get();

    G4ThreeVector shiny_top_position(0, 0, cylinder_half_height - shiny_half_height);

    int top_row = 0;

    std::vector<G4VSolid *> subtraction_solids;

    for(std::string const & pmt_id : pmt_ids) {
        bool pmt_on_cap, coated;
        int pmt_row, pmt_col, pmt_ring, pmt_number;
        ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

        if(pmt_row != top_row)
            continue;

        std::shared_ptr<G4RotationMatrix> rotationMatrix(GetPMTRotationMatrix(pmt_id));
        G4ThreeVector position = GetPMTPosition(pmt_id);

        std::shared_ptr<G4VSolid> bridle_solid = bridle_cap_filled_solid;
        std::shared_ptr<G4VSolid> frill_solid = frill_cap_filled_solid;

        G4double dx = shiny_top_position.x() - position.x();
        G4double dy = shiny_top_position.y() - position.y();
        G4double dr = std::sqrt(dx*dx + dy*dy);
        if(dr < shiny_top_radius + bridle_radius) {
            G4DisplacedSolid * displaced_bridle = new G4DisplacedSolid("DisplacedBridle", bridle_solid.get(), rotationMatrix.get(), position);
            temp_solids.emplace_back(displaced_bridle);
            subtraction_solids.push_back(displaced_bridle);
        }
        if(dr < shiny_top_radius + frill_radius) {
            G4DisplacedSolid * displaced_frill = new G4DisplacedSolid("DisplacedFrill", frill_solid.get(), rotationMatrix.get(), position);
            temp_solids.emplace_back(displaced_frill);
            subtraction_solids.push_back(displaced_frill);
        }
    }

    G4MultiUnion * shiny_top_sub_union = new G4MultiUnion("ShinyTopSubtractionUnion");
    temp_solids.emplace_back(shiny_top_sub_union);
    for(G4VSolid * solid : subtraction_solids) {
        shiny_top_sub_union->AddNode(*solid, G4Transform3D());
    }
    shiny_top_sub_union->Voxelize();

    G4SubtractionSolid * shiny_top_subtraction = new G4SubtractionSolid("ShinyTop", shiny_disk, shiny_top_sub_union, 0, -shiny_top_position);
    temp_solids.emplace_back(shiny_top_subtraction);
    shiny_top_solid = temp_solids.back();

    shiny_solids["ShinyTop"] = shiny_top_solid;
    shiny_rotations["ShinyTop"] = nullptr;
    shiny_positions["ShinyTop"] = shiny_top_position;

    return shiny_top_solid;
}

std::shared_ptr<G4VSolid> J4PMTSolidMaker::CreateShinyBottomSolid() {
    temp_solids.emplace_back(new G4Tubs("ShinyBottomOriginal", 0*cm, shiny_top_radius, shiny_half_height, 0*deg, 360*deg));
    G4VSolid * shiny_disk = temp_solids.back().get();

    G4ThreeVector shiny_bottom_position(0, 0, -cylinder_half_height + shiny_half_height);

    int bottom_row = 6;

    std::vector<G4VSolid *> subtraction_solids;

    for(std::string const & pmt_id : pmt_ids) {
        bool pmt_on_cap, coated;
        int pmt_row, pmt_col, pmt_ring, pmt_number;
        ParsePMTID(pmt_id, pmt_on_cap, pmt_row, pmt_col, pmt_ring, pmt_number, coated);

        if(pmt_row != bottom_row)
            continue;

        std::shared_ptr<G4RotationMatrix> rotationMatrix(GetPMTRotationMatrix(pmt_id));
        G4ThreeVector position = GetPMTPosition(pmt_id);

        std::shared_ptr<G4VSolid> bridle_solid = bridle_cap_filled_solid;
        std::shared_ptr<G4VSolid> frill_solid = frill_cap_filled_solid;

        G4double dx = shiny_bottom_position.x() - position.x();
        G4double dy = shiny_bottom_position.y() - position.y();
        G4double dr = std::sqrt(dx*dx + dy*dy);
        if(dr < shiny_top_radius + bridle_radius) {
            G4DisplacedSolid * displaced_bridle = new G4DisplacedSolid("DisplacedBridle", bridle_solid.get(), rotationMatrix.get(), position);
            temp_solids.emplace_back(displaced_bridle);
            subtraction_solids.push_back(displaced_bridle);
        }
        if(dr < shiny_top_radius + frill_radius) {
            G4DisplacedSolid * displaced_frill = new G4DisplacedSolid("DisplacedFrill", frill_solid.get(), rotationMatrix.get(), position);
            temp_solids.emplace_back(displaced_frill);
            subtraction_solids.push_back(displaced_frill);
        }
    }

    G4MultiUnion * shiny_bottom_sub_union = new G4MultiUnion("ShinyBottomSubtractionUnion");
    temp_solids.emplace_back(shiny_bottom_sub_union);
    for(G4VSolid * solid : subtraction_solids) {
        shiny_bottom_sub_union->AddNode(*solid, G4Transform3D());
    }
    shiny_bottom_sub_union->Voxelize();

    G4SubtractionSolid * shiny_bottom_subtraction = new G4SubtractionSolid("ShinyBottom", shiny_disk, shiny_bottom_sub_union, 0, -shiny_bottom_position);
    temp_solids.emplace_back(shiny_bottom_subtraction);
    shiny_bottom_solid = temp_solids.back();

    shiny_solids["ShinyBottom"] = shiny_bottom_solid;
    shiny_rotations["ShinyBottom"] = nullptr;
    shiny_positions["ShinyBottom"] = shiny_bottom_position;

    return shiny_bottom_solid;
}

void J4PMTSolidMaker::CreateShinySolids() {
    shiny_solid = CreateShinySolid();
    shiny_top_solid = CreateShinyTopSolid();
    shiny_bottom_solid = CreateShinyBottomSolid();
}


/*
void J4PMTSolidMaker::CreateFrillWallSolid(G4double bridle_radius, G4double frill_radius, G4double frill_width, G4double tpb_foil_radius, G4bool coated) {

    G4double max_diameter = 250 * mm;
    G4VSolid * fFrillPunchOut = new G4Tubs("FrillPunchOut", bridle_radius, frill_radius, max_diameter/2.0, 0*deg, 360*deg);

    G4VSolid * fFrillSurface = new G4Tubs("FrillSurface", tpb_foil_radius - frill_width, tpb_foil_radius, frill_radius * 2, 0*deg, 360*deg);

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(M_PI / 2.0); // Rotate 90degrees around the Y axis

    G4double z_center = tpb_foil_radius;
    G4ThreeVector centerOfTub(0, 0, z_center);

    if(coated) {
        fFrillWallCoated = new G4IntersectionSolid("FrillWall", fFrillPunchOut, fFrillSurface, rotationMatrix, centerOfTub);
    } else {
        fFrillWallUncoated = new G4IntersectionSolid("FrillWall", fFrillPunchOut, fFrillSurface, rotationMatrix, centerOfTub);
    }
}

void J4PMTSolidMaker::CreateFrillCapSolid(G4double bridle_radius, G4double frill_radius, G4double frill_width, G4bool coated) {
    if(coated) {
        fFrillCapCoated = new G4Tubs("FrillCap", bridle_radius, frill_radius, frill_width/2.0, 0*deg, 360*deg);
    } else {
        fFrillCasUncoated = new G4Tubs("FrillCap", bridle_radius, frill_radius, frill_width/2.0, 0*deg, 360*deg);
    }
}


void J4PMTSolidMaker::CreateBridleWallSolid(G4double bridle_radius, G4double bridle_width, G4double tpb_foil_radius, G4double protrusion_distance, G4bool coated) {
    G4double max_diameter = 250 * mm;
    G4VSolid * fBridlePunchOut = new G4Tubs("BridlePunchOut", 0*cm, bridle_radius, max_diameter/2.0, 0*deg, 360*deg);

    G4VSolid * pmt;
    G4double radius = f8inchPMTRadius;
    if (coated){
        //GetTPBCoatingSolid();
        //pmt = fPMTAndTPBCoatingSolid;
        //radius += fTPBThickness/mm;
        pmt = Get8inchPMTSolid();
    } else {
        pmt = Get8inchPMTSolid();
    }

    G4VSolid * fBridleSurface = new G4Tubs("BridleSurface", tpb_foil_radius - bridle_width, tpb_foil_radius, bridle_radius * 2, 0*deg, 360*deg);

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(M_PI / 2.0); // Rotate 90degrees around the Y axis

    G4double z_center = tpb_foil_radius;
    G4ThreeVector centerOfTub(0, 0, z_center);

    G4VSolid * fBridleSurfaceBridlePunchOut = new G4IntersectionSolid("BridlePunchOut", fBridlePunchOut, fBridleSurface, rotationMatrix, centerOfTub);

    //G4double z_center_bridle = f8inchPMTRadius - protrusion_distance;
    G4double z_center_bridle = radius - protrusion_distance;
    G4ThreeVector centerofPMT(0, 0, z_center_bridle);

    if(coated) {
        fBridleWallCoated = new G4SubtractionSolid("Bridle", fBridleSurfaceBridlePunchOut, pmt, 0, centerofPMT);
    } else {
        fBridleWallUncoated = new G4SubtractionSolid("Bridle", fBridleSurfaceBridlePunchOut, pmt, 0, centerofPMT);
    }
}

void J4PMTSolidMaker::CreateBridleCapSolid(G4double bridle_radius, G4double bridle_width, G4double protrusion_distance, G4bool coated) {
    G4VSolid * pmt;
    G4double radius = f8inchPMTRadius;
    if (coated){
        //GetTPBCoatingSolid();
        //pmt = fPMTAndTPBCoatingSolid;
        //radius += fTPBThickness/mm;
        pmt = Get8inchPMTSolid();
    } else {
        pmt = Get8inchPMTSolid();
    }

    G4VSolid * fBridleDisk = new G4Tubs("BridleDisk", 0*cm, bridle_radius, bridle_width/2.0, 0*deg, 360*deg);

    //G4double z_offset = -(f8inchPMTRadius - protrusion_distance + bridle_width / 2.0);
    G4double z_offset = -(radius - protrusion_distance + bridle_width / 2.0);
    G4ThreeVector centerofBridlePMT(0, 0, z_offset);

    if(coated) {
        fBridleCapCoated = new G4SubtractionSolid("BridleCaps", fBridleDisk, pmt, 0, centerofBridlePMT);
    } else {
        fBridleCapUncoated = new G4SubtractionSolid("BridleCaps", fBridleDisk, pmt, 0, centerofBridlePMT);
    }

}

void J4PMTSolidMaker::Create8inchPMTCapSolid(G4double protrusion_distance) {
    G4VSolid * pmt = Get8inchPMTSolid();

    G4double max_diameter = 250 * mm;
    G4double max_height = 290 * mm;

    fBoxSolid = new G4Box("box", max_diameter/2.0, max_diameter/2.0, max_height/2.0);

    G4double z_bottom = f8inchPMTRadius - protrusion_distance;
    G4double z_box_center = z_bottom + max_height/2.0;

    G4ThreeVector centerOfBox(0, 0, z_box_center);

    f8inchPMTCapSolid = new G4IntersectionSolid("8inchPMTCaps", pmt, fBoxSolid, 0, centerOfBox);
}

void J4PMTSolidMaker::Create8inchPMTWallSolid(G4double cylinder_radius, G4double protrusion_distance) {
    G4VSolid * pmt = Get8inchPMTSolid();

    G4double max_diameter = 250 * mm;

    fTubSolid = new G4Tubs("InnerFrame", 0*cm, cylinder_radius, max_diameter/2.0, 0*deg, 360*deg);

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(M_PI / 2.0); // Rotate 90degrees around the x axis

    G4double z_center = f8inchPMTRadius - protrusion_distance + cylinder_radius;

    G4ThreeVector centerOfTub(0, 0, z_center);

    f8inchPMTWallSolid = new G4IntersectionSolid("8inchPMTWall", pmt, fTubSolid, rotationMatrix, centerOfTub);
}

void J4PMTSolidMaker::CreateTPBCoatingCapSolid(G4double protrusion_distance) {
    G4VSolid * tpb = GetTPBCoatingSolid();

    G4double max_diameter = 250 * mm;
    G4double max_height = 290 * mm;

    fBoxSolid = new G4Box("box", max_diameter/2.0, max_diameter/2.0, max_height/2.0);

    G4double z_bottom = f8inchPMTRadius - protrusion_distance;
    G4double z_box_center = z_bottom + max_height/2.0;

    G4ThreeVector centerOfBox(0, 0, z_box_center);

    fTPBCoatingCapSolid = new G4IntersectionSolid("TPBCoatingCaps", tpb, fBoxSolid, 0, centerOfBox);
}

void J4PMTSolidMaker::CreateTPBCoatingWallSolid(G4double cylinder_radius, G4double protrusion_distance) {
    G4VSolid * tpb = GetTPBCoatingSolid();

    G4double max_diameter = 250 * mm;

    fTubSolid = new G4Tubs("InnerFrame", 0*cm, cylinder_radius, max_diameter/2.0, 0*deg, 360*deg);

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateY(M_PI / 2.0); // Rotate 90degrees around the x axis

    G4double z_center = f8inchPMTRadius - protrusion_distance + cylinder_radius;

    G4ThreeVector centerOfTub(0, 0, z_center);

    fTPBCoatingWallSolid = new G4IntersectionSolid("TPBCoatingWall", tpb, fTubSolid, rotationMatrix, centerOfTub);
}

///////////////////////////////////////////////////////

void J4PMTSolidMaker::CreatePhotocathodeSolid() {

    // so we need to make photocathode that is smaller than the pmts by an even amounts
    // to do this, we modify the paramters as such:
    // sr0 --> sr0 - delta_photocathode
    // rSmax --> rSmax - delta_photocathode
    // rmax --> rmax - delta_photocathode*1.75
    // rmax2 --> rmax2 - delta_photocathode*1.75
    // centerofConslb --> centerofConslb + delta_photocathode/2.1 (centerofConslb is 70.4284)

    G4ThreeVector centerOfPolycone;
    G4ThreeVector centerOfCons;
    G4ThreeVector centerOfTubs;

    G4double rmin2, rmax2, dz;
    G4double rmin, rmax, sphi, dphi, stheta, dtheta;
    G4double rSmin, rSmax, sSphi, dSphi, sStheta, dStheta;
    G4double centerofConslb;
    centerofConslb = 70.4284 + fPhotocathodeDelta / 2.1;

    rSmin   = 0.0;
    rSmax   = f8inchPMTRadius - fPhotocathodeDelta;
    dSphi   = 2*M_PI;
    dStheta = 37.6392 *degree;

    std::vector<G4double> tempZ, tempInner, tempOuter;

    G4int segment = 100;

    G4double sr0 = 57.0894 - fPhotocathodeDelta;//radius of spindle torus sphere
    G4double centerOfsr0 = 43.9106; // distance from center of torus sphere to z-axis
    G4double rplanet = sr0*2./segment; // z length of each planet

    //our function is a fixed number+the value of the sphere projection
    //should spread across 2R

    G4double rInner[segment+1], rOuter[segment+1], zPlane[segment+1];

    for (G4int j=0; j<=segment; ++j) {
      tempZ.push_back((sr0 - j*rplanet)*mm);
      tempInner.push_back(0.);
      tempOuter.push_back((centerOfsr0 + sqrt(sr0*sr0-(sr0-j*rplanet)*(sr0-j*rplanet)))*mm);
    }

    for (G4int i=0; i<=segment; i++) {
      rInner[i] = tempInner[i];
      rOuter[i] = tempOuter[i];
      zPlane[i] = tempZ[i];
    }

    G4Polycone* polycone1 = new
      G4Polycone("polycone1", 0, 2*M_PI, segment+1, zPlane, rInner, rOuter);

    centerOfPolycone = G4ThreeVector(0, 0, 59.5 *mm);
    centerOfCons = G4ThreeVector(0, 0, (centerofConslb+33.8363/2-89.) *mm);
    centerOfTubs = G4ThreeVector(0, 0, -(89.-centerofConslb/2) *mm);

    sSphi = 0;
    sStheta = 0;
    G4Sphere *sphere = new G4Sphere("sphere",rSmin, rSmax,
				      sSphi, dSphi, sStheta, dStheta);

    //to create two cons

    rmin   = 0;
    rmin2  = 0;
    rmax   = 84.5/2*mm - (1.75*fPhotocathodeDelta)/2;
    rmax2  = 160./2 *mm - (1.75*fPhotocathodeDelta)/2;
    dz     = 33.8363/2 *mm;
    sphi   = 0;
    dphi   = 2*M_PI;

    G4Cons *cons= new G4Cons("cons",rmin, rmax, rmin2,
			       rmax2, dz, sphi, dphi);

    // create two tubs....

    rmin   = 0;
    rmax   = 84.5/2 *mm;
    dz     = centerofConslb/2 *mm;
    sphi   = 0;
    dphi   = 2*M_PI;

    G4Tubs *tubs = new G4Tubs("tubs",rmin, rmax, dz,
				sphi, dphi);

    // to create two PMTs

    G4UnionSolid *solid1 = new G4UnionSolid("solid1", sphere, polycone1, 0, centerOfPolycone);
    fPhotocathodeSolid = solid1;
    //G4UnionSolid *solid2 = new G4UnionSolid("solid2", solid1, tubs, 0, centerOfTubs);

    //fPhotocathodeSolid = new G4UnionSolid("solid", solid2, cons, 0, centerOfCons);
}

*/
