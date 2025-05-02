// $Id: J4PMTSolidMaker.hh,v 1.1.1.1 2004/08/26 07:04:26 hoshina Exp $
#ifndef __J4PMTSolidMaker__hh
#define __J4PMTSolidMaker__hh
//*************************************************************************
//* --------------------
//* J4PMTSolidMaker
//* --------------------
//* (Description)
//* 	J4PMTSolidMaker discribes how to make glass solid for DEgg
//*
//* (Update Record)
//*	2021/05/27  K.Hoshina	Original version.
//*************************************************************************

#include <array>
#include <memory>
#include <vector>

#include <G4VSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4RotationMatrix.hh>

class J4PMTSolidMaker {
    std::vector<G4VSolid *> temp_solids;
public:
    G4VSolid * pmt_vacuum_filled_solid = nullptr;
    G4VSolid * pmt_glass_filled_solid = nullptr;
    G4VSolid * pmt_tpb_filled_solid = nullptr;
    G4VSolid * pmt_wall_vacuum_filled_solid = nullptr;
    G4VSolid * pmt_wall_vacuum_filled_offset_solid = nullptr;
    G4VSolid * pmt_wall_glass_filled_solid = nullptr;
    G4VSolid * pmt_wall_glass_filled_offset_solid = nullptr;
    G4VSolid * pmt_wall_tpb_filled_solid = nullptr;
    G4VSolid * pmt_cap_vacuum_filled_solid = nullptr;
    G4VSolid * pmt_cap_vacuum_filled_offset_solid = nullptr;
    G4VSolid * pmt_cap_glass_filled_solid = nullptr;
    G4VSolid * pmt_cap_glass_filled_offset_solid = nullptr;
    G4VSolid * pmt_cap_tpb_filled_solid = nullptr;

    G4VSolid * cylinder_curvature_solid = nullptr;
    G4VSolid * cylinder_face_solid = nullptr;
    G4RotationMatrix * cylinder_curvature_rotation = nullptr;
    G4ThreeVector cylinder_curvature_position;
    G4RotationMatrix * cylinder_face_rotation = nullptr;
    G4ThreeVector cylinder_face_position;

    G4VSolid * frill_center_solid = nullptr;

    G4VSolid * bridle_filled_solid = nullptr;
    G4VSolid * bridle_cylinder_solid = nullptr;
    G4VSolid * frill_filled_solid = nullptr;
    G4VSolid * frill_cylinder_solid = nullptr;

    G4VSolid * bridle_wall_filled_solid = nullptr;
    G4VSolid * bridle_cap_filled_solid = nullptr;
    G4VSolid * frill_wall_filled_solid = nullptr;
    G4VSolid * frill_cap_filled_solid = nullptr;

    G4VSolid * pmt_glass_solid = nullptr;
    G4VSolid * pmt_tpb_solid = nullptr;
    G4VSolid * pmt_wall_glass_solid = nullptr;
    G4VSolid * pmt_cap_glass_solid = nullptr;
    G4VSolid * pmt_wall_tpb_solid = nullptr;
    G4VSolid * pmt_cap_tpb_solid = nullptr;

    G4VSolid * bridle_wall_coated_solid = nullptr;
    G4VSolid * bridle_wall_uncoated_solid = nullptr;
    G4VSolid * bridle_cap_coated_solid = nullptr;
    G4VSolid * bridle_cap_uncoated_solid = nullptr;

    G4VSolid * frill_wall_solid = nullptr;
    G4VSolid * frill_cap_solid = nullptr;

    G4VSolid * shiny_solid = nullptr;
    G4VSolid * shiny_top_solid = nullptr;
    G4VSolid * shiny_bottom_solid = nullptr;

    std::map<std::string, std::tuple<std::array<G4VSolid *, 5>, std::shared_ptr<G4RotationMatrix>, G4ThreeVector>> base_solids;

    std::map<std::string, G4VSolid *> frill_solids;
    std::map<std::string, G4ThreeVector> frill_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> frill_rotations;

    std::map<std::string, G4VSolid *> bridle_solids;
    std::map<std::string, G4ThreeVector> bridle_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> bridle_rotations;

    std::map<std::string, G4VSolid *> tpb_solids;
    std::map<std::string, G4ThreeVector> tpb_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> tpb_rotations;

    std::map<std::string, G4VSolid *> pmt_solids;
    std::map<std::string, G4ThreeVector> pmt_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> pmt_rotations;

    std::map<std::string, G4VSolid *> vacuum_solids;
    std::map<std::string, G4ThreeVector> vacuum_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> vacuum_rotations;

    std::map<std::string, G4VSolid *> shiny_solids;
    std::map<std::string, G4ThreeVector> shiny_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> shiny_rotations;

    G4double cylinder_radius = 1037.0*mm - 3.15*mm;
    G4double cylinder_half_height = 1239.6*mm / 2.0;
    G4double pmt_radius = 131.0 * mm;
    G4double pmt_protrusion_distance = 61.89 * mm;
    G4double pmt_tpb_thickness = 0.00203892 * mm;;
    G4double pmt_glass_thickness = 1.0 * mm;

    bool solid_pmt = false;
    bool protruding_pmt = false;

    G4double bridle_width = 12.0 * mm;
    G4double bridle_radius = 209.55 / 2.0 * mm;
    G4double frill_width = 0.813 * mm;
    G4double frill_radius = 228.6 / 2.0 * mm;
    G4double shiny_top_radius = 13.335 * cm;
    G4double shiny_radius = 205.0 / 2.0 * mm; // making the assumption that it is the same radius as the bridle
    G4double shiny_half_height = 0.125 * mm;

    G4double glass_offset = 1.0 * um;
    G4double vacuum_offset = 2.0 * um;

    J4PMTSolidMaker() {}

    G4double pmt_radial_distance;
    G4double pmt_vertical_distance;
    G4double frill_radial_distance;
    G4double frill_vertical_distance;
    G4double bridle_radial_distance;
    G4double bridle_vertical_distance;

    void ComputeParameters() {
        pmt_radial_distance = cylinder_radius + (pmt_radius - pmt_protrusion_distance);
        pmt_vertical_distance = cylinder_half_height + (pmt_radius - pmt_protrusion_distance);

        frill_radial_distance =  cylinder_radius;
        frill_vertical_distance = cylinder_half_height - (frill_width / 2.0);

        bridle_radial_distance = cylinder_radius;
        bridle_vertical_distance = cylinder_half_height - (bridle_width / 2.0);
    }

    G4VSolid * CreatePMTPolyconeSolid(double scale_factor);
    G4VSolid * CreatePMTSphereSolid(double scale_factor);
    G4VSolid * CreatePMTConeSolid(double scale_factor);
    G4VSolid * CreatePMTTubsSolid(double scale_factor);
    G4VSolid * CreatePMTSolid(double scale_factor);

    G4VSolid * CreatePMTGlassFilledSolid();
    G4VSolid * CreatePMTTPBFilledSolid();
    G4VSolid * CreatePMTVacuumFilledSolid();

    G4VSolid * CreateCylinderCurvatureSolid();
    G4VSolid * CreateCylinderFaceSolid();

    G4VSolid * CreateWallSolid(G4VSolid * pmt_solid, G4double offset = 0.0);
    G4VSolid * CreateCapSolid(G4VSolid * pmt_solid, G4double offset = 0.0);

    G4VSolid * CreatePMTWallGlassFilledSolid();
    G4VSolid * CreatePMTWallTPBFilledSolid();
    G4VSolid * CreatePMTWallVacuumFilledSolid();
    G4VSolid * CreatePMTCapGlassFilledSolid();
    G4VSolid * CreatePMTCapTPBFilledSolid();
    G4VSolid * CreatePMTCapVacuumFilledSolid();

    G4VSolid * CreateFrillCenterSolid();

    G4VSolid * CreateBridleFilledSolid();
    G4VSolid * CreateBridleCylinderSolid();
    G4VSolid * CreateFrillFilledSolid();
    G4VSolid * CreateFrillCylinderSolid();

    G4VSolid * CreateBridleWallFilledSolid();
    G4VSolid * CreateBridleCapFilledSolid();
    G4VSolid * CreateFrillWallFilledSolid();
    G4VSolid * CreateFrillCapFilledSolid();

    G4VSolid * CreatePMTGlassSolid();
    G4VSolid * CreatePMTTPBSolid();
    G4VSolid * CreatePMTWallGlassSolid();
    G4VSolid * CreatePMTCapGlassSolid();
    G4VSolid * CreatePMTWallTPBSolid();
    G4VSolid * CreatePMTCapTPBSolid();

    G4VSolid * CreateBridleWallCoatedSolid();
    G4VSolid * CreateBridleWallUncoatedSolid();
    G4VSolid * CreateBridleCapCoatedSolid();
    G4VSolid * CreateBridleCapUncoatedSolid();
    G4VSolid * CreateFrillWallSolid();
    G4VSolid * CreateFrillCapSolid();

    G4VSolid * CreateShinySolid();
    G4VSolid * CreateShinyTopSolid();
    G4VSolid * CreateShinyBottomSolid();

    G4RotationMatrix * GetPMTRotationMatrix(std::string pmt_id);
    G4ThreeVector GetPMTPosition(std::string pmt_id);

    static void ParsePMTID(std::string pmt_id, bool & pmt_on_cap, int & pmt_row, int & pmt_col, int & pmt_ring, int & pmt_number, bool & coated);
    void CreateBaseSolids();
    void CreateShinySolids();

    void CreateTPBCoatingSolid();
    void CreatePhotocathodeSolid();

  private:
};

#endif


