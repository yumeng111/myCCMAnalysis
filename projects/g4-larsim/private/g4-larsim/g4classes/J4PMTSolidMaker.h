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
#include <vector>

#include <G4VSolid.hh>
#include <G4SystemOfUnits.hh>
#include <G4RotationMatrix.hh>

class J4PMTSolidMaker {
    std::vector<std::shared_ptr<G4VSolid>> temp_solids;
public:
    std::shared_ptr<G4VSolid> pmt_vacuum_solid;
    std::shared_ptr<G4VSolid> pmt_glass_solid;
    std::shared_ptr<G4VSolid> pmt_tpb_solid;
    std::shared_ptr<G4VSolid> pmt_wall_vacuum_solid;
    std::shared_ptr<G4VSolid> pmt_wall_glass_solid;
    std::shared_ptr<G4VSolid> pmt_wall_tpb_solid;
    std::shared_ptr<G4VSolid> pmt_cap_vacuum_solid;
    std::shared_ptr<G4VSolid> pmt_cap_glass_solid;
    std::shared_ptr<G4VSolid> pmt_cap_tpb_solid;

    std::shared_ptr<G4VSolid> bridle_wall_solid;
    std::shared_ptr<G4VSolid> bridle_cap_solid;
    std::shared_ptr<G4VSolid> frill_wall_solid;
    std::shared_ptr<G4VSolid> frill_cap_solid;

    std::map<std::string, std::tuple<std::array<std::shared_ptr<G4VSolid>, 5>, std::shared_ptr<G4RotationMatrix>, G4ThreeVector>> base_solids;

    std::map<std::string, std::shared_ptr<G4VSolid>> frill_solids;
    std::map<std::string, G4ThreeVector> frill_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> frill_rotations;

    std::map<std::string, std::shared_ptr<G4VSolid>> bridle_solids;
    std::map<std::string, G4ThreeVector> bridle_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> bridle_rotations;

    std::map<std::string, std::shared_ptr<G4VSolid>> tpb_solids;
    std::map<std::string, G4ThreeVector> tpb_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> tpb_rotations;

    std::map<std::string, std::shared_ptr<G4VSolid>> pmt_solids;
    std::map<std::string, G4ThreeVector> pmt_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> pmt_rotations;

    std::map<std::string, std::shared_ptr<G4VSolid>> vacuum_solids;
    std::map<std::string, G4ThreeVector> vacuum_positions;
    std::map<std::string, std::shared_ptr<G4RotationMatrix>> vacuum_rotations;

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

    std::shared_ptr<G4VSolid> CreatePMTPolyconeSolid(double scale_factor);
    std::shared_ptr<G4VSolid> CreatePMTSphereSolid(double scale_factor);
    std::shared_ptr<G4VSolid> CreatePMTConeSolid(double scale_factor);
    std::shared_ptr<G4VSolid> CreatePMTTubsSolid(double scale_factor);
    std::shared_ptr<G4VSolid> CreatePMTSolid(double scale_factor);

    std::shared_ptr<G4VSolid> CreatePMTGlassFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTTPBFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTVacuumFilledSolid();

    std::shared_ptr<G4VSolid> CreateWallSolid(std::shared_ptr<G4VSolid> pmt_solid);
    std::shared_ptr<G4VSolid> CreateCapSolid(std::shared_ptr<G4VSolid> pmt_solid);

    std::shared_ptr<G4VSolid> CreatePMTWallGlassFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTWallTPBFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTWallVacuumFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTCapGlassFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTCapTPBFilledSolid();
    std::shared_ptr<G4VSolid> CreatePMTCapVacuumFilledSolid();

    std::shared_ptr<G4VSolid> CreateBridleWallFilledSolid();
    std::shared_ptr<G4VSolid> CreateBridleCapFilledSolid();
    std::shared_ptr<G4VSolid> CreateFrillWallFilledSolid();
    std::shared_ptr<G4VSolid> CreateFrillCapFilledSolid();

    G4RotationMatrix * GetPMTRotationMatrix(std::string pmt_id);
    G4ThreeVector GetPMTPosition(std::string pmt_id);

    static void ParsePMTID(std::string pmt_id, bool & pmt_on_cap, int & pmt_row, int & pmt_col, int & pmt_ring, int & pmt_number, bool & coated);
    void CreateBaseSolids();

    void CreateTPBCoatingSolid();
    void CreatePhotocathodeSolid();

  private:
};

#endif


