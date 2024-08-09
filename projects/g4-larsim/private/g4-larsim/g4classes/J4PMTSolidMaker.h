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

#include <G4VSolid.hh>

class J4PMTSolidMaker {

  protected :
    J4PMTSolidMaker() {}

  public:

    static void Create8inchPMTSolid();
    static void CreateTPBCoatingSolid();
    static void CreatePhotocathodeSolid();

    static void Create8inchPMTCapsSolid(G4double protrusion_distance);
    static void Create8inchPMTWallSolid(G4double cylinder_radius, G4double protrusion_distance);
    static void CreateTPBCoatingCapsSolid(G4double protrusion_distance);
    static void CreateTPBCoatingWallSolid(G4double cylinder_radius, G4double protrusion_distance);
    static void CreateBridleWallSolid(G4double bridle_radius, G4double bridle_width, G4double tpb_foil_radius, G4double protrusion_distance);
    static void CreateBridleCapsSolid(G4double bridle_radius, G4double bridle_width, G4double protrusion_distance);
    static void CreateFrillWallSolid(G4double bridle_radius, G4double frill_radius, G4double frill_width, G4double tpb_foil_radius);
    static void CreateFrillCapsSolid(G4double bride_radius, G4double frill_radius, G4double frill_width);

    static G4VSolid * Get8inchPMTSolid();
    static G4VSolid * GetTPBCoatingSolid();
    static G4VSolid * GetPhotocathodeSolid();

    static G4double Get8inchPMTRadius() { return f8inchPMTRadius; }
    static G4double GetTPBThickness() { return fTPBThickness; }
    static G4double GetPhotocathodeThickness() { return fPhotocathodeDelta; }

    static G4VSolid * Get8inchPMTCapsSolid(G4double protrusion_distance);
    static G4VSolid * Get8inchPMTWallSolid(G4double cylinder_radius, G4double protrusion_distance);
    static G4VSolid * GetTPBCoatingCapsSolid(G4double protrusion_distance);
    static G4VSolid * GetTPBCoatingWallSolid(G4double cylinder_radius, G4double protrusion_distance);
    static G4VSolid * GetBridleWall(G4double bridle_radius, G4double bridle_width, G4double tpb_foil_radius, G4double protrusion_distance);
    static G4VSolid * GetBridleCaps(G4double bridle_radius, G4double bridle_width, G4double protrusion_distance);
    static G4VSolid * GetFrillWall(G4double bridle_radius, G4double frill_radius,  G4double frill_width, G4double tpb_foil_radius);
    static G4VSolid * GetFrillCaps(G4double bridle_radius, G4double frill_radius, G4double frill_width);

  private:
    static G4double f8inchPMTRadius;
    static G4double fTPBThickness;
    static G4double fCylinderRadius;
    static G4double fPMTProtrusionDistance;

    static G4VSolid * fTubSolid;
    static G4VSolid * fBoxSolid;

    static G4VSolid * f8inchPMTSolid;

    static G4VSolid * f8inchPMTCapsSolid;
    static G4VSolid * f8inchPMTWallSolid;

    static G4VSolid * fTPBCoatingSolid;
    static G4VSolid * fTPBCoatingCapsSolid;
    static G4VSolid * fTPBCoatingWallSolid;

    static G4VSolid * fBridleWall;
    static G4VSolid * fBridleCaps;
    static G4VSolid * fFrillWall;
    static G4VSolid * fFrillCaps;
    
    static G4VSolid * fPhotocathodeSolid;
    static G4double fPhotocathodeDelta;

};

#endif


