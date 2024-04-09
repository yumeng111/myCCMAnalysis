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

    static G4VSolid * Get8inchPMTSolid();
    static G4VSolid * GetTPBCoatingSolid();
    static G4VSolid * GetPhotcathodeSolid();

    static G4double Get8inchPMTRadius() { return f8inchPMTRadius; }
    static G4double GetTPBThickness() { return fTPBThickness; }
    static G4double GetPhotocathodeThickness() { return fPhotocathodeDelta; }

  private:

    static G4VSolid * f8inchPMTSolid;
    static G4VSolid * fTPBCoatingSolid;
    static G4VSolid * fPhotocathodeSolid;
    
    static G4double f8inchPMTRadius;
    static G4double fTPBThickness;
    static G4double fPhotocathodeDelta;

};

#endif


