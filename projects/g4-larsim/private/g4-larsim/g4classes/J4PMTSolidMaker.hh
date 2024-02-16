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

#include "G4VSolid.hh"

class J4PMTSolidMaker {

  protected :
    J4PMTSolidMaker() {}

  public:


    static void Create8inchPMTSolid();
    static G4VSolid * Get8inchPMTSolid();
    static G4double Get8inchPMTRadius() { return f8inchPMTRadius; }

  private:

    static G4VSolid * f8inchPMTSolid;

    static G4double f8inchPMTRadius;

};

#endif


