// $Id: J4PMTSolidMaker.cc,v 1.2 2007/03/13 17:29:14 hoshina Exp $
/**
* @file J4PMTSolidMaker.cc
* @brief namespace that provides geometory functions to generate G4Solids
* @date 2021/04/11
* @note many numbers are hard coded, if you want to change it on the fly move these parameters to J4PartsParameterList.
* (Update Record)
*	2021/04/11  K.Hoshina	collected functions from each detector classes
*/

#include "g4-larsim/g4classes/J4PMTSolidMaker.h"
#include "g4-larsim/g4classes/J4UnionSolid.h"

#include "G4Sphere.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4UserLimits.hh"
#include "G4Polycone.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"

G4VSolid * J4PMTSolidMaker::f8inchPMTSolid;
G4double   J4PMTSolidMaker::f8inchPMTRadius = 131.0 * mm;

//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* J4PMTSolidMaker::Get8inchPMTSolid()
{
   if (!f8inchPMTSolid) {
       Create8inchPMTSolid();
   }

   return f8inchPMTSolid;
}

//////////////////////////////////////////////////////////////////////////////////////
void J4PMTSolidMaker::Create8inchPMTSolid()
{

    G4ThreeVector centerOfPolycone;
    G4ThreeVector centerOfCons; 
    G4ThreeVector centerOfTubs; 
  
    G4double rmin2, rmax2, dz;
    G4double rmin, rmax, sphi, dphi, stheta, dtheta;
    G4double rSmin, rSmax, sSphi, dSphi, sStheta, dStheta;
  
    rSmin   = 0;
    rSmax   = f8inchPMTRadius;
    dSphi   = 2*M_PI;
    dStheta = 37.6392 *degree;
     
    std::vector<G4double> tempZ, tempInner, tempOuter;

    G4int segment = 100;
  
    G4double sr0 = 57.0894;//radius of spindle torus sphere
    G4double centerOfsr0 = 43.9106; // distance from center of torus sphere to z-axis      
    G4double rplanet = sr0*2./segment; // z length of each planet
  
    //our function is a fixed number+the value of the sphere projection
    //should spread across 2R

    G4double rInner[segment+1], rOuter[segment+1], zPlane[segment+1];

    for (G4int j=0; j<=segment; ++j) {
      tempZ.push_back((sr0 - j*rplanet)*mm);
      tempInner.push_back(0.);
      tempOuter.push_back((centerOfsr0 + sqrt(sr0*sr0-(sr0-j*rplanet)*(sr0-j*rplanet)))*mm);
      //std::cout<< "pmt " <<j<<"  ,  "<<sr0-j*rplanet<<", "<<(centerOfsr0+sqrt(sr0*sr0-(sr0-j*rplanet)*(sr0-j*rplanet)))<<std::endl;
    }
    
    for (G4int i=0; i<=segment; i++) {
      rInner[i] = tempInner[i];
      rOuter[i] = tempOuter[i];
      zPlane[i] = tempZ[i];
    }
    
    G4Polycone* polycone1 = new 
      G4Polycone("polycone1", 0, 2*M_PI, segment+1, zPlane, rInner, rOuter);
    
    std::cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxpositive"<<std::endl;
    centerOfPolycone = G4ThreeVector(0, 0, 59.5 *mm);
    centerOfCons = G4ThreeVector(0, 0, (70.4284+33.8363/2-89.) *mm); 
    centerOfTubs = G4ThreeVector(0, 0, -(89.-70.4284/2) *mm); 

    sSphi = 0;
    sStheta = 0;
    G4Sphere *sphere = new G4Sphere("sphere",rSmin, rSmax, 
				      sSphi, dSphi, sStheta, dStheta);
      
    //to create two cons
     
    rmin   = 0;
    rmin2  = 0;
    rmax   = 84.5/2 *mm;
    rmax2  = 160./2 *mm;  
    dz     = 33.8363/2 *mm;
    sphi   = 0;
    dphi   = 2*M_PI;
      
    G4Cons *cons= new G4Cons("cons",rmin, rmax, rmin2, 
			       rmax2, dz, sphi, dphi);
      
    // create two tubs....
      
    rmin   = 0;
    rmax   = 84.5/2 *mm;
    dz     = 70.4284/2 *mm;
    sphi   = 0;
    dphi   = 2*M_PI;
      
    G4Tubs *tubs = new G4Tubs("tubs",rmin, rmax, dz, 
				sphi, dphi);
      
    // to create two PMTs
      
    J4UnionSolid *solid1
	= new J4UnionSolid("solid1", sphere, polycone1, 0, centerOfPolycone);
      
    J4UnionSolid *solid2
	= new J4UnionSolid("solid2", solid1, tubs, 0, centerOfTubs);
      
    f8inchPMTSolid = new J4UnionSolid("solid", solid2, cons, 0, centerOfCons);
     
}

///////////////////////////////////////////////////////

 
