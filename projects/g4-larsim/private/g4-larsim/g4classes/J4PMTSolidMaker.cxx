// $Id: J4PMTSolidMaker.cc,v 1.2 2007/03/13 17:29:14 hoshina Exp $
/**
* @file J4PMTSolidMaker.cc
* @brief namespace that provides geometory functions to generate G4Solids
* @date 2021/04/11
* @note many numbers are hard coded, if you want to change it on the fly move these parameters to J4PartsParameterList.
* (Update Record)
*	2021/04/11  K.Hoshina	collected functions from each detector classes
*/

#include "g4-larsim/g4classes/J4UnionSolid.h"
#include "g4-larsim/g4classes/J4PMTSolidMaker.h"

#include <G4Cons.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Polycone.hh>
#include <G4UserLimits.hh>
#include <G4SystemOfUnits.hh>
#include <G4OpticalSurface.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalSkinSurface.hh>

G4VSolid * J4PMTSolidMaker::f8inchPMTSolid;
G4VSolid * J4PMTSolidMaker::fTPBCoatingSolid;
G4VSolid * J4PMTSolidMaker::fPhotocathodeSolid;
G4double   J4PMTSolidMaker::f8inchPMTRadius = 131.0 * mm;
G4double   J4PMTSolidMaker::fTPBThickness = 0.002 * mm;
G4double   J4PMTSolidMaker::fPhotocathodeDelta = 1.0 * mm; // idk just needs to be greater than the pmt glass thickness

//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* J4PMTSolidMaker::Get8inchPMTSolid()
{
   if (!f8inchPMTSolid) {
       Create8inchPMTSolid();
   }

   return f8inchPMTSolid;
}

G4VSolid* J4PMTSolidMaker::GetTPBCoatingSolid()
{
   if (!fTPBCoatingSolid) {
       CreateTPBCoatingSolid();
   }

   return fTPBCoatingSolid;
}

G4VSolid* J4PMTSolidMaker::GetPhotcathodeSolid()
{
   if (!fPhotocathodeSolid) {
       CreatePhotocathodeSolid();
   }

   return fPhotocathodeSolid;
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
    }
    
    for (G4int i=0; i<=segment; i++) {
      rInner[i] = tempInner[i];
      rOuter[i] = tempOuter[i];
      zPlane[i] = tempZ[i];
    }
    
    G4Polycone* polycone1 = new 
      G4Polycone("polycone1", 0, 2*M_PI, segment+1, zPlane, rInner, rOuter);
    
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
      
    J4UnionSolid *solid1 = new J4UnionSolid("solid1", sphere, polycone1, 0, centerOfPolycone);
      
    J4UnionSolid *solid2 = new J4UnionSolid("solid2", solid1, tubs, 0, centerOfTubs);
      
    f8inchPMTSolid = new J4UnionSolid("solid", solid2, cons, 0, centerOfCons);

}

void J4PMTSolidMaker::CreateTPBCoatingSolid()
{

    // so we need to make tpb coating that is an even thickness
    // to do this, we modify the paramters as such:
    // sr0 --> sr0 + delta_tpb  
    // rSmax --> rSmax + delta_tpb
    // rmax --> rmax + delta_tpb*1.75
    // rmax2 --> rmax2 + delta_tpb*1.75
    // centerofConslb --> centerofConslb - delta_tpb/2.1 (centerofConslb is 70.4284)
    // then after making slightly bigger pmt, we subtract 8in pmt from big pmt and that's our tpb coating!

    G4ThreeVector centerOfPolycone;
    G4ThreeVector centerOfCons; 
    G4ThreeVector centerOfTubs; 
  
    G4double rmin2, rmax2, dz;
    G4double rmin, rmax, sphi, dphi, stheta, dtheta;
    G4double rSmin, rSmax, sSphi, dSphi, sStheta, dStheta;
    G4double centerofConslb;
    centerofConslb = 70.4284 - fTPBThickness / 2.1;

    rSmin   = 0;
    rSmax   = f8inchPMTRadius + fTPBThickness;
    dSphi   = 2*M_PI;
    dStheta = 37.6392 *degree;
     
    std::vector<G4double> tempZ, tempInner, tempOuter;

    G4int segment = 100;
  
    G4double sr0 = 57.0894 + fTPBThickness;//radius of spindle torus sphere
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
    rmax   = 84.5/2*mm + (1.75*fTPBThickness)/2;
    rmax2  = 160./2 *mm + (1.75*fTPBThickness)/2;  
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
      
    J4UnionSolid *solid1 = new J4UnionSolid("solid1", sphere, polycone1, 0, centerOfPolycone);
      
    J4UnionSolid *solid2 = new J4UnionSolid("solid2", solid1, tubs, 0, centerOfTubs);
      
    G4VSolid * biggerPMT = new J4UnionSolid("solid", solid2, cons, 0, centerOfCons);

    // now grab an 8in pmt
    if (!f8inchPMTSolid) {
        Create8inchPMTSolid();
    }

    // now subtract!
    fTPBCoatingSolid = new G4SubtractionSolid("TPBShell", f8inchPMTSolid, biggerPMT);
}

///////////////////////////////////////////////////////

void J4PMTSolidMaker::CreatePhotocathodeSolid()
{

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

    rSmin   = 0;
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
      
    J4UnionSolid *solid1 = new J4UnionSolid("solid1", sphere, polycone1, 0, centerOfPolycone);
      
    J4UnionSolid *solid2 = new J4UnionSolid("solid2", solid1, tubs, 0, centerOfTubs);
      
    fPhotocathodeSolid = new J4UnionSolid("solid", solid2, cons, 0, centerOfCons);
}
 
