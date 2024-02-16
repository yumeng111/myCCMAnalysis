// $Id: J4GeomFunctions.cc,v 1.2 2007/03/13 17:29:14 hoshina Exp $
/**
* @file J4GeomFunction.cc
* @brief namespace that provides geometory functions to generate G4Solids
* @date 2021/04/11
* @note many numbers are hard coded, if you want to change it on the fly move these parameters to J4EggParameterList.
* (Update Record)
*	2021/04/11  K.Hoshina	collected functions from each detector classes
*/

#include "J4GeomFunctions.hh"
#include "G4Sphere.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "J4UnionSolid.hh"
#include "G4UserLimits.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4SystemOfUnits.hh"
#include <algorithm>

namespace J4GeomFunctions {

//////////////////////////////////////////////////////////////////////////////////////
/// Calculate z and radius values for G4Polycone
/// 
/// nsegments : number of segments in z axis
/// step : size of step along with z axis
/// torus_r : radius of torus surface
/// torus_zmax : maximum z value to calculate torus surface
///              (minimum z will be torus_zmax - nsegents * step)
/// centerOfTorus_z : origin of the torus in z-axis, signed
/// centerOfTorus_r : origin of the torus in radius, signed
//////////////////////////////////////////////////////////////////////////////////////
void FillTorusSurface(G4int nsegments, // number 
                      G4double step, 
                      G4double torus_r,
                      G4double torus_zmax, 
                      G4double centerOfTorus_z, 
                      G4double centerOfTorus_r,
                      std::vector<G4double> & tempZ,
                      std::vector<G4double> & tempOuter)
{

/*
    std::cout << " nsegments " << nsegments << " step " << step
           << " torus_r " << torus_r << " torus_zmax " << torus_zmax
           << " centerOfTorus_z " << centerOfTorus_z 
           << " centerOfTorus_r " << centerOfTorus_r << G4endl;
*/

    G4double r;
    G4double torus_relative_zmax = torus_zmax - centerOfTorus_z;
    for (G4int j=0; j<=nsegments; ++j) {
       tempZ.push_back(torus_zmax - j*step);
       // r = sqrt(torus_r*torus_r - (torus_relative_zmax-j*step)*(torus_relative_zmax-j*step))
       r = sqrt((torus_r + torus_relative_zmax - j*step)*(torus_r - torus_relative_zmax + j*step));
       tempOuter.push_back(centerOfTorus_r + r);
    }
}


//////////////////////////////////////////////////////////////////////////////////////
/// Creating sphere part
//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* CreateEggSphere
          (G4double sphere_rmax,      // radius of sphere
           G4double sphere_dtheta    // solid angle in theta
          )   // origin of torus1 in mother's coordinate, signed
{
    G4double rmin   = 0;
    G4double rmax   = sphere_rmax;
    G4double sphi   = 0. *degree;
    G4double dphi   = 2*M_PI;
    G4double dtheta = sphere_dtheta;
    G4double stheta = 0. *degree;
 
    G4Sphere * sphere = new G4Sphere("sphere",rmin, rmax, 
				     sphi, dphi, stheta, dtheta);
    return sphere;   
}

//////////////////////////////////////////////////////////////////////////////////////
/// Creating torus 1 part
//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* CreateEggTorus1
          (G4int    nsegments1,       // nsegments for torus1
           G4double torus1_r,         // radius of torus1
           G4double centerOfTorus1_r // origin of torus1 r
           )   // origin of torus1 in mother's coordinate, signed
{
    // building small polycones. 
    // Here we define full size of torus, but eventually only part of it will be used after
    // G4UnionSolid is defined with another polycone and sphere.

    G4double step = torus1_r / nsegments1; // divide diameter of the torus into segment1 numbers.
    G4double torus1_zmax = torus1_r;

    std::vector<G4double> tempZ, tempOuter;
    J4GeomFunctions::FillTorusSurface(nsegments1,
                        step, torus1_r, torus1_zmax,
                        0, centerOfTorus1_r,
                        tempZ, tempOuter);

    G4double rInner[nsegments1+1], rOuter[nsegments1+1], zPlane[nsegments1+1];
    for (G4int i=0; i<=nsegments1; ++i) {
       zPlane[i] = tempZ[i];
       rInner[i] = 0.;
       rOuter[i] = tempOuter[i];
       //std::cout <<"EggTorus1 "<<i<<", "<<zPlane[i] <<", "<< rOuter[i] << G4endl;
    }
    
    G4Polycone * torus1 = new 
      G4Polycone("torus1", 0, 2*M_PI, nsegments1+1, zPlane, rInner, rOuter);
    return torus1;
}

//////////////////////////////////////////////////////////////////////////////////////
/// Creating torus 2 part
//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* CreateEggTorus2
          (
           G4int    nsegments2,       // nsegments for torus2
           G4double torus2_r,         // radius of torus2
           G4double centerOfTorus2_r, // origin of torus2 r
           G4double centerOfTorus2_z, // origin of torus2 z, signed
           G4double torus2_zmin,      // zmin of torus2
           G4double torus2_zmax,      // zmax of torus2
           G4double torus2_z0         // r of torus2 at z=0
           )   
{
    //
    //building the large sphere revolution
    // this does not work when drawing with segment1=40
    G4double zminrelative = torus2_zmin - centerOfTorus2_z; //minimum z shift from center of torus in positive z direction
    G4double zmaxrelative = torus2_zmax - centerOfTorus2_z; //maximum z shift from center of torus in positive z direction

    std::vector<G4double> tempZ2, tempOuter2;
    G4double step2 = (zmaxrelative-zminrelative)/(nsegments2-1);

/*
    std::cout << "nsegments2 " << nsegments2 << " step2 " << step2 << G4endl;
    std::cout << "zmaxrelative " << zmaxrelative << " zminrelative " << zminrelative << G4endl;
    std::cout << "torus2_zmin " << torus2_zmin << " centerObTorus2_z " << centerOfTorus2_z << G4endl;
    std::cout << "torus2_zmax " << torus2_zmax << " centerObTorus2_z " << centerOfTorus2_z << G4endl;
*/

    // upper half
    J4GeomFunctions::FillTorusSurface(nsegments2-1,
                        step2, torus2_r, torus2_zmax,
                        centerOfTorus2_z, centerOfTorus2_r,
                        tempZ2, tempOuter2);

    // add point on z=0
    tempZ2.push_back(0.);
    tempOuter2.push_back(torus2_z0);

    G4double rInner2[nsegments2+1], rOuter2[nsegments2+1], zPlane2[nsegments2+1];
    for (G4int i=0; i<=nsegments2; i++) {
       rInner2[i] = 0;
       rOuter2[i] = tempOuter2[i];
       zPlane2[i] = tempZ2[i];
       std::cout<<"EggTorus2 "<<i<<" "<<zPlane2[i]<<" "<<rInner2[i]<<" "<<rOuter2[i]<<std::endl;
    }

    G4Polycone * torus2 = new 
      G4Polycone("polycone2", 0, 2*M_PI, nsegments2+1, zPlane2, rInner2, rOuter2);

    return torus2;
}

//////////////////////////////////////////////////////////////////////////////////////
/// Creating half Egg glass shape(solid). 
/// Egg glass solid consists of :
///   two large polycones (outer shape is defined with torus function)
///   two small polycones (outer shape is defined with torus function)
///   two spheres
//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* CreateHalfEggSolid
          (G4int    nsegments1,       // nsegments for torus1
           G4double sphere_rmax,      // radius of sphere
           G4double sphere_dtheta,    // solid angle in theta
           G4double sphere_transform_z,  // origin of sphere z in mother's coordinate, signed
           G4double torus1_r,         // radius of torus1
           G4double centerOfTorus1_r, // origin of torus1 r
           G4int    nsegments2,       // nsegments for torus2
           G4double torus2_r,         // radius of torus2
           G4double centerOfTorus2_r, // origin of torus2 r
           G4double centerOfTorus2_z, // origin of torus2 z, signed
           G4double torus2_zmin,      // zmin of torus2
           G4double torus2_zmax,      // zmax of torus2
           G4double torus2_z0,        // r of torus2 at z=0
           G4double torus1_transform_z)   // origin of torus1 in mother's coordinate, signed
{
    //creating two spheres for top part
    G4VSolid *sphere = CreateEggSphere(sphere_rmax, sphere_dtheta);
    G4ThreeVector centerOfSphereUp(0, 0, sphere_transform_z); 

    // building small polycones. 
    // Here we define full size of torus, but eventually only part of it will be used after
    // G4UnionSolid is defined with another polycone and sphere.
    G4VSolid * torus1 = CreateEggTorus1(nsegments1, torus1_r, centerOfTorus1_r); 

    //
    //building the large sphere revolution
    // this does not work when drawing with segment1=40

    G4VSolid * torus2 = CreateEggTorus2(nsegments2, torus2_r, 
                       centerOfTorus2_r, centerOfTorus2_z,
                       torus2_zmin, torus2_zmax, torus2_z0);

    G4ThreeVector centerOfPolycone1(0, 0, torus1_transform_z); 

    J4UnionSolid *solid1
      = new J4UnionSolid("solid1", torus2, torus1, 0, centerOfPolycone1);
    
    J4UnionSolid *solid
        = new J4UnionSolid("solid", solid1, sphere, 0, centerOfSphereUp);
    
    return solid;

}

//////////////////////////////////////////////////////////////////////////////////////
/// Creating full Egg glass shape(solid). 
/// Egg glass solid consists of :
///   two large polycones (outer shape is defined with torus function)
///   two small polycones (outer shape is defined with torus function)
///   two spheres
//////////////////////////////////////////////////////////////////////////////////////
G4VSolid* CreateEggSolid
          (G4int    nsegments1,       // nsegments for torus1
           G4double sphere_rmax,      // radius of sphere
           G4double sphere_dtheta,    // solid angle in theta
           G4double sphere_transform_z,  // origin of sphere z in mother's coordinate, signed
           G4double torus1_r,         // radius of torus1
           G4double centerOfTorus1_r, // origin of torus1 r
           G4int    nsegments2,       // nsegments for torus2
           G4double torus2_r,         // radius of torus2
           G4double centerOfTorus2_r, // origin of torus2 r
           G4double centerOfTorus2_z, // origin of torus2 z, signed
           G4double torus2_zmin,      // zmin of torus2
           G4double torus2_zmax,      // zmax of torus2
           G4double torus2_z0,        // r of torus2 at z=0
           G4double torus1_transform_z)   // origin of torus1 in mother's coordinate, signed
{

   G4VSolid * deggup = CreateHalfEggSolid (
              nsegments1,
              sphere_rmax,
              sphere_dtheta,
              sphere_transform_z,
              torus1_r,
              centerOfTorus1_r,
              nsegments2,
              torus2_r, 
              centerOfTorus2_r, 
              centerOfTorus2_z, 
              torus2_zmin, 
              torus2_zmax,
              torus2_z0,  
              torus1_transform_z);

   G4VSolid * deggdown = CreateHalfEggSolid (
              nsegments1,
              sphere_rmax,
              sphere_dtheta,
              sphere_transform_z,
              torus1_r,
              centerOfTorus1_r,
              nsegments2,
              torus2_r, 
              centerOfTorus2_r, 
              centerOfTorus2_z, 
              torus2_zmin, 
              torus2_zmax,
              torus2_z0,  
              torus1_transform_z);


   G4RotationMatrix *rot = new G4RotationMatrix();
   rot->rotateY(180.0*deg);
   J4UnionSolid * degg 
      = new J4UnionSolid("degg", deggup, deggdown, rot, G4ThreeVector());

   return degg;

}
}
