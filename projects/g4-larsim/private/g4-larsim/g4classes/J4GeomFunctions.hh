// $Id: J4GeomFunctions.hh,v 1.1.1.1 2004/08/26 07:04:26 hoshina Exp $
#ifndef __J4GeomFunctions__hh
#define __J4GeomFunctions__hh
//*************************************************************************
//* --------------------
//* J4GeomFunctions
//* --------------------
//* (Description)
//* 	J4GeomFunctions is a set of geometory construction function
//*    
//* (Update Record)
//*	2021/5/20  K.Hoshina	Original version.
//*************************************************************************

#include "G4Types.hh"
#include "G4VSolid.hh"


namespace J4GeomFunctions {	

   void FillTorusSurface(G4int nsegments, // number
                      G4double step,
                      G4double torus_r,
                      G4double torus_zmax,
                      G4double centerOfTorus_z,
                      G4double centerOfTorus_r,
                      std::vector<G4double> & tempZ,
                      std::vector<G4double> & tempOuter);

   G4VSolid* CreateEggTorus2
          (G4int    nsegments2,       // nsegments of torus
           G4double torus2_r,         // radius of torus2
           G4double centerOfTorus2_r, // origin of torus2 r
           G4double centerOfTorus2_z, // origin of torus2 z
           G4double torus2_zmin,      // zmin of torus2
           G4double torus2_zmax,      // zmax of torus2
           G4double torus2_z0        // r of torus2 at z=0
           ) ; // origin of torus1 in mother's coordinate



   G4VSolid* CreateHalfEggSolid
          (G4int    nsegments1,       // nsegments of sphere 
           G4double sphere_rmax,      // radius of sphere
           G4double sphere_dtheta,    // solid angle in theta
           G4double sphere_transform_z,  // origin of sphere z in mother's coordinate
           G4double torus1_r,         // radius of torus1
           G4double centerOfTorus1_r, // origin of torus1 r
           G4int    nsegments2,       // nsegments of torus
           G4double torus2_r,         // radius of torus2
           G4double centerOfTorus2_r, // origin of torus2 r
           G4double centerOfTorus2_z, // origin of torus2 z
           G4double torus2_zmin,      // zmin of torus2
           G4double torus2_zmax,      // zmax of torus2
           G4double torus2_z0,        // r of torus2 at z=0
           G4double torus1_transform_z) ; // origin of torus1 in mother's coordinate

   G4VSolid* CreateEggSolid
          (G4int    nsegments1,       // nsegments of sphere 
           G4double sphere_rmax,      // radius of sphere
           G4double sphere_dtheta,    // solid angle in theta
           G4double sphere_transform_z,  // origin of sphere z in mother's coordinate
           G4double torus1_r,         // radius of torus1
           G4double centerOfTorus1_r, // origin of torus1 r
           G4int    nsegments2,       // nsegments of torus
           G4double torus2_r,         // radius of torus2
           G4double centerOfTorus2_r, // origin of torus2 r
           G4double centerOfTorus2_z, // origin of torus2 z
           G4double torus2_zmin,      // zmin of torus2
           G4double torus2_zmax,      // zmax of torus2
           G4double torus2_z0,        // r of torus2 at z=0
           G4double torus1_transform_z) ; // origin of torus1 in mother's coordinate

}

#endif

