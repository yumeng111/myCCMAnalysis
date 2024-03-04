// $Id: J4UnionSolid.cc,v 1.1.1.1 2004/08/26 07:04:26 hoshina Exp $
//*************************************************************************
//* --------------------
//* J4UnionSolid
//* --------------------
//* (Description)
//*     
//* (Update Record)
//*	2001/08/26  K.Hoshina	Original version.
//*************************************************************************

#include "g4-larsim/g4classes/J4UnionSolid.h"

///////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid

J4UnionSolid:: J4UnionSolid( const G4String& pName,
                                     G4VSolid* pSolidA ,
                                     G4VSolid* pSolidB   ):
		G4UnionSolid(pName,pSolidA,pSolidB), fIsOwner(true)
                , fSolidA(pSolidA), fSolidB(pSolidB)
{
   ;
}

/////////////////////////////////////////////////////////////////////
//
//
 
J4UnionSolid:: J4UnionSolid( const G4String& pName,
                                     G4VSolid* pSolidA ,
                                     G4VSolid* pSolidB ,
                                     G4RotationMatrix* rotMatrix,
                               const G4ThreeVector& transVector    ):
		G4UnionSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)
		, fIsOwner(true) , fSolidA(pSolidA), fSolidB(pSolidB)

{
   ;
}

///////////////////////////////////////////////////////////
//
//
 
J4UnionSolid:: J4UnionSolid( const G4String& pName,
                                     G4VSolid* pSolidA ,
                                     G4VSolid* pSolidB ,
                               const G4Transform3D& transform  ):
		G4UnionSolid(pName,pSolidA,pSolidB,transform)
		, fIsOwner(true) , fSolidA(pSolidA), fSolidB(pSolidB)

{
   ;
} 


///////////////////////////////////////////////////////////
//
//

J4UnionSolid::~J4UnionSolid()
{
#if 0
   if (fIsOwner) {
   	if (fSolidA) delete fSolidA;
   	if (fSolidB) delete fSolidB;
   }
#endif
}


