// $Id: J4UnionSolid.hh,v 1.1.1.1 2004/08/26 07:04:26 hoshina Exp $
#ifndef __J4UNIONSOLID__
#define __J4UNIONSOLID__
//*************************************************************************
//* --------------------
//* J4UnionSolid
//* --------------------
//* (Description)
//*     
//* (Update Record)
//*	2001/08/26  K.Hoshina	Original version.
//*************************************************************************

#include "G4UnionSolid.hh"

class J4UnionSolid : public G4UnionSolid
{
  public:  // with description

    J4UnionSolid( const G4String& pName,
                        G4VSolid* pSolidA ,
                        G4VSolid* pSolidB   ) ;

    J4UnionSolid( const G4String& pName,
                        G4VSolid* pSolidA ,
                        G4VSolid* pSolidB ,
                        G4RotationMatrix* rotMatrix,
                  const G4ThreeVector& transVector    ) ;

    J4UnionSolid( const G4String& pName,
                        G4VSolid* pSolidA ,
                        G4VSolid* pSolidB ,
                  const G4Transform3D& transform   ) ;

    virtual ~J4UnionSolid() ;

  public:  // without description

    inline virtual void  SetOwner(G4bool isowner = true)
    					{ fIsOwner = isowner; }
    
private:

    G4bool	fIsOwner;
    G4VSolid   *fSolidA;
    G4VSolid   *fSolidB; 
};

#endif



