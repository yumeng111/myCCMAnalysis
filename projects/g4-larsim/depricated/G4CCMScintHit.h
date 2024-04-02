//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file optical/LXe/include/LXeScintHit.hh
/// \brief Definition of the LXeScintHit class
//
//
#ifndef G4CCMScintHit_h
#define G4CCMScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

class G4CCMScintHit : public G4VHit
{
 public:
  G4CCMScintHit() = default;
  G4CCMScintHit(G4VPhysicalVolume* pVol);
  ~G4CCMScintHit() override = default;

  G4CCMScintHit(const G4CCMScintHit& right);
  const G4CCMScintHit& operator=(const G4CCMScintHit& right);
  G4bool operator==(const G4CCMScintHit& right) const;

  inline void* operator new(size_t);
  inline void operator delete(void* aHit);

  inline void SetEdep(G4double de) { fEdep = de; }
  inline void AddEdep(G4double de) { fEdep += de; }
  inline G4double GetEdep() { return fEdep; }

  inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }
  inline G4ThreeVector GetPos() { return fPos; }

  inline const G4VPhysicalVolume* GetPhysV() { return fPhysVol; }

 private:
  G4double fEdep = 0.;
  G4ThreeVector fPos;
  const G4VPhysicalVolume* fPhysVol = nullptr;
};

typedef G4THitsCollection<G4CCMScintHit> G4CCMScintHitsCollection;

extern G4ThreadLocal G4Allocator<G4CCMScintHit>* G4CCMScintHitAllocator;

inline void* G4CCMScintHit::operator new(size_t)
{
  if(!G4CCMScintHitAllocator)
    G4CCMScintHitAllocator = new G4Allocator<G4CCMScintHit>;
  return (void*) G4CCMScintHitAllocator->MallocSingle();
}

inline void G4CCMScintHit::operator delete(void* aHit)
{
  G4CCMScintHitAllocator->FreeSingle((G4CCMScintHit*) aHit);
}

#endif
