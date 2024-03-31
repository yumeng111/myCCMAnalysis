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
/// \file optical/LXe/include/G4CCMPMTHit.hh
/// \brief Definition of the G4CCMPMTHit class
//
//
#ifndef G4CCMPMTHit_h
#define G4CCMPMTHit_h 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"

class G4CCMPMTHit : public G4VHit {
    public:
        G4CCMPMTHit() = default;
        G4CCMPMTHit(const G4CCMPMTHit& right);
        ~G4CCMPMTHit() override = default;

        const G4CCMPMTHit& operator=(const G4CCMPMTHit& right);
        G4bool operator==(const G4CCMPMTHit& right) const;

        inline void* operator new(size_t);
        inline void operator delete(void* aHit);

        void Draw() override;
        void Print() override;

        inline void SetDrawit(G4bool b) { fDrawit = b; }
        inline G4bool GetDrawit() { return fDrawit; }

        inline void IncPhotonCount() { ++fPhotons; }
        inline G4int GetPhotonCount() { return fPhotons; }

        inline void SetPMTNumber(G4int n) { fPmtNumber = n; }
        inline G4int GetPMTNumber() { return fPmtNumber; }

        inline void SetPMTPhysVol(G4VPhysicalVolume* physVol) { fPhysVol = physVol; }
        inline G4VPhysicalVolume* GetPMTPhysVol() { return fPhysVol; }

        inline void SetPMTPos(G4double x, G4double y, G4double z) { fPos = G4ThreeVector(x, y, z); }

        inline G4ThreeVector GetPMTPos() { return fPos; }

        inline void SetCCMPMTKeyRow(G4int r) { fCCMPMTKeyRow = r; }        
        inline void SetCCMPMTKeyNumber(G4int n) { fCCMPMTKeyNumber = n; }

        inline G4int GetCCMPMTKeyRow() { return fCCMPMTKeyRow; }
        inline G4int GetCCMPMTKeyNumber() { return fCCMPMTKeyNumber; }

        inline void AppendPhotonTime(G4double t) { photonTime.push_back(t); }
        inline void AppendPhotonEnergy(G4double w) { photonEnergy.push_back(w); }
        inline void AppendPhotonPositionX(G4double x) { photonPositionX.push_back(x); }
        inline void AppendPhotonPositionY(G4double y) { photonPositionY.push_back(y); }
        inline void AppendPhotonPositionZ(G4double z) { photonPositionZ.push_back(z); }
        inline void AppendPhotonDirectionX(G4double x) { photonDirectionX.push_back(x); }
        inline void AppendPhotonDirectionY(G4double y) { photonDirectionY.push_back(y); }
        inline void AppendPhotonDirectionZ(G4double z) { photonDirectionZ.push_back(z); }
        inline void AppendPhotonCreationProcess(G4double c) { photonCreationProcess.push_back(c); }

        inline G4DataVector GetPhotonTime() { return photonTime;} 
        inline G4DataVector GetPhotonEnergy() { return photonEnergy; } 
        inline G4DataVector GetPhotonPositionX() { return photonPositionX; }
        inline G4DataVector GetPhotonPositionY() { return photonPositionY; }   
        inline G4DataVector GetPhotonPositionZ() { return photonPositionZ; }   
        inline G4DataVector GetPhotonDirectionX() { return photonDirectionX; }   
        inline G4DataVector GetPhotonDirectionY() { return photonDirectionX; }   
        inline G4DataVector GetPhotonDirectionZ() { return photonDirectionX; }   
        inline G4DataVector GetPhotonCreationProcess() { return photonCreationProcess; } 

    private:
        G4int fPmtNumber = -1;
        G4int fPhotons = 0;
        G4ThreeVector fPos;
        G4VPhysicalVolume* fPhysVol = nullptr;
        G4bool fDrawit = false;
        // we also need to keep track of things necessary to make a CCMPMTKey and information for CCMMCPE 
        // ideally, I would like to have a std::vector<CCMMCPE> but type needs to be G4...
        // there is a G4FastVector that is a template class so could have G4FastVector<CCMMCPE>
        // but the issue I ran into was you need to initialize it containing a specific data type and of a certain size
        // then you can add elements at a specific index (that's it, initialize, and add elements)
        // so since we don't know the size, so to .push_back equivalent, i would have to make a new G4FastVector of size N+1 than previous size,
        // add all previous elements, then add new element
        // so instead I am using G4DataVector (which inherit from std::vector<double>) to keep track of all elements necessary to make CCMMCPE
        // then the run manager will make the CCMMCPE objects and save appropriately
        G4int fCCMPMTKeyRow;
        G4int fCCMPMTKeyNumber;
        G4DataVector photonTime; 
        G4DataVector photonEnergy; 
        G4DataVector photonPositionX; 
        G4DataVector photonPositionY; 
        G4DataVector photonPositionZ; 
        G4DataVector photonDirectionX; 
        G4DataVector photonDirectionY; 
        G4DataVector photonDirectionZ; 
        G4DataVector photonCreationProcess; 
};

typedef G4THitsCollection<G4CCMPMTHit> G4CCMPMTHitsCollection;

extern G4ThreadLocal G4Allocator<G4CCMPMTHit>* G4CCMPMTHitAllocator;

inline void* G4CCMPMTHit::operator new(size_t) {
    if(!G4CCMPMTHitAllocator)
        G4CCMPMTHitAllocator = new G4Allocator<G4CCMPMTHit>;
    return (void*) G4CCMPMTHitAllocator->MallocSingle();
}

inline void G4CCMPMTHit::operator delete(void* aHit) {
    G4CCMPMTHitAllocator->FreeSingle((G4CCMPMTHit*) aHit);
}

#endif

