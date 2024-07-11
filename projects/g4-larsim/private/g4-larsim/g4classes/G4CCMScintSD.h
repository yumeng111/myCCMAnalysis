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
/// \file optical/LXe/include/LXeScintSD.hh
/// \brief Definition of the LXeScintSD class
//
//
#ifndef G4CCMScintSD_h
#define G4CCMScintSD_h 1

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"
#include "dataclasses/physics/I3MCTreeUtils.h"

#include <vector>
#include <string>
#include <sstream>

#include <G4SystemOfUnits.hh>
#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;

class G4CCMScintSD : public G4VSensitiveDetector {
    public:
        G4CCMScintSD(G4String name);
        ~G4CCMScintSD() override = default;

        void Initialize(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;
        
        // return updated MCTree
        I3MCTreePtr GetUpdatedMCTree(){ return mcTree; }
        I3MCTreePtr GetPhotonSummary(){ return photon_summary; }
        
        bool GetPMTSDStatus() { return PMTSDStatus_; }
        void SetPMTSDStatus(bool PMTSDStatus) { PMTSDStatus_ = PMTSDStatus; }

        void SetPrimaryParticle(I3Particle primary) {
            // now let's set the primary particle
            ClearResults();

            primaryParticleType_ = primary.GetType();
            DaughterParticleMap[1] = primary.GetID();
            I3MCTreeUtils::AddPrimary(*mcTree, primary);
        } 
        
        void ClearResults() {
            mcTree->clear();
            photon_summary->clear();
            DaughterParticleMap.clear();
            DaughterOpticalPhotonMap.clear();
        } 

    private:
        G4CCMScintHitsCollection* fScintCollection = nullptr;
        G4int fHitsCID = -1;
        
        // controls to turn SD on/off (set by our response service)
        // we just need to know if PMTSD is on --> if so, we do NOT kill tracks, but if it is off, we DO kill tracks after registering one hit
        bool PMTSDStatus_; 

        // let's also grab the primary particle information
        I3Particle::ParticleType primaryParticleType_;

        // our mc tree we will save energy depositions to
        I3MCTreePtr mcTree = boost::make_shared<I3MCTree>();

        // and photon summary that keeps track of optical photons
        I3MCTreePtr photon_summary = boost::make_shared<I3MCTree>();
        
        static const std::unordered_map<std::string, int> energyLossToI3ParticlePDGCode;

        std::map<int, I3ParticleID> DaughterParticleMap; // map between track_id and I3ParticleID
        std::map<int, I3ParticleID> DaughterOpticalPhotonMap; // map between track_id and I3ParticleID
        bool firstOpPh = true;
        int firstOpPhParentID;
};

#endif
