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
/// \file optical/G4CCM/include/G4CCMPMTSD.hh
/// \brief Definition of the G4CCMPMTSD class
//
//
#ifndef G4CCMPMTSD_h
#define G4CCMPMTSD_h 1

#include "icetray/I3Units.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Direction.h"

#include "g4-larsim/g4classes/G4CCMPMTHit.h"
#include "g4-larsim/g4classes/G4CCMMainVolume.h"
#include "g4-larsim/g4classes/G4CCMReadout.h"

#include "icetray/CCMPMTKey.h"

#include "simclasses/CCMMCPE.h"

#include <vector>
#include <string>
#include <sstream>

#include <G4SystemOfUnits.hh>
#include <G4VSensitiveDetector.hh>


class G4DataVector;
class G4HCofThisEvent;
class G4Step;
class G4CCMMainVolume;

class G4CCMPMTSD : public G4VSensitiveDetector {
    public:
        G4CCMPMTSD(G4String name);
        ~G4CCMPMTSD() override;

        void SetReadout(G4CCMReadout * readout) { readout_ = readout; }

        void Initialize(G4HCofThisEvent*) override;
        void EndOfEvent(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;

        // Initialize the arrays to store pmt positions
        inline void InitPMTs() {
            if(fPMTPositionsX)
                delete fPMTPositionsX;
            if(fPMTPositionsY)
                delete fPMTPositionsY;
            if(fPMTPositionsZ)
                delete fPMTPositionsZ;
            fPMTPositionsX = new G4DataVector();
            fPMTPositionsY = new G4DataVector();
            fPMTPositionsZ = new G4DataVector();
        }

        // Store a pmt position
        void SetPmtPositions(const std::vector<G4ThreeVector>& positions);

        // map from volume names to CCMPMTKeys
        std::map<std::string, CCMPMTKey> volumeToKey;

        // return CCMMCPEMap
        boost::shared_ptr<CCMMCPESeriesMap> GetCCMMCPEMap(){ return CCMMCPEMap; }

        void Reset() {
            CCMMCPEMap = nullptr;
        }
        
        void SetPhotonTracking(bool FullPhotonTracking) { FullPhotonTracking_ = FullPhotonTracking; }

    private:
        bool FullPhotonTracking_;
        int event_id = -1;
        G4CCMReadout * readout_ = nullptr;
        G4CCMPMTHitsCollection* fPMTHitCollection = nullptr;

        G4DataVector* fPMTPositionsX = nullptr;
        G4DataVector* fPMTPositionsY = nullptr;
        G4DataVector* fPMTPositionsZ = nullptr;

        // define a few things for converting energy to wavelength
        const G4double hbarc = 197.326; //eV * nm
        const G4double hc =  hbarc * (2 * 3.14159265358979323846); // eV * nm

        G4int fHitCID = -1;

        static const std::unordered_map<std::string, CCMMCPE::PhotonSource> processNameToPhotonSource;

        boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap = nullptr;
};

#endif

