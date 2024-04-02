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

#include "g4-larsim/g4classes/G4CCMScintHit.h"

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"

#include <vector>
#include <sstream>
#include <string>
#include <simclasses/CCMMCPE.h>
#include <icetray/CCMPMTKey.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>

class G4Step;
class G4HCofThisEvent;

class G4CCMScintSD : public G4VSensitiveDetector {
    public:
        G4CCMScintSD(G4String name);
        ~G4CCMScintSD() override = default;

        void Initialize(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;
        
        // return list of CCMMCPE
        boost::shared_ptr<CCMMCPESeries> GetCCMMCPEList(){ return CCMMCPEList; }
        
        bool GetPMTSDStatus() { return PMTSDStatus_; }
        void SetPMTSDStatus(bool PMTSDStatus) { PMTSDStatus_ = PMTSDStatus; }

    private:
        G4CCMScintHitsCollection* fScintCollection = nullptr;
        G4int fHitsCID = -1;
        static const std::unordered_map<std::string, CCMMCPE::PhotonSource> processNameToPhotonSource;
        // define a few things for converting energy to wavelength
        const G4double h_Planck = 6.62607015e-34 * joule * second;
        const G4double c_light = 2.99792458e8 * meter / second;
        
        // controls to turn SD on/off (set by our response service)
        // we just need to know if PMTSD is on --> if so, we do NOT kill tracks, but if it is off, we DO kill tracks after registering one hit
        bool PMTSDStatus_; 
        
        boost::shared_ptr<CCMMCPESeries> CCMMCPEList = boost::make_shared<CCMMCPESeries> ();
};

#endif
