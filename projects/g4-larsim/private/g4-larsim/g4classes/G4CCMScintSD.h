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
#include "g4-larsim/g4classes/G4CCMReadout.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "dataclasses/I3Map.h"
#include "simclasses/PhotonSummary.h"

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

        void SetReadout(G4CCMReadout * readout) { readout_ = readout; }

        void Initialize(G4HCofThisEvent*) override;
        void EndOfEvent(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;

        // return updated MCTree
        PhotonSummarySeriesPtr GetPhotonSummarySeries(){ return photon_summary; }
        boost::shared_ptr<I3Map<int, size_t>> GetPhotonSummaryMap() { return optical_photon_map; }

        bool GetPMTSDStatus() { return PMTSDStatus_; }
        void SetPMTSDStatus(bool PMTSDStatus) { PMTSDStatus_ = PMTSDStatus; }

        bool GetTimeCutStatus() { return TimeCut_; }
        void SetTimeCutStatus(bool TimeCut) { TimeCut_ = TimeCut; }

        bool GetKillCherenkovStatus() { return KillCherenkov_; }
        void SetKillCherenkovStatus(bool KillCherenkov) { KillCherenkov_ = KillCherenkov; }

        void Reset() {
            DaughterParticleMap.clear();
            photon_summary = boost::make_shared<PhotonSummarySeries>();
            optical_photon_map = boost::make_shared<I3Map<int, size_t>>();
            mcTree = nullptr;
        }

        void AddEntryToPhotonSummary(int parent_id, int track_id, double g4_uv_distance, double g4_vis_distance,
                                     double calculated_uv_distance, double calculated_vis_distance,
                                     double g4_time, double calculated_time, std::string creationProcessName);
        void UpdatePhotonSummary(int parent_id, int track_id, double g4_uv_distance, double g4_vis_distance,
                                 double calculated_uv_distance, double calculated_vis_distance,
                                 double g4_time, double calculated_time, std::string creationProcessName,
                                 std::map<int, size_t>::iterator it, bool new_process);

        double InterpolateRindex(double wavelength);

    private:
        int event_id = -1;
        G4CCMReadout * readout_ = nullptr;
        G4CCMScintHitsCollection* fScintCollection = nullptr;
        G4int fHitsCID = -1;

        // controls to turn SD on/off (set by our response service)
        // we just need to know if PMTSD is on --> if so, we do NOT kill tracks, but if it is off, we DO kill tracks after registering one hit
        bool PMTSDStatus_;

        bool TimeCut_; // true kills tracks after 200 nsec
        bool KillCherenkov_; // true kills cerenkov light

        // our mc tree we will save energy depositions to
        I3MCTreePtr mcTree = nullptr;
        I3Particle primary_;

        // and photon summary that keeps track of optical photons
        // note this is in two parts -- map between track id and idx
        // and vector containing photon summaries (idx in vector corresponds to track id)
        boost::shared_ptr<I3Map<int, size_t>> optical_photon_map = boost::make_shared<I3Map<int, size_t>>(); // map between track id and idx in vector
        PhotonSummarySeriesPtr photon_summary = boost::make_shared<PhotonSummarySeries>();

        // define a few things for converting energy to wavelength
        const G4double hbarc = 197.326; //eV * nm
        const G4double hc =  hbarc * (2 * 3.14159265358979323846); // eV * nm

        static const std::unordered_map<std::string, int> energyLossToI3ParticlePDGCode;
        static const std::unordered_map<std::string, PhotonSummary::PhotonSource> processNameToPhotonSource;
        static const std::unordered_map<PhotonSummary::PhotonSource, std::string> photonSourceToProcessName;

        std::map<int, I3ParticleID> DaughterParticleMap; // map between track_id and I3ParticleID

        double c_mm_per_nsec = 299.792458 * (I3Units::mm / I3Units::ns); // speed of light in mm/nsec

        std::vector<double> rindex_wavelength = {99.99909859236817, 109.99897296368657, 119.99889895379226, 123.99888385471958,
                                                 127.99884206872682, 133.99878562916103, 139.99873802931543, 159.99855258590853,
                                                 179.99837456270487, 199.99819718473634, 299.9973199734223, 399.99632984597423,
                                                 499.9954929618409, 599.9944947689613, 699.9936901465773}; // nm

        std::vector<double> rindex = {1.7898800000000006, 1.6199600000000003, 1.4500400000000002, 1.40284, 1.358,
                                      1.335167, 1.315166, 1.28071, 1.263128, 1.25475, 1.24, 1.23, 1.225, 1.222, 1.22};
};

#endif
