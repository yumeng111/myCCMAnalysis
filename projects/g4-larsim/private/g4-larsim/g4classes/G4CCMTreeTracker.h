#ifndef G4CCMTreeTracker_H
#define G4CCMTreeTracker_H

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
#include <unordered_map>
#include <unordered_set>

#include <G4SystemOfUnits.hh>
#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;
class G4CCMEDepSD;

class G4CCMTreeTracker : public G4VSensitiveDetector {
    public:
        friend G4CCMEDepSD;
        G4CCMTreeTracker(G4String name);
        ~G4CCMTreeTracker() override = default;

        void SetReadout(G4CCMReadout * readout) {readout_ = readout;}

        void Initialize(G4HCofThisEvent*) override;
        void EndOfEvent(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;

        void SetTracking(bool do_tracking);

        bool GetKillNeutrinos() { return KillNeutrinos_; }
        void SetKillNeutrinos(bool KillNeutrinos) { KillNeutrinos_ = KillNeutrinos; }

        bool GetKillPhotons() { return KillPhotons_; }
        void SetKillPhotons(bool KillPhotons) { KillPhotons_ = KillPhotons; }

        double GetTimeCutValue() { return time_cut_value_; }
        void SetTimeCutValue(double time_cut) { time_cut_value_ = time_cut; }

        bool GetTimeCut() { return TimeCut_; }
        void SetTimeCut(bool TimeCut) { TimeCut_ = TimeCut; }

        bool GetKillCherenkov() { return KillCherenkov_; }
        void SetKillCherenkov(bool KillCherenkov) { KillCherenkov_ = KillCherenkov; }

        bool GetKillScintillation() { return KillScintillation_; }
        void SetKillScintillation(bool KillScintillation) { KillScintillation_ = KillScintillation; }

        bool GetDetailedPhotonTracking() { return DetailedPhotonTracking_; }
        void SetDetailedPhotonTracking(bool DetailedPhotonTracking) { DetailedPhotonTracking_ = DetailedPhotonTracking; }

        bool GetTrackParticles() { return TrackParticles_; }
        void SetTrackParticles(bool TrackParticles) { TrackParticles_ = TrackParticles; }

        bool GetTrackEnergyLosses() { return TrackEnergyLosses_; }
        void SetTrackEnergyLosses(bool TrackEnergyLosses) { TrackEnergyLosses_ = TrackEnergyLosses; }

        G4double GetG4RangeCut() { return G4RangeCut_ / I3Units::MeV * MeV; }
        void SetG4RangeCut(G4double G4RangeCut) { G4RangeCut_ = G4RangeCut / MeV * I3Units::MeV; }
        G4double GetG4EDepMin() { return G4EDepMin_ / I3Units::MeV * MeV; }
        void SetG4EDepMin(G4double G4EDepMin) { G4EDepMin_ = G4EDepMin / MeV * I3Units::MeV; }
        G4double GetG4ETrackingMin() { return G4ETrackingMin_ / I3Units::MeV * MeV; }
        void SetG4ETrackingMin(G4double G4ETrackingMin) { G4ETrackingMin_ = G4ETrackingMin / MeV * I3Units::MeV; }


        void Reset() {
            DaughterParticleMap.clear();
            photon_summary = boost::make_shared<I3Map<int, PhotonSummary>>();
            wls_parent_daughter_map = boost::make_shared<I3Map<int, std::set<int>>>();
            sub_threshold_losses.clear();
            parent_map.clear();
            mcTree = nullptr;
        }

        void AddNewPhoton(int parent_id, int track_id, double time, double distance, double wavelength, std::string creation_process_name);
        void AddPhotonTrack(int parent_id, int track_id, double delta_time, double delta_distance, std::string creation_process_name);
        void UpdatePhoton(int parent_id, int track_id, double delta_time, double delta_distance, std::string creation_process_name);
        void AddWLSPhotonTrack(int parent_id, int track_id, double delta_time, double delta_distance, G4Step* aStep);


    private:
        int event_id = -1;
        G4CCMReadout * readout_ = nullptr;

        bool KillNeutrinos_ = false;
        bool KillPhotons_ = false;
        bool KillScintillation_ = false;
        bool KillCherenkov_ = false;
        bool TimeCut_ = false;
        bool DetailedPhotonTracking_ = false;
        bool TrackParticles_ = false;
        bool TrackEnergyLosses_ = false;

        double time_cut_value_ = 200.0;

        double G4RangeCut_ = 0.0; // range cut for all particles
        double G4EDepMin_ = 0.0; // minimum energy deposition for all particles
        double G4ETrackingMin_ = 0.0; // minimum energy for tracking

        // our mc tree we will save energy depositions to
        I3MCTreePtr mcTree = nullptr;
        I3Particle primary_;

        // and photon summary that keeps track of optical photons
        // note this is in two parts -- map between track id and idx
        // and vector containing photon summaries (idx in vector corresponds to track id)
        boost::shared_ptr<I3Map<int, std::set<int>>> wls_parent_daughter_map = boost::make_shared<I3Map<int, std::set<int>>>(); // map between parent id and group of track ids
        int last_photon_summary_idx = -1;
        boost::shared_ptr<I3Map<int, PhotonSummary>> photon_summary = boost::make_shared<I3Map<int, PhotonSummary>>();

        // define a few things for converting energy to wavelength
        G4double const hbarc = 197.326; // eV * nm
        G4double const hc = hbarc * 2 * M_PI; // eV * nm

        static const std::unordered_map<std::string, int> energyLossToI3ParticlePDGCode;
        static const std::unordered_map<std::string, PhotonSummary::PhotonSource> processNameToPhotonSource;
        static const std::unordered_map<PhotonSummary::PhotonSource, std::string> photonSourceToProcessName;
        static const std::unordered_map<WLSLocation::WLSLoc, std::string> wlsLocationToProcessName;
        static const std::unordered_set<int> energyLossPDGCodes;
        static const std::unordered_set<std::string> knownProcessNames;

        std::map<int, I3ParticleID> DaughterParticleMap; // map between track_id and I3ParticleID
        std::map<int, int> parent_map;
        std::map<int, std::tuple<double, std::vector<I3Particle>>> sub_threshold_losses;

        std::map<int, size_t> num_children;

        double c_mm_per_nsec = 299.792458 * (I3Units::mm / I3Units::ns); // speed of light in mm/nsec

        std::vector<double> rindex_wavelength = {99.99909859236817, 109.99897296368657, 119.99889895379226, 123.99888385471958,
                                                 127.99884206872682, 133.99878562916103, 139.99873802931543, 159.99855258590853,
                                                 179.99837456270487, 199.99819718473634, 299.9973199734223, 399.99632984597423,
                                                 499.9954929618409, 599.9944947689613, 699.9936901465773}; // nm

        std::vector<double> rindex = {1.7898800000000006, 1.6199600000000003, 1.4500400000000002, 1.40284, 1.358,
                                      1.335167, 1.315166, 1.28071, 1.263128, 1.25475, 1.24, 1.23, 1.225, 1.222, 1.22};
};

#endif // G4CCMTreeTracker_H
