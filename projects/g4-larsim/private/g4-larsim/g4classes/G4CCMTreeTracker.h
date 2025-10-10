#ifndef G4CCMTreeTracker_H
#define G4CCMTreeTracker_H

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMReadout.h"
#include "dataclasses/I3Map.h"
#include "simclasses/PhotonSummary.h"
#include "simclasses/CCMSimulationSettings.h"
#include "simclasses/DetectorResponseConfig.h"

#include <tuple>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <G4SystemOfUnits.hh>
#include <G4VSensitiveDetector.hh>

class G4Step;
class G4HCofThisEvent;
class G4CCMEDepSD;

class G4CCMTreeTracker : public G4VSensitiveDetector {
    typedef G4CCMReadout::ParentInfo ParentInfo;
    public:
        friend G4CCMEDepSD;
        G4CCMTreeTracker(G4String name, bool trackParticles = false, bool trackEnergyLosses = false);
        ~G4CCMTreeTracker() override = default;

        void SetReadout(G4CCMReadout * readout) {readout_ = readout;}

        void Initialize(G4HCofThisEvent*) override;
        void EndOfEvent(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;

        void SetDetectorResponseConfig(DetectorResponseConfig const & config) {
            detectorConfig_ = config;
            detectorConfig_.to_g4_units();
        }

        void SetSimulationSettings(CCMSimulationSettings const & settings) {
            simulationSettings_ = settings;
            simulationSettings_.to_g4_units();
        }

        bool GetTrackParticles() const { return TrackParticles_; }
        bool GetTrackEnergyLosses() const { return TrackEnergyLosses_; }

        void Reset() {
            DaughterParticleMap.clear();
            summary_map = boost::make_shared<I3Map<int, std::tuple<ParentInfo, PhotonSummary>>>();
            source_map = boost::make_shared<I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>>();
            wls_parent_daughter_map = boost::make_shared<I3Map<int, std::set<int>>>();
            sub_threshold_losses.clear();
            parent_map.clear();
            num_children.clear();
            num_siblings.clear();
            mcTree = nullptr;
        }

        void AddNewPhoton(int parent_id, int track_id, double time, double distance, double wavelength, std::string creation_process_name);
        void AddPhotonTrack(int parent_id, int track_id, double delta_time, double delta_distance, std::string creation_process_name);
        void UpdatePhoton(int parent_id, int track_id, double delta_time, double delta_distance, std::string creation_process_name);
        void AddWLSPhotonTrack(int parent_id, int track_id, double delta_time, double delta_distance, G4Step* aStep);


    private:
        int event_id = -1;
        G4CCMReadout * readout_ = nullptr;

        DetectorResponseConfig detectorConfig_;
        CCMSimulationSettings simulationSettings_;

        bool TrackParticles_ = false;
        bool TrackEnergyLosses_ = false;

        // our mc tree we will save energy depositions to
        I3MCTreePtr mcTree = nullptr;
        I3Particle primary_;

        // and photon summary that keeps track of optical photons
        // note this is in two parts -- map between track id and idx
        // and vector containing photon summaries (idx in vector corresponds to track id)
        boost::shared_ptr<I3Map<int, std::set<int>>> wls_parent_daughter_map = boost::make_shared<I3Map<int, std::set<int>>>(); // map between parent id and group of track ids

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
        boost::shared_ptr<I3Map<int, std::tuple<ParentInfo, PhotonSummary>>> summary_map = boost::make_shared<I3Map<int, std::tuple<ParentInfo, PhotonSummary>>>();
        boost::shared_ptr<I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>> source_map = boost::make_shared<I3Map<int, std::tuple<ParentInfo, PhotonSummary::PhotonSource>>>();
        std::map<int, std::tuple<double, std::vector<I3Particle>>> sub_threshold_losses;

        std::map<int, size_t> num_children;
        std::map<int, std::shared_ptr<size_t>> num_siblings;

        double c_mm_per_nsec = 299.792458 * (I3Units::mm / I3Units::ns); // speed of light in mm/nsec

        std::vector<double> rindex_wavelength = {99.99909859236817, 109.99897296368657, 119.99889895379226, 123.99888385471958,
                                                 127.99884206872682, 133.99878562916103, 139.99873802931543, 159.99855258590853,
                                                 179.99837456270487, 199.99819718473634, 299.9973199734223, 399.99632984597423,
                                                 499.9954929618409, 599.9944947689613, 699.9936901465773}; // nm

        std::vector<double> rindex = {1.7898800000000006, 1.6199600000000003, 1.4500400000000002, 1.40284, 1.358,
                                      1.335167, 1.315166, 1.28071, 1.263128, 1.25475, 1.24, 1.23, 1.225, 1.222, 1.22};
};

#endif // G4CCMTreeTracker_H
