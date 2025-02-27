#ifndef G4CCMVetoSD_H
#define G4CCMVetoSD_H

#include "icetray/I3Units.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"
#include "g4-larsim/g4classes/G4CCMReadout.h"
#include "g4-larsim/g4classes/G4CCMTreeTracker.h"
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

class G4CCMVetoSD : public G4VSensitiveDetector {
    public:
        G4CCMVetoSD(G4String name);
        ~G4CCMVetoSD() override = default;

        void SetReadout(G4CCMReadout * readout) { readout_ = readout; }

        void Initialize(G4HCofThisEvent*) override;
        void EndOfEvent(G4HCofThisEvent*) override;
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*) override;

        void Reset() {}

        void SetSaveEnergyLossesVector(bool save) { SaveEnergyLossesVector_ = save; }
        void SetSaveEnergyLossesTree(bool save) { SaveEnergyLossesTree_ = save; }
        void SetPruneTree(bool prune) { PruneTree_ = prune; }
        bool GetSaveEnergyLossesVector() { return SaveEnergyLossesVector_; }
        bool GetSaveEnergyLossesTree() { return SaveEnergyLossesTree_; }
        bool GetPruneTree() { return PruneTree_; }

    private:

        G4CCMTreeTracker * tree_tracker;

        bool SaveEnergyLossesVector_ = false;
        bool SaveEnergyLossesTree_ = false;
        bool PruneTree_ = false;

        I3VectorI3ParticlePtr output_energy_losses_vector = nullptr;
        I3MCTreePtr output_energy_losses_tree = nullptr;

        std::set<I3ParticleID> energy_loss_ids;

        int event_id = -1;
        G4CCMReadout * readout_ = nullptr;
        G4CCMScintHitsCollection* fScintCollection = nullptr;
        G4int fHitsCID = -1;

        // and photon summary that keeps track of optical photons
        // note this is in two parts -- map between track id and idx
        // and vector containing photon summaries (idx in vector corresponds to track id)
        boost::shared_ptr<I3Map<int, size_t>> optical_photon_map = boost::make_shared<I3Map<int, size_t>>(); // map between track id and idx in vector
        boost::shared_ptr<I3Map<int, std::vector<int>>> wls_parent_daughter_map = boost::make_shared<I3Map<int, std::vector<int>>>(); // map between parent id and group of track ids
        PhotonSummarySeriesPtr photon_summary = boost::make_shared<PhotonSummarySeries>();
};

#endif // G4CCMVetoSD_H
