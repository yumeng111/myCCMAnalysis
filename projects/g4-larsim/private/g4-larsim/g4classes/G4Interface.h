#ifndef G4INTERFACE_H
#define G4INTERFACE_H

#include <dataclasses/I3Map.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

#include <g4-larsim/g4classes/G4CCMRunManager.h>

#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>

#include <simclasses/CCMMCPE.h>
#include "simclasses/PhotonSummary.h"

#ifdef G4VIS_USE
class G4VisManager;
#endif

class I3Particle;
class G4CCMDetectorConstruction;

/**
 * Top-level class to handle Geant4. All global things are initialized here (run manager, visualization manager, detector construction, physics list and user actions).
 */

class G4Interface {
    private:
        G4Interface(const std::string& visMacro="");
    public:
        ~G4Interface();

        //  Static method which returns the singleton pointer to this class
        static std::shared_ptr<G4Interface> GetInstance() {
            if(g4Interface_ == nullptr) {
                g4Interface_ = std::shared_ptr<G4Interface>(new G4Interface());
            }
            return g4Interface_;
        };
        static void DestroyInstance() {
            g4Interface_ = std::shared_ptr<G4Interface>(nullptr);
        };

        /// Add the detector to the geometry. Should not be called after initialized.
        void InstallDetector(bool PMTSDStatus, bool LArSDStatus, bool SourceRodIn, double SourceRodLocation, bool CobaltSourceRun, bool SodiumSourceRun,
                             double SingletTau, double TripletTau, double Rayleigh128,
                             double UVAbsLength, bool TimeCut, bool KillCherenkov, long RandomSeed);
        /// Initialize event. Most Geant4 global things are initialized the first time this is called.
        void InitializeRun();
        /// To be called after simulating each IceTray event.
        void TerminateEvent();
        void TerminateRun();
        /// Simulate a single particle (InitializeRun must be called first)
        void InjectParticle(const I3Particle& particle, I3MCTreePtr edep_tree, CCMMCPESeriesMapPtr mcpeseries);

        static void MergeMCPESeries(CCMMCPESeriesMapPtr mcpeseries_dest, CCMMCPESeriesMapPtr mcpeseries_source, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map, PhotonSummarySeriesPtr photon_summary_series);
        static void MergeMCPESeries(CCMMCPESeriesMapPtr mcpeseries_dest, CCMMCPESeriesMapPtr mcpeseries_source);

        boost::shared_ptr<CCMMCPESeriesMap> GetCCMMCPEMap() { return mcpeseries_result_ ; }

    private:
        void Initialize();

        static std::shared_ptr<G4Interface> g4Interface_;

        std::shared_ptr<G4CCMRunManager> runManager_ = nullptr;

        #ifdef G4VIS_USE
        G4VisManager* visManager_;
        #endif

        CCMMCPESeriesMapPtr mcpeseries_result_;

        G4CCMDetectorConstruction* detector_;
        bool initialized_;
        bool runInitialized_;
        std::string visMacro_;

        // controls to turn SD on/off (set by our response service)
        bool PMTSDStatus_; // turn PMT SD on/off
        bool LArSDStatus_; // turn fiducial LAr SD on/off

        static const std::unordered_map<PhotonSummary::PhotonSource, CCMMCPE::PhotonSource> PhotonSummarytoCCMMCPEPhotonSource;

        SET_LOGGER("G4Interface");
};

#endif
