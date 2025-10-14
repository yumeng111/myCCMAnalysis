#ifndef G4INTERFACE_H
#define G4INTERFACE_H

#include <dataclasses/I3Map.h>
#include <dataclasses/I3Vector.h>
#include <dataclasses/physics/I3MCTree.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/physics/I3MCTreeUtils.h>

#include "g4-larsim/g4classes/G4CCMPrimaryGeneratorAction.h"
#include "g4-larsim/g4classes/G4CCMReadout.h"

#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>

#include <simclasses/CCMMCPE.h>
#include <simclasses/CCMSimulationSettings.h>
#include <simclasses/DetectorResponseConfig.h>

#include "G4MTRunManager.hh"

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
        void InstallDetector(CCMSimulationSettings const & settings, DetectorResponseConfig const & config);

        void SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries, I3MCTreePtr veto_tree, I3MCTreePtr inner_tree, I3VectorI3ParticlePtr veto_vector, I3VectorI3ParticlePtr inner_vector);
        void SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> veto_trees, std::vector<I3MCTreePtr> inner_trees, std::vector<I3VectorI3ParticlePtr> veto_vectors, std::vector<I3VectorI3ParticlePtr> inner_vectors);
        void SetNumberOfThreads(size_t n_threads) {
            n_cores_ = n_threads;
        }

    private:
        void Initialize();
        void InitializeRun();

        static std::shared_ptr<G4Interface> g4Interface_;

        size_t n_cores_ = 0;
        std::shared_ptr<G4CCMReadout> readout_;
        std::shared_ptr<G4MTRunManager> runManager_;
        std::shared_ptr<G4CCMParticleList> particle_list_;

        CCMSimulationSettings simulationSettings_;
        DetectorResponseConfig detectorConfig_;

        #ifdef G4VIS_USE
        G4VisManager* visManager_;
        #endif

        G4CCMDetectorConstruction * detector_ = nullptr;
        bool initialized_;
        std::string visMacro_;

        // controls to turn SD on/off (set by our response service)
        bool RecordHits_; // turn hit recording on/off

        SET_LOGGER("G4Interface");
};

#endif
