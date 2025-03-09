

#include "icetray/I3Logging.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4Interface.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "g4-larsim/g4classes/G4CCMPhysicsList.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"
#include "g4-larsim/g4classes/G4CCMActionInitialization.h"

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#include <G4Run.hh>
#include <G4Event.hh>
#include <FTFP_BERT.hh>
#include <G4IonTable.hh>
#include <G4UImanager.hh>
#include <G4SDManager.hh>
#include <G4RunManager.hh>
#include <G4ParticleGun.hh>
#include <G4DecayPhysics.hh>
#include <G4EventManager.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4OpticalPhysics.hh>
#include <G4OpticalParameters.hh>
#include <G4ParticleDefinition.hh>
#include <G4VModularPhysicsList.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4DecayTable.hh>
#include <G4RadioactiveDecay.hh>
#include <G4BetaPlusDecay.hh>
#include <G4GenericIon.hh>
#include <G4Geantino.hh>

std::shared_ptr<G4Interface> G4Interface::g4Interface_ = std::shared_ptr<G4Interface>(nullptr);

G4Interface::G4Interface(const std::string& visMacro):
    detector_(NULL), initialized_(false), visMacro_(visMacro) {

    // Visualization manager
    #ifdef G4VIS_USE
    visManager_ = NULL;
    if(!visMacro_.empty()) {
        visManager_ = new G4VisExecutive();
        visManager_->Initialize();
    }
    #endif
}


G4Interface::~G4Interface() {
    #ifdef G4VIS_USE
    if(visManager_) delete visManager_;
    #endif
}


void G4Interface::InstallDetector(
                                  bool SaveAllEnergyLossesTree,
                                  bool VetoSDSaveEnergyLossesVector, bool VetoSDSaveEnergyLossesTree, bool VetoSDPruneTree,
                                  bool InteriorSDSaveEnergyLossesVector, bool InteriorSDSaveEnergyLossesTree, bool InteriorSDPruneTree,
                                  bool KillNeutrinos, bool KillPhotons, bool KillScintillation, bool KillCherenkov,
                                  bool TimeCut, bool DetailedPhotonTracking, bool TrackParticles, bool TrackEnergyLosses,
                                  bool RecordHits, bool SourceRodIn, double SourceRodLocation,
                                  bool CobaltSourceRun, bool SodiumSourceRun, bool TrainingSource, 
                                  double DecayX, double DecayY, double DecayZ,
                                  double SingletTau, double TripletTau, double Rayleigh128,
                                  double UVAbsLength1, double UVAbsLength2, double UVAbsScaling, 
                                  double WLSNPhotonsEndCapFoil, double WLSNPhotonsSideFoil, double WLSNPhotonsPMT, 
                                  double EndCapFoilTPBThickness, double SideFoilTPBThickness, double PMTTPBThickness, 
                                  double TPBAbsTau, double TPBAbsNorm, double TPBAbsScale, double Mie_GG, double Mie_Ratio,
                                  double Normalization, double PhotonSampling, long RandomSeed) {
    if(initialized_) {
        log_fatal("G4Interface already initialized. Cannot install detector!");
        return;
    }

    RecordHits_ = RecordHits;

    // add random seed for geant4
    G4Random::setTheSeed(RandomSeed);

    if(readout_ == nullptr) {
        readout_ = std::make_shared<G4CCMReadout>();
    }

    if(runManager_ == nullptr) {
        runManager_ = std::make_shared<G4MTRunManager>();
    }

    if(particle_list_ == nullptr) {
        particle_list_ = std::make_shared<G4CCMParticleList>();
    }

    if(detector_ == nullptr) {
        detector_ = new G4CCMDetectorConstruction(SingletTau / I3Units::nanosecond * CLHEP::ns, TripletTau / I3Units::nanosecond * CLHEP::ns,
                                                  UVAbsLength1 / I3Units::cm * CLHEP::cm, UVAbsLength2 / I3Units::cm * CLHEP::cm, UVAbsScaling,
                                                  WLSNPhotonsEndCapFoil, WLSNPhotonsSideFoil, WLSNPhotonsPMT,
                                                  EndCapFoilTPBThickness / I3Units::mm * CLHEP::mm, SideFoilTPBThickness / I3Units::mm * CLHEP::mm, PMTTPBThickness / I3Units::mm * CLHEP::mm,
                                                  Rayleigh128 / I3Units::cm * CLHEP::cm, TPBAbsTau, TPBAbsNorm, TPBAbsScale, Mie_GG, Mie_Ratio, Normalization, PhotonSampling);
        // set readout
        detector_->SetReadout(readout_.get());

        detector_->SetSaveAllEnergyLossesTree(SaveAllEnergyLossesTree);

        detector_->SetVetoSDSaveEnergyLossesVector(VetoSDSaveEnergyLossesVector);
        detector_->SetVetoSDSaveEnergyLossesTree(VetoSDSaveEnergyLossesTree);
        detector_->SetVetoSDPruneTree(VetoSDPruneTree);

        detector_->SetInteriorSDSaveEnergyLossesVector(InteriorSDSaveEnergyLossesVector);
        detector_->SetInteriorSDSaveEnergyLossesTree(InteriorSDSaveEnergyLossesTree);
        detector_->SetInteriorSDPruneTree(InteriorSDPruneTree);

        // set SD status
        detector_->SetRecordHits(RecordHits_);

        // set time cut and cerenkov control
        detector_->SetTimeCut(TimeCut);
        detector_->SetKillNeutrinos(KillNeutrinos);
        detector_->SetKillPhotons(KillPhotons);
        detector_->SetKillCherenkov(KillCherenkov);
        detector_->SetKillScintillation(KillScintillation);
        detector_->SetDetailedPhotonTracking(DetailedPhotonTracking);
        // set sodium rod status
        detector_->InitializeSodiumSourceRun(SourceRodIn, SourceRodLocation / I3Units::cm * CLHEP::cm, CobaltSourceRun, SodiumSourceRun,
                                             TrainingSource, DecayX / I3Units::cm * CLHEP::cm, DecayY / I3Units::cm * CLHEP::cm, DecayZ / I3Units::cm * CLHEP::cm);
        // Force reinitializatiion
        //runManager_->ReinitializeGeometry(true);
    }
}

void G4Interface::InitializeRun() {
    if(!initialized_) {
        Initialize();
    }

    readout_->Reset();
}

//void SetOutputs(std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> edep_trees, std::vector<I3MCTreePtr> veto_edep_trees, std::vector<I3MCTreePtr> inner_edep_trees, std::vector<I3VectorI3ParticlePtr> veto_edep_vector, std::vector<I3VectorI3ParticlePtr> inner_edep_vector);

void G4Interface::SimulateEvent(const I3Particle& particle, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries, I3MCTreePtr veto_tree, I3MCTreePtr inner_tree, I3VectorI3ParticlePtr veto_vector, I3VectorI3ParticlePtr inner_vector) {
    InitializeRun();

    std::vector<I3Particle> particles = {particle};

    std::vector<CCMMCPESeriesMapPtr> mcpeseries_list = {mcpeseries};
    std::vector<I3MCTreePtr> trees = {tree};
    std::vector<I3MCTreePtr> veto_trees = {veto_tree};
    std::vector<I3MCTreePtr> inner_trees = {inner_tree};
    std::vector<I3VectorI3ParticlePtr> veto_vectors = {veto_vector};
    std::vector<I3VectorI3ParticlePtr> inner_vectors = {inner_vector};

    // Set the particle list used for the primary generator user action
    particle_list_->SetParticles(particles);

    // Set readout information for sensitive detectors
    readout_->SetInput(particles);
    readout_->SetOutputs(mcpeseries_list, trees, veto_trees, inner_trees, veto_vectors, inner_vectors);

    // Run the event
    runManager_->BeamOn(1);

    // I3MCTree and CCMMCPESeriesMap were passed as shared pointers to readout_,
    // and so have already been updated
}

void G4Interface::SimulateEvents(std::vector<I3Particle> const & particles, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries, std::vector<I3MCTreePtr> veto_trees, std::vector<I3MCTreePtr> inner_trees, std::vector<I3VectorI3ParticlePtr> veto_vectors, std::vector<I3VectorI3ParticlePtr> inner_vectors) {
    InitializeRun();

    // Set the particle list used for the primary generator user action
    particle_list_->SetParticles(particles);

    // Set readout information for sensitive detectors
    readout_->SetInput(particles);
    readout_->SetOutputs(mcpeseries, trees, veto_trees, inner_trees, veto_vectors, inner_vectors);

    // Run the event
    runManager_->BeamOn(particles.size());

    // I3MCTree and CCMMCPESeriesMap were passed as shared pointers to readout_,
    // and so have already been updated
}

void G4Interface::Initialize() {
    if(initialized_) {
        log_error("G4Interface has already been initialized. Ignoring this call!");
        return;
    }

    // set number of threads
    if(n_cores_ == 0) {
        n_cores_ = std::thread::hardware_concurrency();
    }
    runManager_->SetNumberOfThreads(n_cores_);
    readout_->SetNumberOfThreads(n_cores_);

    // Set verbosity
    int32_t verboseLevel = 0;

    log_debug("Init geometry ...");
    runManager_->SetUserInitialization(detector_);

    log_debug("Init physics list ...");

    // adding physics list
    G4CCMPhysicsList* physics_list = new G4CCMPhysicsList(verboseLevel);
    runManager_->SetUserInitialization(physics_list);

    // Set user action initialization
    G4CCMActionInitialization* actionInitialization = new G4CCMActionInitialization(particle_list_.get());
    runManager_->SetUserInitialization(actionInitialization);

    // Initialize G4 kernel
    log_debug("Init run manager ...");
    runManager_->Initialize();


    switch (GetIcetrayLogger()->LogLevelForUnit("G4Interface")) {
        case I3LOG_FATAL:
        case I3LOG_ERROR:
        case I3LOG_WARN:
        case I3LOG_INFO:
        case I3LOG_NOTICE:
        default:
            verboseLevel = 0;
            break;
        case I3LOG_DEBUG:
            verboseLevel = 1;
            break;
        case I3LOG_TRACE:
            verboseLevel = 2;
            break;
    }

    runManager_->SetVerboseLevel(verboseLevel);
    G4EventManager::GetEventManager()->SetVerboseLevel(verboseLevel);
    G4EventManager::GetEventManager()->GetStackManager()->SetVerboseLevel(verboseLevel);
    G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(verboseLevel);
#ifdef G4VIS_USE
    if(visManager_) visManager_->SetVerboseLevel(verboseLevel);
#endif

    // Execute visualization macro (if specified)
    if(!visMacro_.empty()) {
        G4UImanager* uim = G4UImanager::GetUIpointer();

        // Checking geometry
        uim->ApplyCommand("/geometry/test/grid_test");

        // Execute visualization macro
        std::string visCmd = "/control/execute " + visMacro_;
        uim->ApplyCommand(visCmd.c_str());
    }

    initialized_ = true;
}
