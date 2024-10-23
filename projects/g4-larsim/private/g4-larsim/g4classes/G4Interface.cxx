

#include "icetray/I3Logging.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4Interface.h"
#include "g4-larsim/g4classes/G4CCMScintSD.h"
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


void G4Interface::InstallDetector(bool PMTSDStatus, bool LArSDStatus, bool SourceRodIn, double SourceRodLocation, bool CobaltSourceRun, bool SodiumSourceRun, 
                                  double SingletTau, double TripletTau, double Rayleigh128, double UVAbsLength,
                                  double WLSNPhotonsEndCapFoil, double WLSNPhotonsSideFoil, double WLSNPhotonsPMT, 
                                  double EndCapFoilTPBThickness, double SideFoilTPBThickness, double PMTTPBThickness, 
                                  double TPBAbsTau, double TPBAbsNorm, double TPBAbsScale, double Mie_GG, double Mie_Ratio,
                                  bool TimeCut, bool KillCherenkov, bool FullPhotonTracking, long RandomSeed) {
    if(initialized_) {
        log_fatal("G4Interface already initialized. Cannot install detector!");
        return;
    }

    PMTSDStatus_ = PMTSDStatus;
    LArSDStatus_ = LArSDStatus;

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
                                                  UVAbsLength / I3Units::cm * CLHEP::cm, WLSNPhotonsEndCapFoil, WLSNPhotonsSideFoil, WLSNPhotonsPMT,
                                                  EndCapFoilTPBThickness / I3Units::mm * CLHEP::mm, SideFoilTPBThickness / I3Units::mm * CLHEP::mm, PMTTPBThickness / I3Units::mm * CLHEP::mm,
                                                  Rayleigh128 / I3Units::cm * CLHEP::cm, TPBAbsTau, TPBAbsNorm, TPBAbsScale, Mie_GG, Mie_Ratio);
        // set readout
        detector_->SetReadout(readout_.get());
        // set SD status
        detector_->SetPMTSDStatus(PMTSDStatus_);
        detector_->SetLArSDStatus(LArSDStatus_);
        // set time cut and cerenkov control
        detector_->SetTimeCut(TimeCut);
        detector_->SetKillCherenkov(KillCherenkov);
        detector_->SetPhotonTracking(FullPhotonTracking);
        // set sodium rod status
        detector_->InitializeSodiumSourceRun(SourceRodIn, SourceRodLocation / I3Units::cm * CLHEP::cm, CobaltSourceRun, SodiumSourceRun);
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

void G4Interface::SimulateEvent(const I3Particle& particle, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) {
    InitializeRun();

    std::vector<I3Particle> particles = {particle};
    std::vector<I3MCTreePtr> trees = {tree};
    std::vector<CCMMCPESeriesMapPtr> mcpeseries_list = {mcpeseries};

    // Set the particle list used for the primary generator user action
    particle_list_->SetParticles(particles);

    // Set readout information for sensitive detectors
    readout_->SetInput(particles, mcpeseries_list, trees);

    // Run the event
    runManager_->BeamOn(1);

    // I3MCTree and CCMMCPESeriesMap were passed as shared pointers to readout_,
    // and so have already been updated
}

void G4Interface::SimulateEvents(std::vector<I3Particle> const & particles, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries) {
    InitializeRun();

    // Set the particle list used for the primary generator user action
    particle_list_->SetParticles(particles);

    // Set readout information for sensitive detectors
    readout_->SetInput(particles, mcpeseries, trees);

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
