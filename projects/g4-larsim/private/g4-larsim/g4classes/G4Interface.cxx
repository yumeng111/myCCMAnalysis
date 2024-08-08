

#include "icetray/I3Logging.h"
#include "dataclasses/physics/I3MCTree.h"
#include "dataclasses/physics/I3Particle.h"
#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4Interface.h"
#include "g4-larsim/g4classes/G4CCMScintSD.h"
#include "dataclasses/physics/I3MCTreeUtils.h"
#include "g4-larsim/g4classes/G4CCMPhysicsList.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"

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
    detector_(NULL), initialized_(false), runInitialized_(false), visMacro_(visMacro) {

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


void G4Interface::InstallDetector(bool PMTSDStatus, bool LArSDStatus, bool SodiumSourceRun, double SodiumSourceLocation, double SingletTau, double TripletTau, double Rayleigh128,
                                  bool UVAbsStatus, bool TimeCut, bool CerenkovControl, long RandomSeed) {
    if(initialized_) {
        log_fatal("G4Interface already initialized. Cannot install detector!");
        return;
    }

    PMTSDStatus_ = PMTSDStatus;
    LArSDStatus_ = LArSDStatus;

    // add random seed for geant4
    G4Random::setTheSeed(RandomSeed);

    if (runManager_ == nullptr){
        runManager_ = std::make_shared<G4CCMRunManager>();
    }

    if(!detector_) {
        detector_ = new G4CCMDetectorConstruction(SingletTau / I3Units::nanosecond * CLHEP::ns, TripletTau / I3Units::nanosecond * CLHEP::ns, UVAbsStatus, Rayleigh128 / I3Units::cm * CLHEP::cm);
        // set SD status
        detector_->SetPMTSDStatus(PMTSDStatus_);
        detector_->SetLArSDStatus(LArSDStatus_);
        // set time cut and cerenkov control
        detector_->SetTimeCut(TimeCut);
        detector_->SetCerenkovControl(CerenkovControl);
        // set sodium rod status
        detector_->InitializeSodiumSourceRun(SodiumSourceRun, SodiumSourceLocation / I3Units::cm * CLHEP::cm);
        // Force reinitializatiion
        runManager_->ReinitializeGeometry(true);
    }

}

void G4Interface::InitializeRun()
{
    if(!initialized_) {
        Initialize();
    }

    if(!runInitialized_) {
        runManager_->InitializeRun();
        runInitialized_ = true;
    }
}


void G4Interface::InjectParticle(const I3Particle& particle, I3MCTreePtr edep_tree, CCMMCPESeriesMapPtr mcpeseries)
{
    mcpeseries_result_ = mcpeseries;

    if(!runInitialized_) {
        log_fatal("No run initialized. Cannot inject particle!");
        return;
    }

    // if we are tracking LAr energy deposition, let's pass on the primary particle information
    if (LArSDStatus_) {
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdNameScint = "/LAr/scintSD";
        G4CCMScintSD* scintSD = (G4CCMScintSD*) SDman->FindSensitiveDetector(sdNameScint);
        scintSD->SetPrimaryParticle(particle, edep_tree);
    }
    if (PMTSDStatus_) {
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdNamePMT = "/LAr/pmtSD";
        G4CCMPMTSD* pmtSD = (G4CCMPMTSD*) SDman->FindSensitiveDetector(sdNamePMT);
        pmtSD->Reset();
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDef = NULL;
    bool sodium_run = false;
    switch(particle.GetType())
    {
    case I3Particle::Gamma:
       particleDef = particleTable->FindParticle("opticalphoton");
       break;
    case I3Particle::EMinus:
       particleDef = particleTable->FindParticle("e-");
       break;
    case I3Particle::EPlus:
       particleDef = particleTable->FindParticle("e+");
       break;
    case I3Particle::MuMinus:
       particleDef = particleTable->FindParticle("mu-");
       break;
    case I3Particle::MuPlus:
       particleDef = particleTable->FindParticle("mu+");
       break;
    case I3Particle::PPlus:
       particleDef = particleTable->FindParticle("proton");
       break;
    case I3Particle::PMinus:
       particleDef = particleTable->FindParticle("anti_proton");
       break;
    case I3Particle::Neutron:
       particleDef = particleTable->FindParticle("neutron");
       break;
#ifdef I3PARTICLE_SUPPORTS_PDG_ENCODINGS
    case I3Particle::NeutronBar:
#else
    case 25:
#endif
       particleDef = particleTable->FindParticle("anti_neutron");
       break;
    case I3Particle::PiPlus:
       particleDef = particleTable->FindParticle("pi+");
       break;
    case I3Particle::PiMinus:
       particleDef = particleTable->FindParticle("pi-");
       break;
    case I3Particle::Pi0:
       particleDef = particleTable->FindParticle("pi0");
       break;
    case I3Particle::KPlus:
       particleDef = particleTable->FindParticle("kaon+");
       break;
    case I3Particle::KMinus:
       particleDef = particleTable->FindParticle("kaon-");
       break;
    case I3Particle::K0_Long:
       particleDef = particleTable->FindParticle("kaon0L");
       break;
    case I3Particle::K0_Short:
       particleDef = particleTable->FindParticle("kaon0S");
       break;
    case I3Particle::NuE:
       particleDef = particleTable->FindParticle("nu_e");
       break;
    case I3Particle::NuEBar:
       particleDef = particleTable->FindParticle("anti_nu_e");
       break;
    case I3Particle::NuMu:
       particleDef = particleTable->FindParticle("nu_mu");
       break;
    case I3Particle::NuMuBar:
       particleDef = particleTable->FindParticle("anti_nu_mu");
       break;
    case I3Particle::NuTau:
       particleDef = particleTable->FindParticle("nu_tau");
       break;
    case I3Particle::NuTauBar:
       particleDef = particleTable->FindParticle("anti_nu_tau");
       break;
    case I3Particle::Lambda:
       particleDef = particleTable->FindParticle("lambda");
       break;
    //~ new particles added!!!
    case I3Particle::LambdaBar:
       particleDef = particleTable->FindParticle("anti_lambda");
       break;
    case I3Particle::SigmaMinusBar:
       particleDef = particleTable->FindParticle("anti_sigma+");
       break;
    case I3Particle::Xi0Bar:
       particleDef = particleTable->FindParticle("anti_xi0");
       break;
    case I3Particle::Xi0:
       particleDef = particleTable->FindParticle("xi0");
       break;
     case I3Particle::SigmaPlus:
       particleDef = particleTable->FindParticle("sigma+");
       break;
    case I3Particle::SigmaMinus:
       particleDef = particleTable->FindParticle("sigma-");
       break;
    case I3Particle::XiMinus:
       particleDef = particleTable->FindParticle("xi-");
       break;
    case I3Particle::SigmaPlusBar:
       particleDef = particleTable->FindParticle("anti_sigma-");
       break;
    case I3Particle::XiPlusBar:
       particleDef = particleTable->FindParticle("anti_xi-");
       break;
    case I3Particle::OmegaPlusBar:
       particleDef = particleTable->FindParticle("anti_omega-");
       break;
    case I3Particle::H2Nucleus:
       particleDef = particleTable->FindParticle("deuteron");
       break;
    case I3Particle::H3Nucleus:
       particleDef = particleTable->FindParticle("triton");
       break;
    case I3Particle::He3Nucleus:
       particleDef = particleTable->FindParticle("He3");
       break;
    case I3Particle::He4Nucleus:
       particleDef = particleTable->FindParticle("He4");
       break;
    case I3Particle::He5Nucleus:
       particleDef = particleTable->FindParticle("He5");
       break;
    case I3Particle::He6Nucleus:
       particleDef = particleTable->FindParticle("He6");
       break;
    case I3Particle::Na22Nucleus:
       sodium_run = true;
       break;
    default:
      log_warn("Man, check out that strange particle \"%s\" ?!", particle.GetTypeString().c_str());
      return;
    }

    G4ParticleGun gun(1);

    // special logic for adding sodium particle
    if (sodium_run){

        G4int Z = 11, A = 22;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;

        particleDef = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
        gun.SetParticleCharge(ionCharge);
    }

    if (!particleDef){
        log_warn("You passed NULL particleDef \"%s\" ?!", particle.GetTypeString().c_str());
        return;
    }


    // Particle position in G4 units
    G4ThreeVector position((particle.GetX() / I3Units::m) * CLHEP::m,
                           (particle.GetY() / I3Units::m) * CLHEP::m,
                           (particle.GetZ() / I3Units::m) * CLHEP::m);

    G4ThreeVector direction(particle.GetDir().GetX(),
                            particle.GetDir().GetY(),
                            particle.GetDir().GetZ());

    gun.SetParticleDefinition(particleDef);
    gun.SetParticlePosition(position);
    gun.SetParticleEnergy((particle.GetEnergy() / I3Units::MeV) * CLHEP::MeV);
    gun.SetParticleMomentumDirection(direction);

    log_trace("Injecting %s: x=%.2f m, y=%.2f m, z=%.2f m, E=%.3f MeV",
              particle.GetTypeString().c_str(),
              position.x() / CLHEP::m,
              position.y() / CLHEP::m,
              position.z() / CLHEP::m,
              gun.GetParticleEnergy() / CLHEP::MeV);

    runManager_->InjectParticle(&gun);
}


void G4Interface::TerminateEvent() {
    // let's grab the CCMMCPE map from G4Interface
    // now let's grab SD information
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap;
    boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map;
    PhotonSummarySeriesPtr photon_summary_series;

    if (PMTSDStatus_) {
        G4String sdNamePMT = "/LAr/pmtSD";
        G4CCMPMTSD* pmtSD = (G4CCMPMTSD*) SDman->FindSensitiveDetector(sdNamePMT);
        CCMMCPEMap = pmtSD->GetCCMMCPEMap();
    }

    if (LArSDStatus_) {
        G4String sdNameScint = "/LAr/scintSD";
        G4CCMScintSD* scintSD = (G4CCMScintSD*) SDman->FindSensitiveDetector(sdNameScint);
        photon_summary_series = scintSD->GetPhotonSummarySeries();
        photon_summary_series_map = scintSD->GetPhotonSummaryMap();
    }

    if(PMTSDStatus_ and LArSDStatus_) {
        MergeMCPESeries(mcpeseries_result_, CCMMCPEMap, photon_summary_series_map, photon_summary_series);
    } else if(PMTSDStatus_) {
        MergeMCPESeries(mcpeseries_result_, CCMMCPEMap);
    }
}

void G4Interface::MergeMCPESeries(CCMMCPESeriesMapPtr mcpeseries_dest, CCMMCPESeriesMapPtr mcpeseries_source, boost::shared_ptr<I3Map<int, size_t>> photon_summary_series_map, PhotonSummarySeriesPtr photon_summary_series) {
    // Iterate over PMTs in source map
    for (CCMMCPESeriesMap::iterator it = mcpeseries_source->begin(); it != mcpeseries_source->end(); ++it) {

        // Find the corresponding PMT in the destination map
        CCMMCPESeriesMap::iterator it_dest = mcpeseries_dest->find(it->first);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == mcpeseries_dest->end()) {
            mcpeseries_dest->insert(std::make_pair(it->first, CCMMCPESeries()));
            // Update the iterator so it points to our new entry
            it_dest = mcpeseries_dest->find(it->first);
        }

        // Reference to the destination
        CCMMCPESeries & dest_series = it_dest->second;

        // Iterate over the vector of CCMMCPE in the source map for this PMT
        for (CCMMCPESeries::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            CCMMCPE const & pe = *it2;

            // Check if the photon has everything properly recorded
            I3Map<int, size_t>::iterator it_map = photon_summary_series_map->find(pe.track_id);
            // If not trash it as an edge case
            if(it_map == photon_summary_series_map->end())
                continue;

            // Grab the summary information for this photon track
            PhotonSummary const & this_photon_summary = photon_summary_series->at(it_map->second);

            // Copy the CCMMCPE from source to destination
            dest_series.emplace_back(pe);
            // Grab a reference to the CCMMCPE in the destination
            CCMMCPE & dest = dest_series.back();

            // Update the destination CCMMCPE with the summary information
            dest.g4_time = this_photon_summary.g4_time;
            dest.calculated_time = this_photon_summary.calculated_time;
            dest.g4_distance_uv = this_photon_summary.g4_distance_uv;
            dest.g4_distance_visible = this_photon_summary.g4_distance_visible;
            dest.calculated_distance_uv = this_photon_summary.calculated_distance_uv;
            dest.calculated_distance_visible = this_photon_summary.calculated_distance_visible;
        }
    }
}

void G4Interface::MergeMCPESeries(CCMMCPESeriesMapPtr mcpeseries_dest, CCMMCPESeriesMapPtr mcpeseries_source) {
    // Iterate over PMTs in source map
    for (CCMMCPESeriesMap::iterator it = mcpeseries_source->begin(); it != mcpeseries_source->end(); ++it) {

        // Find the corresponding PMT in the destination map
        CCMMCPESeriesMap::iterator it_dest = mcpeseries_dest->find(it->first);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == mcpeseries_dest->end()) {
            mcpeseries_dest->insert(std::make_pair(it->first, CCMMCPESeries()));
            // Update the iterator so it points to our new entry
            it_dest = mcpeseries_dest->find(it->first);
        }

        // Reference to the destination
        CCMMCPESeries & dest_series = it_dest->second;

        // Iterate over the vector of CCMMCPE in the source map for this PMT
        for (CCMMCPESeries::iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            CCMMCPE const & pe = *it2;
            // Copy the CCMMCPE from source to destination
            dest_series.emplace_back(pe);
        }
    }
}

void G4Interface::TerminateRun()
{
    runManager_->TerminateRun();
    runInitialized_ = false;
}

void G4Interface::Initialize()
{
    if(initialized_) {
        log_error("G4Interface has already been initialized. Ignoring this call!");
        return;
    }

    // set number of threads
    //int num_threads = std::thread::hardware_concurrency();
    //runManager_.SetNumberOfThreads(num_threads);

    // Set verbosity
    int32_t verboseLevel = 0;

    log_debug("Init geometry ...");
    runManager_->SetUserInitialization(detector_);

    log_debug("Init physics list ...");

    // adding physics list
    G4CCMPhysicsList* physics_list = new G4CCMPhysicsList(verboseLevel);
    runManager_->SetUserInitialization(physics_list);

    // Initialize G4 kernel
    log_debug("Init run manager ...");
    runManager_->Initialize();


    switch (GetIcetrayLogger()->LogLevelForUnit("G4Interface"))
    {
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
    if(!visMacro_.empty())
    {
    G4UImanager* uim = G4UImanager::GetUIpointer();

    // Checking geometry
    uim->ApplyCommand("/geometry/test/grid_test");

    // Execute visualization macro
    std::string visCmd = "/control/execute " + visMacro_;
    uim->ApplyCommand(visCmd.c_str());
    }

    initialized_ = true;
}
