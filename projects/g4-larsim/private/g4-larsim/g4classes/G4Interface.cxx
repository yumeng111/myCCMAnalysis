
#include <g4-larsim/g4classes/G4Interface.h>
#include <g4-larsim/g4classes/G4CCMDetectorConstruction.h>
#include <g4-larsim/g4classes/G4CCMPhysicsList.h>
#include <g4-larsim/g4classes/G4CCMUserTrackingAction.h>
#include <g4-larsim/g4classes/G4CCMUserSteppingAction.h>
#include <g4-larsim/g4classes/G4CCMPMTSD.h>

#include <icetray/I3Logging.h>
#include <dataclasses/physics/I3Particle.h>

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4UImanager.hh>
#include <G4SDManager.hh>

G4Interface* G4Interface::g4Interface_ = NULL;

G4Interface::G4Interface(const std::string& visMacro):
    detector_(NULL), initialized_(false), eventInitialized_(false), visMacro_(visMacro) {
    g4Interface_ = this;

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
    g4Interface_ = NULL;

    #ifdef G4VIS_USE
    if(visManager_) delete visManager_;
    #endif
}


void G4Interface::InstallDetector() {
    if(initialized_) {
        log_fatal("G4Interface already initialized. Cannot install detector!");
        return;
    }

    if(!detector_) {
        detector_ = new G4CCMDetectorConstruction();
    }
    
}

void G4Interface::InitializeEvent()
{
    if(!initialized_) {
        Initialize();
    }

    if(!eventInitialized_) {
        runManager_.InitializeRun();
        eventInitialized_ = true;
    }
}


void G4Interface::InjectParticle(const I3Particle& particle)
{
    if(!eventInitialized_) {
        log_fatal("No event initialized. Cannot inject particle!");
        return;
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDef = NULL;
    switch(particle.GetType())
    {
    case I3Particle::Gamma:
       particleDef = particleTable->FindParticle("gamma");
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
    default:
      log_warn("Man, check out that strange particle \"%s\" ?!", particle.GetTypeString().c_str());
      return;
    }
  
    // Particle position in G4 units
    G4ThreeVector position((particle.GetX() / I3Units::m) * CLHEP::m,
                           (particle.GetY() / I3Units::m) * CLHEP::m,
                           (particle.GetZ() / I3Units::m) * CLHEP::m);

    // Transform I3 coordinates to world system
    //position -= detector_->GetWorldOrigin();

    G4ThreeVector direction(particle.GetDir().GetX(),
                            particle.GetDir().GetY(),
                            particle.GetDir().GetZ());

    if (!particleDef){
        log_warn("You passed NULL particleDef \"%s\" ?!", particle.GetTypeString().c_str());
        return;
    }

    G4ParticleGun gun(1);
    gun.SetParticleDefinition(particleDef);
    gun.SetParticleEnergy((particle.GetEnergy() / I3Units::MeV) * CLHEP::MeV);
    gun.SetParticlePosition(position);
    gun.SetParticleMomentumDirection(direction);

    log_trace("Injecting %s: x=%.2f m, y=%.2f m, z=%.2f m, E=%.3f MeV",
              particle.GetTypeString().c_str(),
              position.x() / CLHEP::m,
              position.y() / CLHEP::m,
              position.z() / CLHEP::m,
              gun.GetParticleEnergy() / CLHEP::MeV);

    runManager_.InjectParticle(&gun);
}


void G4Interface::TerminateEvent()
{
    // call sd function to get vector of ccmmmcpe
    if(eventInitialized_) {
        // ok so we are done with this event
        // let's grab the map between CCMPMTKey and std::vector<CCMMCPE> from out sensitve detector to save to the frame
        G4SDManager* sdManager = G4SDManager::GetSDMpointer();
        G4CCMPMTSD* PMTSD = dynamic_cast<G4CCMPMTSD*>(sdManager->FindSensitiveDetector("/LAr/pmtSD"));
        if (PMTSD) {
            CCMMCPEMap = PMTSD->GetCCMMCPEMap();
        }

        runManager_.TerminateRun();
        eventInitialized_ = false;
    }
}


void G4Interface::Initialize()
{
    if(initialized_) {
        log_error("G4Interface has already been initialized. Ignoring this call!");
        return;
    }

    log_debug("Init geometry ...");
    runManager_.SetUserInitialization(detector_);

    log_debug("Init physics list ...");
    //runManager_.SetUserInitialization(new G4CCMPhysicsList(useScintillator_));

    log_debug("Init UserTrackingAction ...");
    runManager_.SetUserAction(new G4CCMUserTrackingAction());

    log_debug("Init UserSteppingAction ...");
    runManager_.SetUserAction(new G4CCMUserSteppingAction());

    // Initialize G4 kernel
    log_debug("Init run manager ...");
    runManager_.Initialize();

    // Set verbosity
    int32_t verboseLevel = 0;
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

    runManager_.SetVerboseLevel(verboseLevel);
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
