#include <cmath>

#include "g4-larsim/g4classes/G4CCMPrimaryGeneratorAction.h"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonTable.hh"

G4CCMParticleList::G4CCMParticleList() : particle_list_(), current_particle_ (0) {
}

void G4CCMParticleList::AddParticle(I3Particle p) {
    particle_list_.push_back(p);
}

void G4CCMParticleList::AddParticles(std::vector<I3Particle> p) {
    particle_list_.insert(particle_list_.end(), p.begin(), p.end());
}

void G4CCMParticleList::SetParticles(std::vector<I3Particle> p) {
    particle_list_ = p;
    current_particle_ = 0;
}

void G4CCMParticleList::Clear() {
    particle_list_.clear();
    current_particle_ = 0;
}

I3Particle G4CCMParticleList::GetNextParticle() {
    std::lock_guard<std::mutex> const lock(mutex);
    if(current_particle_ >= particle_list_.size())
        return I3Particle();
    if(current_particle_ >= particle_list_.size())
        log_fatal("current_particle_ >= particle_list_.size()");
    return particle_list_[current_particle_++];
}

G4CCMPrimaryGeneratorAction::G4CCMPrimaryGeneratorAction(G4CCMParticleList * particle_list) {
    particle_list_ = particle_list;
}

void G4CCMPrimaryGeneratorAction::GeneratePrimaries(G4Event * evt) {
    I3Particle p = particle_list_->GetNextParticle();

    if( std::isnan(p.GetTime()) ) {
        p.SetTime(0.);
    }

    G4PrimaryVertex * vertex = CreatePrimaryVertex(p);

    evt->AddPrimaryVertex(vertex);
}

G4PrimaryVertex * CreatePrimaryVertex(I3Particle const & particle) {


    G4double x = (particle.GetX() / I3Units::m) * CLHEP::m;
    G4double y = (particle.GetY() / I3Units::m) * CLHEP::m;
    G4double z = (particle.GetZ() / I3Units::m) * CLHEP::m;

    G4ThreeVector direction(particle.GetDir().GetX(),
                            particle.GetDir().GetY(),
                            particle.GetDir().GetZ());

    G4double time = particle.GetTime() / I3Units::ns * CLHEP::ns;

    G4PrimaryVertex * vertex = new G4PrimaryVertex(x, y, z, time);
    G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition * particleDefinition = GetParticleDefinition(particle);

    G4double mass = particleDefinition->GetPDGMass();
    G4double energy = particle.GetEnergy() / I3Units::MeV * CLHEP::MeV;

    G4PrimaryParticle * g4_particle = new G4PrimaryParticle(particleDefinition);

    g4_particle->SetMass(mass);
    g4_particle->SetKineticEnergy(energy);
    g4_particle->SetMomentumDirection(direction);

    if(particle.GetType() == I3Particle::Na22Nucleus) {
        g4_particle->SetCharge(0.);
    }

    vertex->SetPrimary(g4_particle);

    log_trace("Injecting %s: x=%.2f m, y=%.2f m, z=%.2f m, E=%.3f MeV",
              particle.GetTypeString().c_str(),
              x / CLHEP::m,
              y / CLHEP::m,
              z / CLHEP::m,
              energy / CLHEP::MeV);

    return vertex;
}

G4ParticleDefinition * GetParticleDefinition(I3Particle const & particle) {
    G4ParticleTable * particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDef = NULL;

    bool sodium_run = false;

    switch(particle.GetType()) {
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
      return nullptr;
    }

    // special logic for adding sodium particle
    if(sodium_run) {
        G4int Z = 11, A = 22;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;

        particleDef = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    }

    if(!particleDef) {
        log_warn("You passed NULL particleDef \"%s\" ?!", particle.GetTypeString().c_str());
        return nullptr;
    }

    return particleDef;
}

