#ifndef G4CCMPrimaryGeneratorAction_h
#define G4CCMPrimaryGeneratorAction_h

#include <mutex>
#include <vector>

#include "dataclasses/physics/I3Particle.h"

#include "G4PrimaryVertex.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;

class G4CCMParticleList {
    std::mutex mutex;
    std::vector<I3Particle> particle_list_;
    unsigned int current_particle_;
public:
    G4CCMParticleList();
    ~G4CCMParticleList() = default;

    void AddParticle(I3Particle p);
    void AddParticles(std::vector<I3Particle> p);
    void SetParticles(std::vector<I3Particle> p);
    void Clear();
    I3Particle GetNextParticle();
};

class G4CCMPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    G4CCMParticleList * particle_list_;
public:
    G4CCMPrimaryGeneratorAction(G4CCMParticleList * particle_list);
    ~G4CCMPrimaryGeneratorAction() = default;

    void GeneratePrimaries(G4Event* anEvent) override;
};

G4PrimaryVertex * CreatePrimaryVertex(I3Particle const & primary);
G4ParticleDefinition * GetParticleDefinition(I3Particle const & primary);

#endif // G4CCMPrimaryGeneratorAction_h

