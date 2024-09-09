#ifndef G4CCMActionInitialization_h
#define G4CCMActionInitialization_h 1

#include "g4-larsim/g4classes/G4CCMPrimaryGeneratorAction.h"

#include "G4VUserActionInitialization.hh"

class G4CCMActionInitialization : public G4VUserActionInitialization {
public:
    G4CCMActionInitialization(G4CCMParticleList * particle_list);
    ~G4CCMActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;

private:
    G4CCMParticleList * particle_list_;
};

#endif

