#include "g4-larsim/g4classes/G4CCMActionInitialization.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"
#include "g4-larsim/g4classes/G4CCMPrimaryGeneratorAction.h"

G4CCMActionInitialization::G4CCMActionInitialization(G4CCMParticleList * particle_list)
    : G4VUserActionInitialization(), particle_list_(particle_list)
{}

void G4CCMActionInitialization::BuildForMaster() const {
}

void G4CCMActionInitialization::Build() const {
    SetUserAction(new G4CCMPrimaryGeneratorAction(particle_list_));
}

