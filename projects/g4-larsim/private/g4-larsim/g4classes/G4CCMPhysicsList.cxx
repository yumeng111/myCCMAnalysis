#include <globals.hh>
#include <G4Version.hh>

#include "G4CCMPhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4OpticalPhysics.hh"

#include <G4ProcessManager.hh>
#include <G4ParticleTypes.hh>
#include <G4UserSpecialCuts.hh>


G4CCMPhysicsList::G4CCMPhysicsList(G4int ver)
    : G4VUserPhysicsList() {

    defaultCutValue = 0.7*CLHEP::mm;
    SetVerboseLevel(ver);

    // EM Physics
    RegisterPhysics( new G4EmStandardPhysics_option4(ver));

    // Synchroton Radiation & GN Physics
    RegisterPhysics( new G4EmExtraPhysics(ver) );

    // Decays
    RegisterPhysics( new G4DecayPhysics(ver) );
    RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );

    // Hadron Elastic scattering
    RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );

    // Hadron Physics
    RegisterPhysics( new G4HadronPhysicsFTFP_BERT_HP(ver) );

    // Stopping Physics
    RegisterPhysics( new G4StoppingPhysics(ver) );

    // Ion Physics
    RegisterPhysics( new G4IonPhysics(ver) );

    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics(ver);

    //The Following lines set more specific parameters for scintillation.
    opticalPhysics->SetWLSTimeProfile("delta");
    opticalPhysics->SetScintillationYieldFactor(1.0);//for e/m
    //opticalPhysics->SetScintillationYieldFactor(.25);//for nucleon
    opticalPhysics->SetScintillationExcitationRatio(0.0);

    opticalPhysics->SetMaxNumPhotonsPerStep(100);
    opticalPhysics->SetMaxBetaChangePerStep(10.0);

    opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
    opticalPhysics->SetTrackSecondariesFirst(kScintillation, true);

    RegisterPhysics(opticalPhysics);
}

G4CCMPhysicsList::~G4CCMPhysicsList() {
}


void G4CCMPhysicsList::SetCuts() {
    if (verboseLevel >1){
        G4cout << "G4CCMPhysicsList::SetCuts:";
    }
    //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
    //   the default cut value for all particle types

    SetCutsWithDefault();

    //Set proton cut value to 0 for producing low energy recoil nucleus
    SetCutValue(0.0, "proton");
}
