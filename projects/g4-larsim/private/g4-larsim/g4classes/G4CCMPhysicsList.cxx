
#include "g4-larsim/g4classes/G4CCMPhysicsList.h"

#include <G4VUserPhysicsList.hh>
#include <globals.hh>
#include <G4Version.hh>
#include <G4IonPhysics.hh>
#include <G4DecayPhysics.hh>
#include <G4OpticalPhysics.hh>
#include <G4StoppingPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4Radioactivation.hh>
#include <G4EmStandardPhysics_option4.hh>

#include <G4SystemOfUnits.hh>
#include <G4OpticalParameters.hh>
#include <G4GenericIon.hh>
#include <G4PhysicsListHelper.hh>

#include <G4UnitsTable.hh>
#include <G4IonConstructor.hh>
#include <G4LossTableManager.hh>
#include <G4UAtomicDeexcitation.hh>
#include <G4NuclideTable.hh>
#include <G4NuclearLevelData.hh>
#include <G4DeexPrecoParameters.hh>
#include <G4PhysListUtil.hh>
#include <G4EmBuilder.hh>
#include <G4WrapperProcess.hh>

#include <G4EmExtraPhysics.hh>
#include <G4HadronElasticPhysicsHP.hh>
#include <G4HadronPhysicsFTFP_BERT_HP.hh>

#include <G4ProcessManager.hh>
#include <G4ParticleTypes.hh>
#include <G4UserSpecialCuts.hh>


#include <G4EmStandardPhysics.hh>
#include <G4EmParameters.hh>
#include <G4HadronElasticPhysics.hh>
#include <G4HadronPhysicsFTFP_BERT.hh>
#include <G4HadronInelasticQBBC.hh>
#include <G4HadronPhysicsINCLXX.hh>
#include <G4IonElasticPhysics.hh>
#include <G4IonINCLXXPhysics.hh>
#include <FTFP_BERT.hh>

G4CCMPhysicsList::G4CCMPhysicsList(G4int ver)  : G4VUserPhysicsList() {

    defaultCutValue = 0.7*CLHEP::mm;
    SetVerboseLevel(ver);

    // EM Physics
    RegisterPhysics( new G4EmStandardPhysics_option4(ver));

    // Decays
    RegisterPhysics( new G4DecayPhysics(ver) );
    RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );

    // Optical Physics
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics(ver);

#if G4VERSION_NUMBER >= 1100
    RegisterPhysics(opticalPhysics);
    G4OpticalParameters* params = G4OpticalParameters::Instance();
#elif G4VERSION_NUMBER >= 1000
    G4OpticalPhysics* params = opticalPhysics;
#endif

    //The Following lines set more specific parameters for scintillation.
    params->SetWLSTimeProfile("delta");
#if G4VERSION_NUMBER >= 1100
#else
    // In geant4.11.1 and beyond, this needs to be set in the material properties table
    // These can also be set on a per-particle basis and a per-component basis
    params->SetScintillationYieldFactor(1.0);//for e/m
    //params->SetScintillationYieldFactor(.25);//for nucleon
    params->SetScintillationExcitationRatio(0.0);
#endif

#if G4VERSION_NUMBER >= 1100
    params->SetCerenkovMaxPhotonsPerStep(100);
    params->SetCerenkovMaxBetaChange(10.0);

    params->SetCerenkovTrackSecondariesFirst(true);
    params->SetScintTrackSecondariesFirst(true);
#else
    params->SetMaxNumPhotonsPerStep(100);
    params->SetMaxBetaChangePerStep(10.0);

    params->SetTrackSecondariesFirst(kCerenkov, true);
    params->SetTrackSecondariesFirst(kScintillation, true);
#endif

#if G4VERSION_NUMBER >= 1100
#elif G4VERSION_NUMBER >= 1000
    RegisterPhysics(opticalPhysics);
#endif
    std::cout << "loaded optical physics" << std::endl;

    // Synchroton Radiation & GN Physics
    RegisterPhysics( new G4EmExtraPhysics(ver) );

//    // Hadron Elastic scattering
//    RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );
//
//    // Hadron Physics
//    RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(ver) );

    // Stopping Physics
    RegisterPhysics( new G4StoppingPhysics(ver) );
}

G4CCMPhysicsList::~G4CCMPhysicsList() {}

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

