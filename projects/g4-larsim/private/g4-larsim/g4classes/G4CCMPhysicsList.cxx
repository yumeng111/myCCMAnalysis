
#include "g4-larsim/g4classes/G4CCMPhysicsList.h"

#include <globals.hh>
#include <G4Version.hh>
#include <G4IonPhysics.hh>
#include <G4DecayPhysics.hh>
#include <G4ParticleTypes.hh>
#include <G4OpticalPhysics.hh>
#include <G4EmExtraPhysics.hh>
#include <G4ProcessManager.hh>
#include <G4UserSpecialCuts.hh>
#include <G4StoppingPhysics.hh>
#include <G4HadronElasticPhysicsHP.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4HadronPhysicsQGSP_BERT_HP.hh>
#include <G4EmStandardPhysics_option4.hh>

#include <FTFP_BERT.hh>
#include <G4IonTable.hh>
#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4OpticalParameters.hh>
#include <G4VModularPhysicsList.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4DecayTable.hh>
#include <G4RadioactiveDecay.hh>
#include <G4BetaPlusDecay.hh>
#include <G4GenericIon.hh>
#include <G4PhysicsListHelper.hh>
#include <G4WrapperProcess.hh>

G4CCMPhysicsList::G4CCMPhysicsList(G4int ver)  : G4VUserPhysicsList() {

    // EM Physics
    RegisterPhysics( new G4EmStandardPhysics_option4(ver));

    // Decays
    RegisterPhysics( new G4DecayPhysics(ver) );
    RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );

    // Add beta plus decay
    //AddBetaPlusDecay();

    // Ion Physics
    RegisterPhysics( new G4IonPhysics(ver) );

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
    
//    defaultCutValue = 0.7*CLHEP::mm;
//    SetVerboseLevel(ver);
//
//    // EM Physics
//    RegisterPhysics( new G4EmStandardPhysics_option4(ver));
//
//    // Synchroton Radiation & GN Physics
//    RegisterPhysics( new G4EmExtraPhysics(ver) );
//
//    // Decays
//    RegisterPhysics( new G4DecayPhysics(ver) );
//    RegisterPhysics( new G4RadioactiveDecayPhysics(ver) );
//
//    // Hadron Elastic scattering
//    RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );
//
//    // Hadron Physics
//    RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(ver) );
//
//    // Stopping Physics
//    RegisterPhysics( new G4StoppingPhysics(ver) );
//
//    // Ion Physics
//    RegisterPhysics( new G4IonPhysics(ver) );
//
//    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics(ver);
//
//#if G4VERSION_NUMBER >= 1100
//    RegisterPhysics(opticalPhysics);
//    G4OpticalParameters* params = G4OpticalParameters::Instance();
//#elif G4VERSION_NUMBER >= 1000
//    G4OpticalPhysics* params = opticalPhysics;
//#endif
//
//    //The Following lines set more specific parameters for scintillation.
//    params->SetWLSTimeProfile("delta");
//#if G4VERSION_NUMBER >= 1100
//#else
//    // In geant4.11.1 and beyond, this needs to be set in the material properties table
//    // These can also be set on a per-particle basis and a per-component basis
//    params->SetScintillationYieldFactor(1.0);//for e/m
//    //params->SetScintillationYieldFactor(.25);//for nucleon
//    params->SetScintillationExcitationRatio(0.0);
//#endif
//
//#if G4VERSION_NUMBER >= 1100
//    params->SetCerenkovMaxPhotonsPerStep(100);
//    params->SetCerenkovMaxBetaChange(10.0);
//
//    params->SetCerenkovTrackSecondariesFirst(true);
//    params->SetScintTrackSecondariesFirst(true);
//#else
//    params->SetMaxNumPhotonsPerStep(100);
//    params->SetMaxBetaChangePerStep(10.0);
//
//    params->SetTrackSecondariesFirst(kCerenkov, true);
//    params->SetTrackSecondariesFirst(kScintillation, true);
//#endif
//
//#if G4VERSION_NUMBER >= 1100
//#elif G4VERSION_NUMBER >= 1000
//    RegisterPhysics(opticalPhysics);
//#endif
}

G4CCMPhysicsList::~G4CCMPhysicsList() {
}

void G4CCMPhysicsList::AddBetaPlusDecay() {
    
    // Add beta plus decay for sodium-22
    G4double branchingRatio = 1.0;  // Branching ratio for this decay mode
    G4double endpointEnergy = 2.842 * MeV;  // Endpoint energy for the beta-plus decay of Sodium-22
    G4double daughterExcitation = 0.0;  // Excitation energy of the daughter nucleus
    G4Ions::G4FloatLevelBase flb = G4Ions::G4FloatLevelBase::no_Float;  // Float level base
    G4BetaDecayType decayType = G4BetaDecayType::allowed;  // Decay type

    // Get the definition of Sodium-22
    G4ParticleDefinition* sodium22Ion = G4IonTable::GetIonTable()->GetIon(11, 22, 0); // Z=11, A=22

    // Create the beta plus decay process for sodium-22
    G4BetaPlusDecay* betaPlusDecay = new G4BetaPlusDecay(sodium22Ion, branchingRatio, endpointEnergy, daughterExcitation, flb, decayType);
    //G4WrapperProcess* wrapperProcess = new G4WrapperProcess("WrapperForBetaPlusDecay");
    //wrapperProcess->RegisterProcess(betaPlusDecay);

    //G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

    //auto particleIterator = GetParticleIterator();
    //particleIterator->reset();
    //while ((*particleIterator)()) {
    //    G4ParticleDefinition* particle = particleIterator->value();
    //    //if (betaPlusDecay->IsApplicable(*particle) && !particle->IsShortLived()) {
    //    ph->RegisterProcess(wrapperProcess, particle);
    //    //}
    //}
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
