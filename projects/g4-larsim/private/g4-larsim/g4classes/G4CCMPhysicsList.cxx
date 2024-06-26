
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
#include <G4Ions.hh>
#include <G4ParticleDefinition.hh>

G4CCMPhysicsList::G4CCMPhysicsList(G4int ver)  : G4VUserPhysicsList() {

    defaultCutValue = 0.7*CLHEP::mm;
    SetVerboseLevel(ver);

    // EM Physics
    RegisterPhysics( new G4EmStandardPhysics_option4(ver));

    // Decays
    RegisterPhysics( new G4DecayPhysics(ver) );
    RegisterPhysics( new G4IonPhysics(ver) );
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
    
//    // Synchroton Radiation & GN Physics
//    RegisterPhysics( new G4EmExtraPhysics(ver) );
//
//    // Hadron Elastic scattering
//    RegisterPhysics( new G4HadronElasticPhysicsHP(ver) );
//
//    // Hadron Physics
//    RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP(ver) );
//
//    // Stopping Physics
//    RegisterPhysics( new G4StoppingPhysics(ver) );

}

G4CCMPhysicsList::~G4CCMPhysicsList() {
}

void G4CCMPhysicsList::AddSodiumDecay() {

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4IonTable* ionTable = particleTable->GetIonTable();
    
    // Define properties of Sodium-22
    G4double sodium22Mass = 22.0 * CLHEP::amu; // Mass of Sodium-22
    G4int sodium22Charge = 11; // Charge of Sodium-22 (atomic number)
    G4int sodium22Spin = 3; // Spin of Sodium-22 
    G4int sodium22Parity = 0; // Parity of Sodium-22 
    G4int sodium22Isospin = 1; // Isospin of Sodium-22...idk check 
    G4int sodium22IsospinZ = 1; // Isospin Z of Sodium-22...idk check 
    G4String sodium22ParticleType = "ion"; // Particle type of Sodium-22
    G4int sodium22Lepton = 0; // Lepton number of Sodium-22 (arbitrary value)
    G4int sodium22Baryon = 1; // Baryon number of Sodium-22 (arbitrary value)
    G4int sodium22Encoding = 1000110220; // PDG encoding of Sodium-22 (arbitrary value)
    G4bool sodium22Stable = false; // Stability of Sodium-22 (arbitrary value)
    G4double sodium22Lifetime = 2.6019 * 365.25 * 24 * 60 * 60 * s; // Lifetime of Sodium-22
    G4DecayTable* sodium22DecayTable = nullptr; // Decay table of Sodium-22 (null for stable particles)
    
    // Create Sodium-22 particle definition
    G4ParticleDefinition* sodium22Particle = new G4ParticleDefinition("sodium22", sodium22Mass, 0.0, sodium22Charge, sodium22Spin, sodium22Parity,
                                                                      0, sodium22Isospin, sodium22IsospinZ, 0, sodium22ParticleType,
                                                                      sodium22Lepton, sodium22Baryon, sodium22Encoding, sodium22Stable,
                                                                      sodium22Lifetime, sodium22DecayTable);
    
    ionTable->Insert(sodium22Particle);
    
    // Add beta plus decay for sodium-22
    G4ParticleDefinition* sodium22 = ionTable->GetIon(11, 22, 3); 
    G4double branchingRatio = 1.0;  // Branching ratio for this decay mode
    G4double endpointEnergy = 2.842 * MeV;  // Endpoint energy for the beta-plus decay of Sodium-22
    G4double daughterExcitation = 1.275 * MeV;  // Excitation energy of the daughter nucleus
    G4Ions::G4FloatLevelBase flb = G4Ions::G4FloatLevelBase::no_Float;  // Float level base
    G4BetaDecayType decayType = G4BetaDecayType::allowed;  // Decay type
    G4BetaPlusDecay* betaPlusDecay = new G4BetaPlusDecay(sodium22, branchingRatio, endpointEnergy, daughterExcitation, flb, decayType);
    G4DecayTable* decayTable = new G4DecayTable();
    decayTable->Insert(betaPlusDecay); // Insert the beta plus decay process into the decay table
    ionTable->GetIon(11, 22, 3)->SetDecayTable(decayTable);

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
