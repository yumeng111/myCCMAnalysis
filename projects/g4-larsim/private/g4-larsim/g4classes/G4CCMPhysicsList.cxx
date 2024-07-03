//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file G4CCMDecayPhysics.cc
/// \brief Implementation of the G4CCMDecayPhysics class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "g4-larsim/g4classes/G4CCMPhysicsList.h"

#include <G4Version.hh>
#include <G4UnitsTable.hh>
#include <G4EmStandardPhysics_option4.hh>
#include <G4EmStandardPhysics_option3.hh>
#include <G4ParticleTypes.hh>
#include <G4IonConstructor.hh>
#include <G4PhysicsListHelper.hh>
#include <G4Radioactivation.hh>
#include <G4SystemOfUnits.hh>
#include <G4NuclideTable.hh>
#include <G4LossTableManager.hh>
#include <G4UAtomicDeexcitation.hh>
#include <G4NuclideTable.hh>
#include <G4NuclearLevelData.hh>
#include <G4DeexPrecoParameters.hh>
#include <G4PhysListUtil.hh>
#include <G4EmBuilder.hh>
#include <globals.hh>
#include <G4DecayPhysics.hh>
#include <G4RadioactiveDecayPhysics.hh>
#include <G4OpticalPhysics.hh>
#include <G4OpticalParameters.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMPhysicsList::G4CCMPhysicsList(G4int ver):  G4VModularPhysicsList()
{
    defaultCutValue = 0.1*CLHEP::mm;
    SetVerboseLevel(ver);

	// EM physics
	emPhysicsList = new G4EmStandardPhysics_option4(ver);

    // Decays 
    decayPhysicsList = new G4DecayPhysics(ver);

    // Optical Physics
    opticalPhysicsList = new G4OpticalPhysics(ver);

#if G4VERSION_NUMBER >= 1100
    //RegisterPhysics(opticalPhysicsList);
    G4OpticalParameters* params = G4OpticalParameters::Instance();
#elif G4VERSION_NUMBER >= 1000
    G4OpticalPhysics* params = opticalPhysicsList;
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
    //RegisterPhysics(opticalPhysicsList);
#endif

    // instantiate Physics List infrastructure 
    //
    G4PhysListUtil::InitialiseParameters();

    // update G4NuclideTable time limit
    //
    const G4double meanLife = 1*picosecond;  
    G4NuclideTable::GetInstance()->SetMeanLifeThreshold(meanLife);  
    G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);

    // define flags for the atomic de-excitation module
    //
    G4EmParameters::Instance()->SetDefaults();
    G4EmParameters::Instance()->SetAugerCascade(true);
    G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true);    

    // define flags for nuclear gamma de-excitation model
    //
    G4DeexPrecoParameters* deex = 
    G4NuclearLevelData::GetInstance()->GetParameters();
    deex->SetCorrelatedGamma(false);
    deex->SetStoreAllLevels(true);
    deex->SetInternalConversionFlag(true);	  
    deex->SetIsomerProduction(true);  
    deex->SetMaxLifeTime(meanLife);

}

G4CCMPhysicsList::~G4CCMPhysicsList() {  
	
    delete emPhysicsList;
    delete opticalPhysicsList;
    delete decayPhysicsList;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMPhysicsList::ConstructParticle() {
    // minimal set of particles for EM physics and radioactive decay
    G4EmBuilder::ConstructMinimalEmSet();
    
    // Register optical photon
    G4OpticalPhoton::Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMPhysicsList::ConstructProcess() {
    AddTransportation();

    // electromagnetic physics list
    emPhysicsList->ConstructProcess();

    // decay physics list
    decayPhysicsList->ConstructProcess();

    // optical physics list
    opticalPhysicsList->ConstructProcess();

    // radioactive decay things
    G4Radioactivation* radioactiveDecay = new G4Radioactivation("Radioactivation", 1.0e+60*CLHEP::year);

    G4bool ARMflag = false;
    radioactiveDecay->SetARM(ARMflag);        //Atomic Rearangement
      
    // EM physics constructor is not used in this example, so
    // it is needed to instantiate and to initialize atomic deexcitation
    //
    G4LossTableManager* man = G4LossTableManager::Instance();
    G4VAtomDeexcitation* deex = man->AtomDeexcitation();
    if (nullptr == deex) {
     deex = new G4UAtomicDeexcitation();
     man->SetAtomDeexcitation(deex);
    }
    deex->InitialiseAtomicDeexcitation();

    // register radioactiveDecay
    //
    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

    // quick aside, printing decay table for sodium22
    G4int Z = 11; // Atomic number for Sodium
    G4int A = 22; // Mass number for Sodium
    G4double E = 0.0 * keV; // Excitation energy
    G4ParticleDefinition* sodium22 = G4IonTable::GetIonTable()->GetIon(Z, A, E);

    G4DecayTable* decayTable = radioactiveDecay->GetDecayTable(sodium22);

    G4cout << "Decay table for " << sodium22->GetParticleName() << ":" << G4endl;
    for (G4int i = 0; i < decayTable->entries(); ++i) {
        G4VDecayChannel* decayChannel = decayTable->GetDecayChannel(i);
        if (decayChannel) {
            G4cout << "Decay channel " << i + 1 << ":" << G4endl;
            decayChannel->DumpInfo();
        }
    } 

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
