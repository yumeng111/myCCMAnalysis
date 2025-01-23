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
//
/// \file optical/LXe/src/LXePrimaryGeneratorAction.cc
/// \brief Implementation of the LXePrimaryGeneratorAction class
//
//
#include "G4CCMVisPrimaryGeneratorAction.h"

#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMVisPrimaryGeneratorAction::G4CCMVisPrimaryGeneratorAction(){
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("geantino");

    G4double inset = 0.1 * cm;
    G4ThreeVector pos(0.,0.,0. + inset);
    G4ThreeVector mom(0.,0.,-1.);

    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(mom);
    fParticleGun->SetParticleMomentum(0.*GeV);
    fParticleGun->SetParticleDefinition(particle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMVisPrimaryGeneratorAction::~G4CCMVisPrimaryGeneratorAction() { 
    delete fParticleGun; 
    //delete fGeneralParticleSource;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMVisPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
        G4int Z = 11, A = 22;
        G4double ionCharge   = 0.*eplus;
        G4double excitEnergy = 0.*keV;

        G4ParticleDefinition* ion
           = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticleCharge(ionCharge);
    }    
    //create vertex
    //   
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

