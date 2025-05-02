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
/// \file G4CCMPhysicsList.hh
/// \brief Definition of the G4CCMPhysicsList class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4CCMPhysicsList_h
#define G4CCMPhysicsList_h 1

#include <G4VUserPhysicsList.hh>
#include "G4VModularPhysicsList.hh"
#include "G4EmConfigurator.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//class G4CCMPhysicsList: public G4VUserPhysicsList {
class G4CCMPhysicsList: public G4VModularPhysicsList {
public:
    G4CCMPhysicsList(G4int ver=1);
    ~G4CCMPhysicsList();
    double PhotonSamplingFactor_ = 1.0;
    bool SimulateNuclearRecoils_ = false;
    double G4RangeCut_ = 1e-6; // 1 um
    double G4EDepMin_ = 1e-8; // 0.01 keV

protected:
    // Construct particle and physics
    void ConstructParticle() override;
    void ConstructProcess()  override;
    void SetCuts() override;
private:
    G4VPhysicsConstructor* emPhysicsList;
    G4VPhysicsConstructor* emExtraPhysicsList;
    G4VPhysicsConstructor* decayPhysicsList;
    //G4VPhysicsConstructor* radioactiveDecayPhysicsList;
    G4VPhysicsConstructor* hadronElasticPhysicsList;
    G4VPhysicsConstructor* hadronPhysicsList;
    G4VPhysicsConstructor* stoppingPhysicsList;
    G4VPhysicsConstructor* ionPhysicsList;
    G4VPhysicsConstructor* opticalPhysicsList;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif



