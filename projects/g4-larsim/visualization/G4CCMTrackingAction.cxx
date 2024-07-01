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
/// \file optical/LXe/src/LXeTrackingAction.cc
/// \brief Implementation of the LXeTrackingAction class
//
//
#include "G4CCMTrackingAction.h"
#include "g4-larsim/g4classes/G4CCMDetectorConstruction.h"
#include "G4CCMTrajectory.h"
#include "G4CCMUserTrackInformation.h"

#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include <G4ParticleTypes.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMTrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
    // Let this be up to the user via vis.mac
    //  fpTrackingManager->SetStoreTrajectory(true);

    // Use custom trajectory class
    fpTrackingManager->SetTrajectory(new G4CCMTrajectory(aTrack));

    // This user track information is only relevant to the photons
    fpTrackingManager->SetUserTrackInformation(new G4CCMUserTrackInformation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4CCMTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
    auto trajectory = (G4CCMTrajectory*) fpTrackingManager->GimmeTrajectory();
    auto trackInformation = (G4CCMUserTrackInformation*) aTrack->GetUserInformation();

    if ((trackInformation->GetTrackStatus()) == 16){ // primary lar deposition
        trajectory->LArScint();
        trajectory->SetForceDrawTrajectory(true);
    }
    // Let's choose to draw only the photons that hit the sphere and a pmt
    if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
        const G4VProcess* creator = aTrack->GetCreatorProcess();
        if(trackInformation->GetTrackStatus() & hitPMT)
            trajectory->SetDrawTrajectory(true);
    }
    // draw all other (not optical photon) trajectories
    else
        trajectory->SetDrawTrajectory(true);

    if(trackInformation->GetForceDrawTrajectory())
        trajectory->SetDrawTrajectory(true);
    

}

