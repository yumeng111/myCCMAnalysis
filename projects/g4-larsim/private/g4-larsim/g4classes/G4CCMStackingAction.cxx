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
/// \file optical/LXe/src/LXeStackingAction.cc
/// \brief Implementation of the LXeStackingAction class
//
//
#include "g4-larsim/g4classes/G4CCMStackingAction.h"

#include "g4-larsim/g4classes/G4CCMEventAction.h"

#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CCMStackingAction::G4CCMStackingAction(G4CCMEventAction* ea)
  : fEventAction(ea)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack G4CCMStackingAction::ClassifyNewTrack(
  const G4Track* aTrack)
{
  // Count what process generated the optical photons
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    // particle is optical photon
    if(aTrack->GetParentID() > 0)
    {
      // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation")
        fEventAction->IncPhotonCount_Scint();
      else if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov")
        fEventAction->IncPhotonCount_Ceren();
    }
  }
  return fUrgent;
}

