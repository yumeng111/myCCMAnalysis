#include <g4-larsim/g4classes/G4CCMUserTrackingAction.h>

#include "G4Track.hh"
#include "G4UserLimits.hh"
#include "G4TrackVector.hh"
#include "G4TrackingManager.hh"

G4CCMUserTrackingAction::G4CCMUserTrackingAction(){}

void G4CCMUserTrackingAction::PreUserTrackingAction(const G4Track*){}

void G4CCMUserTrackingAction::PostUserTrackingAction(const G4Track* track)
{
  const G4LogicalVolume *volume = track->GetLogicalVolumeAtVertex();
  G4UserLimits *limit = volume->GetUserLimits();
  if(!limit) G4cout << "----> G4LogicalVolume: " << volume->GetName() << " has no defined G4UserLimit" << G4endl;
  G4double threshold = limit->GetUserMinEkine(*track);
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      {
        //check if secondary particle is a gamma
        if((*secondaries)[i]->GetDefinition()->GetParticleName() == "gamma")
        {
          //check if particle energy is below threshold; if true, kill the particle
          G4double energy = (*secondaries)[i]->GetTotalEnergy();
          if(energy < threshold) (*secondaries)[i]->SetTrackStatus(fStopAndKill);
        }
      }
    }
  }
}
