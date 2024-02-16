#include <g4-larsim/g4classes/G4CCMSD.hh>

#include <G4Step.hh>
#include <G4HCofThisEvent.hh>
#include <G4TouchableHistory.hh>
#include <G4ios.hh>
#include <iostream>
#include <ios>
#include <fstream>

#include <G4VProcess.hh>
#include <G4OpticalPhoton.hh>
#include "G4Poisson.hh"

#include <icetray/I3Units.h>
#include <simclasses/CCMMCPE.h>
#include "G4Version.hh"

#if G4VERSION_NUMBER >= 950
// The G4MaterialPropertyVector is gone since 4.9.5.
// It has been typedef'd to a G4UnorderedPhysicsVector
// with a different interface. Try to support both
// versions with an ifdef.
#define MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
#endif

// we G4CCMSD to make one senstive detector aka pmt
// so to make a pmt we need to know : pmt location, pmt orientation (maybe facing dir is good enough for this ?), coating and that should be it
// maybe at some point we need to tell G4 how far the pmts are inserted into the steel array but that's later me's problem
G4CCMSD::G4CCMSD(G4String name, const G4ThreeVector pmtPosition, const G4ThreeVector pmt_facing_dir, const G4bool coatingFlag) 
  : G4VSensitiveDetector(name), pmtPosition_(pmtPosition), pmt_facing_dir_(pmt_facing_dir), coatingFlag_(coatingFlag)
{}


G4CCMSD::~G4CCMSD() {}


void G4CCMSD::Initialize(G4HCofThisEvent* HCE) {
    // idk what we need to initilaize ? 
}


G4bool G4CCMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  
  G4Track *track = aStep->GetTrack();
  if(track->GetDefinition() != G4OpticalPhoton::Definition()) return true; // check to see if this is an optical photon

  // ok so we have a photon! let's save it!
  G4ThreeVector photon_loc = aStep->GetPreStepPoint()->GetPosition(); 
  G4double time = aStep->GetPreStepPoint()->GetGlobalTime(); 
  G4int parent_id = track->GetParentID(); 
  G4ThreeVector momentum_direction = track->GetMomentumDirection(); 

  // now let's save
  // just saving a few things for now for testing
  CCMMCPE this_photon;
  this_photon.npe = 1.0;
  this_photon.time = time;
  photon_info_.push_back(this_photon);  

  // and then extinguish photon!
  track->SetTrackStatus(fStopAndKill);

  return true;
}


void G4CCMSD::EndOfEvent(G4HCofThisEvent*)
{
    // idk how to end event
}

