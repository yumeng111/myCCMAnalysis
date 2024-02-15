#include <g4-larsim/g4classes/G4CCMSD.h>

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

#include "G4Version.hh"

#if G4VERSION_NUMBER >= 950
// The G4MaterialPropertyVector is gone since 4.9.5.
// It has been typedef'd to a G4UnorderedPhysicsVector
// with a different interface. Try to support both
// versions with an ifdef.
#define MATERIAL_PROPERTY_VECTOR_IS_PHYSICS_VECTOR
#endif


// i'm assuming scPositions == scintillation positions so kinda like pmt for detecting just scintillation photons ??
// so I will change this to be uncoated pmt positions and coated pmt positions and keep track of both cherenkov and scintillation photons
// i will make LightKey store both scintillation and cherenkov photons as well as total energy deposited and cogTime if that's relevant to ccm
G4CCMSD::G4CCMSD(G4String name, const std::map<LightKey, G4ThreeVector>& coatedPositions, const std::map<LightKey, G4ThreeVector>& uncoatedPositions, G4double orientation, G4bool applyTimeCut) 
  : G4VSensitiveDetector(name), coatedPositions_(coatedPositions), uncoatedPositions_(uncoatedPositions), orientation_(orientation), applyTimeCut_(applyTimeCut)
{}


G4CCMSD::~G4CCMSD() {}


void G4CCMSD::Initialize(G4HCofThisEvent* HCE)
{
  // first coated pmts
  std::map<LightKey, G4ThreeVector>::const_iterator coated_pmt_iter;
  for(coated_pmt_iter=coatedPositions_.begin(); coated_pmt_iter!=coatedPositions_.end(); ++coated_pmt_iter)
  {
    sumEdep_[coated_pmt_iter->first] = 0;
    cogTime_[coated_pmt_iter->first] = 0; // idk why we are keeping track of this

    scintillatorCounter_[coated_pmt_iter->first] = 0;
    cherenkovCounter_[coated_pmt_iter->first] = 0; // let's also add in a cherenkov counter
  }
  
  // now uncoated pmts
  std::map<LightKey, G4ThreeVector>::const_iterator uncoated_pmt_iter;
  for(uncoated_pmt_iter=uncoatedPositions_.begin(); uncoated_pmt_iter!=uncoatedPositions_.end(); ++uncoated_pmt_iter)
  {
    sumEdep_[uncoated_pmt_iter->first] = 0;
    cogTime_[uncoated_pmt_iter->first] = 0; // idk why we are keeping track of this

    scintillatorCounter_[uncoated_pmt_iter->first] = 0;
    cherenkovCounter_[uncoated_pmt_iter->first] = 0; // let's also add in a cherenkov counter
  }
  
}


G4bool G4CCMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit(); // what is aStep??? saw some comments calling it angular step??? does not explain anything 
  G4double time = aStep->GetPreStepPoint()->GetGlobalTime();

  G4double scintillatorNumber = GetScintillationPhotons(aStep,edep);
  if (edep<=0 && scintillatorNumber<=0) return false;
  if ( applyTimeCut_ && (time > 1.0*CLHEP::ms) ) {
    return false;  
  }
  
  // analogous function for cherenkov photons
  G4double cherenkovNumber = GetCherenkovPhotons(aStep,edep);
  if (edep<=0 && cherenkovNumber<=0) return false;
  if ( applyTimeCut_ && (time > 1.0*CLHEP::ms) ) {
    return false;  
  }
  
  // now iterate over coated pmts
  std::map<LightKey, G4ThreeVector>::const_iterator coated_pmt_iter;
  for(coated_pmt_iter=coatedPositions_.begin(); coated_pmt_iter!=coatedPositions_.end(); ++coated_pmt_iter)
  {

    G4ThreeVector delta = aStep->GetPreStepPoint()->GetPosition() - coated_pmt_iter->second;
    delta.rotate(orientation_ ,G4ThreeVector(0,0,1) ); // idk why the direction is hard coded in....
    
    
    // seems not necessary for ccm
    /*
    G4double delta_X = delta.getX();

    //see DOI: 10.5445/IR/1000131545
    G4double x0 = -8.08952046e+02;
    G4double x1 =  8.11856108e+02;
    G4double con=  1.01318943+00;
    G4double k1 =  5.07342487e-05;
    G4double k2 = -1.68039952e-05;
    G4double k3 = -1.71503981e-05;


    G4double correction_factor_light = 1.0;

    if (delta_X < x0)
    {
        correction_factor_light =  k1 * delta_X + con + k2 * (delta_X - x0) * (delta_X - x0);
    }
    else if (x0 <= delta_X && delta_X <= x1)
    {
        correction_factor_light = k1 * delta_X + con;
    }
    else
    {
        correction_factor_light =  k1 * delta_X + con + k3 * (delta_X - x1) * (delta_X - x1);
    }

    G4double offset_factor_time = -0.00611  * delta_X + 15.73;*/

    sumEdep_[coated_pmt_iter->first] += edep;
    cogTime_[coated_pmt_iter->first] += edep*(time + offset_factor_time);
    scintillatorCounter_[coated_pmt_iter->first] += scintillatorNumber; //* correction_factor_light;
    cherenkovCounter_[coated_pmt_iter->first] += cherenkovNumber;
  }
  
  // now iterate over uncoated pmts
  std::map<LightKey, G4ThreeVector>::const_iterator uncoated_pmt_iter;
  for(uncoated_pmt_iter=uncoatedPositions_.begin(); uncoated_pmt_iter!=uncoatedPositions_.end(); ++uncoated_pmt_iter)
  {

    G4ThreeVector delta = aStep->GetPreStepPoint()->GetPosition() - uncoated_pmt_iter->second;
    delta.rotate(orientation_ ,G4ThreeVector(0,0,1) ); // idk why the direction is hard coded in....

    sumEdep_[uncoated_pmt_iter->first] += edep;
    cogTime_[uncoated_pmt_iter->first] += edep*(time + offset_factor_time);
    scintillatorCounter_[uncoated_pmt_iter->first] += scintillatorNumber; //* correction_factor_light;
    cherenkovCounter_[uncoated_pmt_iter->first] += cherenkovNumber;
  }

  return true;
}


void G4CCMSD::EndOfEvent(G4HCofThisEvent*)
{
  // end event for coated pmts
  std::map<LightKey, G4ThreeVector>::const_iterator coated_pmt_iter;
  for(coated_pmt_iter=coatedPositions_.begin(); coated_pmt_iter!=coatedPositions_.end(); ++coated_pmt_iter)
  {
    if(sumEdep_[coated_pmt_iter->first]>0)
    {
      cogTime_[coated_pmt_iter->first] /= sumEdep_[coated_pmt_iter->first]; // cogTime just time ??? 
    }
  }
  // end event for uncoated pmts
  std::map<LightKey, G4ThreeVector>::const_iterator uncoated_pmt_iter;
  for(uncoated_pmt_iter=uncoatedPositions_.begin(); uncoated_pmt_iter!=uncoatedPositions_.end(); ++uncoated_pmt_iter)
  {
    if(sumEdep_[uncoated_pmt_iter->first]>0)
    {
      cogTime_[uncoated_pmt_iter->first] /= sumEdep_[uncoated_pmt_iter->first]; // cogTime just time ??? 
    }
  }
}


G4double G4CCMSD::GetScintillationPhotons(G4Step* aStep, G4double deposit_energy)
{
  G4Track* aTrack = aStep->GetTrack();
  const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
  const G4Material* aMaterial = aTrack->GetMaterial();

  G4StepPoint* pPreStepPoint  = aStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();
  G4TouchableHandle touch1 = pPreStepPoint->GetTouchableHandle();
  G4TouchableHandle touch2 = pPostStepPoint->GetTouchableHandle();
  G4VPhysicalVolume* volume1 = touch1->GetVolume();
  G4VPhysicalVolume* volume2 = touch2->GetVolume();
  G4String name1 = volume1->GetName();
  G4String name2 = volume2->GetName();

  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable) {
    return 0.0;
  }

  G4double scintyield = aMaterialPropertiesTable->GetConstProperty("SCINTILLATIONYIELD");

  if (!scintyield) {
    return 0.0;
  }

  G4double birks = aMaterial->GetIonisation()->GetBirksConstant(); // birks const.
  G4double step_length = aStep->GetStepLength()/CLHEP::mm;
  G4double number_photons;
  if (step_length==0) {
    number_photons = 0;
  }
  else {
    number_photons = scintyield*deposit_energy/(1.+birks*(deposit_energy/step_length));
  }

  return number_photons;

}


G4double G4CCMSD::GetCherenkovPhotons(G4Step* aStep, G4double deposit_energy)
{
  G4Track* aTrack = aStep->GetTrack();
  const G4DynamicParticle* aParticle = aTrack->GetDynamicParticle();
  const G4Material* aMaterial = aTrack->GetMaterial();

  G4StepPoint* pPreStepPoint  = aStep->GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep->GetPostStepPoint();
  G4TouchableHandle touch1 = pPreStepPoint->GetTouchableHandle();
  G4TouchableHandle touch2 = pPostStepPoint->GetTouchableHandle();
  G4VPhysicalVolume* volume1 = touch1->GetVolume();
  G4VPhysicalVolume* volume2 = touch2->GetVolume();
  G4String name1 = volume1->GetName();
  G4String name2 = volume2->GetName();

  G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable) {
    return 0.0;
  }

  // sorry if "chkv" as abbreviation for "cherenkov" is annoying...can be changed
  G4double chkvyield = aMaterialPropertiesTable->GetConstProperty("CHERENKOVYIELD");

  if (!chkvyield) {
    return 0.0;
  }

  // do we need birks for cherenkov ??? who is to say
  G4double birks = aMaterial->GetIonisation()->GetBirksConstant(); // birks const.
  G4double step_length = aStep->GetStepLength()/CLHEP::mm;
  G4double number_photons;
  if (step_length==0) {
    number_photons = 0;
  }
  else {
    number_photons = chkvyield*deposit_energy/(1.+birks*(deposit_energy/step_length));
  }

  return number_photons;

}


