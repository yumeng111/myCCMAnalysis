#ifndef _G4LARSIM_G4CCMSD_H
#define _G4LARSIM_G4CCMSD_H

#include <simclasses/CCMMCPE.h>
#include <G4VSensitiveDetector.hh>
#include <G4ThreeVector.hh>

#include <map>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/*
 * this sensitive detector represents a single pmt
 *
 * this class keeps track of photon hits on a pmt 
 */
class G4CCMSD : public G4VSensitiveDetector
{
 public:
  G4CCMSD(G4String name, const G4ThreeVector pmtPosition, const G4ThreeVector pmt_facing_dir, const G4bool coatingFlag); 
  ~G4CCMSD();

  /// Methods called by Geant4 framework
  void Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
  void EndOfEvent(G4HCofThisEvent *HCE);

  /// Get relevant photon information for a given pmt 
  std::vector<CCMMCPE> GetPEInfo(){return photon_info_;}

 private:

  const G4ThreeVector pmtPosition_;
  const G4ThreeVector pmt_facing_dir_;
  const G4bool coatingFlag_;

  G4double GetScintillationPhotons(G4Step* aStep, G4double deposit_energy);

  std::vector<CCMMCPE> photon_info_;

};

#endif
