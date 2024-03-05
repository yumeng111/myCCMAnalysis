#ifndef G4CCMRUNMANAGER_H
#define G4CCMRUNMANAGER_H

#include <G4RunManager.hh>

class G4ParticleGun;

/**
 * Implementation of G4RunManager
 */
class G4CCMRunManager: public G4RunManager
{
 public:
  G4CCMRunManager();

  static G4CCMRunManager* GetCCMRunManager() {return (G4CCMRunManager*)GetRunManager();}

  void InitializeRun();
  void InjectParticle(G4ParticleGun* particleGun);
  void TerminateRun();

 protected:
  G4Event* GenerateEvent(G4int i_event);

 private:
  // This method is an exact copy of UpdateScoring which is private in the G4RunManager
  void Update_Scoring();

};

#endif
