#ifndef G4CCMPhysicsList_H
#define G4CCMPhysicsList_H

#include <G4VModularPhysicsList.hh>

/**
 * Implementation of G4VModularPhysicsList. The top-level physics list. It combines all the other physics lists:
 *
 *  \li G4CCMGeneralPhysics
 *  \li G4CCMHadronPhysics
 *  \li G4CCMIonPhysics
 *  \li G4CCMMuonPhysics
 *  \li G4CCMEMPhysics
 */
class G4CCMPhysicsList: public G4VModularPhysicsList {
public:
    G4CCMPhysicsList(G4int ver=1);
    ~G4CCMPhysicsList();
private:
    void SetCuts();
    void AddBetaPlusDecay();
};

#endif // g4-larsim_G4CCMPhysicsList_H
