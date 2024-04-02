#ifndef G4CCMRUNMANAGER_H
#define G4CCMRUNMANAGER_H

#include <vector>

#include <G4RunManager.hh>
#include "G4SystemOfUnits.hh"
#include "G4DataVector.hh"

#include <simclasses/CCMMCPE.h>
#include <dataclasses/I3Vector.h>
#include <icetray/CCMPMTKey.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include "g4-larsim/g4classes/G4CCMPMTHit.h"
#include "g4-larsim/g4classes/G4CCMScintHit.h"

class G4ParticleGun;

/**
 * Implementation of G4RunManager
 */
class G4CCMRunManager: public G4RunManager {
    public:
        G4CCMRunManager();

        static G4CCMRunManager* GetCCMRunManager() {return (G4CCMRunManager*)GetRunManager();}

        void InitializeRun();
        void InjectParticle(G4ParticleGun* particleGun);
        void TerminateRun();
        void GetFinalScores(const G4Event* anEvent);
    
        boost::shared_ptr<CCMMCPESeriesMap> GetCCMMCPEMap() { return CCMMCPEMap; }
        boost::shared_ptr<CCMMCPESeries> GetCCMMCPEList(){ return CCMMCPEList; }

    protected:
        G4Event* GenerateEvent(G4int i_event);

    private:
        // This method is an exact copy of UpdateScoring which is private in the G4RunManager
        void Update_Scoring();
        
        G4int fScintCollID = -1;
        G4int fPMTCollID = -1;
        G4int fVerbose = 2;
        G4int fPMTThreshold = 1;

        CCMPMTKey fKey;
        std::vector<CCMMCPE> fCCMMCPEVector;

        G4int fHitCount = 0;
        G4int fPhotonCount_Scint = 0;
        G4int fPhotonCount_Ceren = 0;
        G4int fAbsorptionCount = 0;
        G4int fBoundaryAbsorptionCount = 0;

        G4double fTotE = 0.;

        // These only have meaning if totE > 0
        // If totE = 0 then these won't be set by EndOfEventAction

        G4ThreeVector fEWeightPos;
        G4ThreeVector fReconPos;  // Also relies on hitCount>0
        G4ThreeVector fConvPos;   // true (initial) converstion position
        G4bool fConvPosSet = false;
        G4ThreeVector fPosMax;
        G4double fEdepMax = 0.;

        G4int fPMTsAboveThreshold = 0;

        boost::shared_ptr<CCMMCPESeriesMap> CCMMCPEMap = boost::make_shared<CCMMCPESeriesMap> ();
        boost::shared_ptr<CCMMCPESeries> CCMMCPEList = boost::make_shared<CCMMCPESeries> ();
};

#endif
