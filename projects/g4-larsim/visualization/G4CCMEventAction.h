#ifndef G4CCMEventAction_h
#define G4CCMEventAction_h 1

#include "G4CCMEventMessenger.h"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4UserEventAction.hh"

class G4Event;
class G4CCMDetectorConstruction;

class G4CCMEventAction : public G4UserEventAction
{
    public:

        G4CCMEventAction(const G4CCMDetectorConstruction*);
        ~G4CCMEventAction() override;

        void BeginOfEventAction(const G4Event*) override;
        void EndOfEventAction(const G4Event*) override;

        void SetEventVerbose(G4int v) { fVerbose = v; }

        void SetPMTThreshold(G4int t) { fPMTThreshold = t; }

        void SetForceDrawPhotons(G4bool b) { fForcedrawphotons = b; }
        void SetForceDrawNoPhotons(G4bool b) { fForcenophotons = b; }

        void IncPhotonCount_Scint() { ++fPhotonCount_Scint; }
        void IncPhotonCount_Ceren() { ++fPhotonCount_Ceren; }
        void IncEDep(G4double dep) { fTotE += dep; }
        void IncAbsorption() { ++fAbsorptionCount; }
        void IncBoundaryAbsorption() { ++fBoundaryAbsorptionCount; }
        void IncHitCount(G4int i = 1) { fHitCount += i; }

        void SetEWeightPos(const G4ThreeVector& p) { fEWeightPos = p; }
        void SetReconPos(const G4ThreeVector& p) { fReconPos = p; }
        void SetConvPos(const G4ThreeVector& p)

        {
            fConvPos    = p;
            fConvPosSet = true;
        }
        void SetPosMax(const G4ThreeVector& p, G4double edep)
        {
            fPosMax  = p;
            fEdepMax = edep;
        }

        G4int GetPhotonCount_Scint() const { return fPhotonCount_Scint; }
        G4int GetPhotonCount_Ceren() const { return fPhotonCount_Ceren; }
        G4int GetHitCount() const { return fHitCount; }
        G4double GetEDep() const { return fTotE; }
        G4int GetAbsorptionCount() const { return fAbsorptionCount; }
        G4int GetBoundaryAbsorptionCount() const { return fBoundaryAbsorptionCount; }

        G4ThreeVector GetEWeightPos() { return fEWeightPos; }
        G4ThreeVector GetReconPos() { return fReconPos; }
        G4ThreeVector GetConvPos() { return fConvPos; }
        G4ThreeVector GetPosMax() { return fPosMax; }
        G4double GetEDepMax() { return fEdepMax; }
        G4double IsConvPosSet() { return fConvPosSet; }

        // Gets the total photon count produced
        G4int GetPhotonCount() { return fPhotonCount_Scint + fPhotonCount_Ceren; }
        void IncPMTSAboveThreshold() { ++fPMTsAboveThreshold; }
        G4int GetPMTSAboveThreshold() { return fPMTsAboveThreshold; }


    private:
        G4CCMEventMessenger* fEventMessenger = nullptr;
        const G4CCMDetectorConstruction* fDetector = nullptr;

        G4int fScintCollID = -1;
        G4int fPMTCollID = -1;

        G4int fVerbose = 0;

        G4int fPMTThreshold = 1;

        G4bool fForcedrawphotons = false;
        G4bool fForcenophotons = false;

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
};
#endif

