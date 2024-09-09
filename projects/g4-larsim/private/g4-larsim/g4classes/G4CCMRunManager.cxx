// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...

#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include "g4-larsim/g4classes/G4CCMScintSD.h"
#include "g4-larsim/g4classes/G4CCMRunManager.h"

#include <G4Run.hh>
#include <G4Timer.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4ParticleGun.hh>
#include <G4EventManager.hh>

G4CCMRunManager::G4CCMRunManager(): G4RunManager() {}

void G4CCMRunManager::SimulateEvent(const I3Particle& primary, I3MCTreePtr tree, CCMMCPESeriesMapPtr mcpeseries) {
    // Reset the event counter
    fakeRun = false;
    numberOfEventToBeProcessed = 1;
    numberOfEventProcessed = 0;
    ConstructScoringWorlds();
    RunInitialization();

    if(!currentRun) {
        G4String text = "Run needs to be initialized before simulating an event.";
        G4Exception("G4CCMRunManager::SimulateEvent()", "G4CCMRunManager001", FatalException, text);
    }
    assert(currentRun); // the G4Exception() above calls abort(). This assert() silences the clang static analyzer

    DoEventLoop(numberOfEventToBeProcessed, nullptr, -1);
    RunTermination();
}

void G4CCMRunManager::SimulateEvents(std::vector<I3Particle> const & primaries, std::vector<I3MCTreePtr> trees, std::vector<CCMMCPESeriesMapPtr> mcpeseries) {
    // Reset the event counter
    fakeRun = false;
    numberOfEventToBeProcessed = primaries.size();
    numberOfEventProcessed = 0;
    ConstructScoringWorlds();
    RunInitialization();

    if(!currentRun) {
        G4String text = "Run needs to be initialized before simulating an event.";
        G4Exception("G4CCMRunManager::SimulateEvent()", "G4CCMRunManager001", FatalException, text);
    }
    assert(currentRun); // the G4Exception() above calls abort(). This assert() silences the clang static analyzer

    DoEventLoop(numberOfEventToBeProcessed, nullptr, -1);
    RunTermination();
}

#include <G4ScoringManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4VHitsCollection.hh>

void G4CCMRunManager::Update_Scoring() {
    G4ScoringManager* ScM = G4ScoringManager::GetScoringManagerIfExist();
    if(!ScM) return;
    G4int nPar = ScM->GetNumberOfMesh();
    if(nPar<1) return;

    G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
    if(!HCE) return;
    G4int nColl = HCE->GetCapacity();
    for(G4int i=0;i<nColl;i++)
    {
    G4VHitsCollection* HC = HCE->GetHC(i);
    if(HC) ScM->Accumulate(HC);
    }
}
