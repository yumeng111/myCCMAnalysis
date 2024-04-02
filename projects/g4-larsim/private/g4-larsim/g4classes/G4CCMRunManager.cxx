// On Sun, to prevent conflict with ObjectSpace, G4Timer.hh has to be
// loaded *before* globals.hh...
#include "G4Timer.hh"

#include <g4-larsim/g4classes/G4CCMRunManager.h>
#include "g4-larsim/g4classes/G4CCMPMTSD.h"
#include <G4ParticleGun.hh>
#include <G4Run.hh>
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

G4CCMRunManager::G4CCMRunManager(): G4RunManager() {}

void G4CCMRunManager::InitializeRun()
{
    // Reset the event counter 
    numberOfEventToBeProcessed = 0;
    ConstructScoringWorlds();
    RunInitialization();

    if(verboseLevel>0) timer->Start();
}


void G4CCMRunManager::InjectParticle(G4ParticleGun* particleGun)
{
    if(!currentRun){
        G4String text = "Run needs to be initialized before injecting a particle.";
        G4Exception("G4CCMRunManager::InjectParticle()", "G4CCMRunManager002", FatalException, text);
    }
    assert(currentRun); // the G4Exception() above calls abort(). This assert() silences the clang static analyzer
    
    numberOfEventToBeProcessed++;
    currentRun->SetNumberOfEventToBeProcessed(numberOfEventToBeProcessed);

    currentEvent = GenerateEvent(numberOfEventToBeProcessed);
    particleGun->GeneratePrimaryVertex(currentEvent);

    G4EventManager* eventManager = G4EventManager::GetEventManager();
    eventManager->ProcessOneEvent(currentEvent);

    AnalyzeEvent(currentEvent);
    Update_Scoring();
    //StackPreviousEvent(currentEvent);
    //currentEvent = 0;


    // let's do one other thing
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String sdName = "/LAr/pmtSD";
    G4CCMPMTSD* pmtSD = (G4CCMPMTSD*) SDman->FindSensitiveDetector(sdName);
    CCMMCPEMap = pmtSD->GetCCMMCPEMap();

    if(runAborted) TerminateRun();
}


G4Event* G4CCMRunManager::GenerateEvent(G4int i_event)
{
    G4Event* anEvent = new G4Event(i_event);
    return anEvent;
}


void G4CCMRunManager::TerminateRun()
{
    if(verboseLevel>0){
        timer->Stop();
        G4cout << "Run terminated." << G4endl;
        G4cout << "Run Summary" << G4endl;
    if(runAborted){
        G4cout << "  Run Aborted after " << numberOfEventToBeProcessed << " events processed." << G4endl;
    }
    else{
        G4cout << "  Number of events processed : " << numberOfEventToBeProcessed << G4endl;
    }

    G4cout << "  "  << *timer << G4endl;
    }

    RunTermination();
}

//
// The following method is an exact copy of 
// UpdateScoring which is private in the G4RunManager
//

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
