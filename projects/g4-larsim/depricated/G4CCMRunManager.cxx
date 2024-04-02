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

const std::unordered_map<G4double, CCMMCPE::PhotonSource> G4CCMRunManager::G4doubletoPhotonSource = {{0., CCMMCPE::PhotonSource::Unknown},
                                                                                                     {1., CCMMCPE::PhotonSource::Scintillation},
                                                                                                     {2., CCMMCPE::PhotonSource::Cerenkov}};

G4CCMRunManager::G4CCMRunManager(): G4RunManager() {}


void G4CCMRunManager::InitializeRun()
{
    // Reset the event counter 
    numberOfEventToBeProcessed = 0;
    ConstructScoringWorlds();
    RunInitialization();

    if(verboseLevel>0) timer->Start();
    // some stuff from EventAction
    fHitCount                = 0;
    fPhotonCount_Scint       = 0;
    fPhotonCount_Ceren       = 0;
    fAbsorptionCount         = 0;
    fBoundaryAbsorptionCount = 0;
    fTotE                    = 0.0;

    fConvPosSet = false;
    fEdepMax    = 0.0;

    fPMTsAboveThreshold = 0;

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //if(fScintCollID < 0)
    //  fScintCollID = SDman->GetCollectionID("scintCollection");
    if(fPMTCollID < 0)
        fPMTCollID = SDman->GetCollectionID("pmtHitCollection");
    // end of stuff from EventAction
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
    GetFinalScores(currentEvent);
    //StackPreviousEvent(currentEvent);
    //currentEvent = 0;


    // let's do one other thing
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String sdName = "/LAr/pmtSD";
    G4CCMPMTSD* pmtSD = (G4CCMPMTSD*) SDman->FindSensitiveDetector(sdName);
    std::vector<CCMPMTKey> pmt_keys = pmtSD->GetPMTKeys();
    std::cout << "pmt keys = " << pmt_keys << std::endl;



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

void G4CCMRunManager::GetFinalScores(const G4Event* anEvent) {
    //G4CCMScintHitsCollection* scintHC = nullptr;
    G4CCMPMTHitsCollection* pmtHC     = nullptr;
    G4HCofThisEvent* hitsCE         = anEvent->GetHCofThisEvent();

    // Get the hit collections
    if(hitsCE) {
        //if(fScintCollID >= 0) {
        //    scintHC = (G4CCMScintHitsCollection*) (hitsCE->GetHC(fScintCollID));
        //}
        if(fPMTCollID >= 0) {
            pmtHC = (G4CCMPMTHitsCollection*) (hitsCE->GetHC(fPMTCollID));
        }
    }

    // Hits in scintillator
    //if(scintHC) {
    //    size_t n_hit = scintHC->entries();
    //    G4ThreeVector eWeightPos(0.);
    //    G4double edep;
    //    G4double edepMax = 0;

    //    for(size_t i = 0; i < n_hit; ++i) {  // gather info on hits in scintillator
    //        edep = (*scintHC)[i]->GetEdep();
    //        fTotE += edep;
    //        eWeightPos += (*scintHC)[i]->GetPos() * edep;  // calculate energy weighted pos
    //        
    //        if(edep > edepMax) {
    //            edepMax = edep;  // store max energy deposit
    //            G4ThreeVector posMax = (*scintHC)[i]->GetPos();
    //            fPosMax              = posMax;
    //            fEdepMax             = edep;
    //        }
    //    }

    //    if(fTotE == 0.) {
    //        if(fVerbose > 0)
    //            G4cout << "No hits in the scintillator this event." << G4endl;
    //    } else {
    //        // Finish calculation of energy weighted position
    //        eWeightPos /= fTotE;
    //        fEWeightPos = eWeightPos;
    //        if(fVerbose > 0) {
    //            G4cout << "\tEnergy weighted position of hits in G4CCM : " << eWeightPos / mm << G4endl;
    //        }
    //    }
    //    if(fVerbose > 0) {
    //        G4cout << "\tTotal energy deposition in scintillator : " << fTotE / keV << " (keV)" << G4endl;
    //    }
    //}

    // hits in pmt
    if(pmtHC) {
        G4ThreeVector reconPos(0., 0., 0.);
        size_t pmts = pmtHC->entries();
        // Gather info from all PMTs
        for(size_t i = 0; i < pmts; ++i) {
            // let's grab information necessary to make a CCMPMTKey
            int pmt_row = static_cast<int>((*pmtHC)[i]->GetCCMPMTKeyRow());
            int pmt_number = static_cast<int>((*pmtHC)[i]->GetCCMPMTKeyNumber());
            fKey = CCMPMTKey(pmt_row, pmt_number);
            // now let's grab information necessary to make a CCMMCPE
            std::vector<CCMMCPE> all_CCMMCPE;

            for (size_t i = 0; i < (*pmtHC)[i]->GetPhotonTime().size(); i++){
                double time = static_cast<double>((*pmtHC)[i]->GetPhotonTime()[i] * CLHEP::ns);          
                double energy = static_cast<double>((*pmtHC)[i]->GetPhotonEnergy()[i] * CLHEP::eV);
                double position_x = static_cast<double>((*pmtHC)[i]->GetPhotonPositionX()[i] * CLHEP::cm);          
                double position_y = static_cast<double>((*pmtHC)[i]->GetPhotonPositionY()[i] * CLHEP::cm);          
                double position_z = static_cast<double>((*pmtHC)[i]->GetPhotonPositionZ()[i] * CLHEP::cm);          
                double direction_x = static_cast<double>((*pmtHC)[i]->GetPhotonDirectionX()[i] * CLHEP::cm);          
                double direction_y = static_cast<double>((*pmtHC)[i]->GetPhotonDirectionY()[i] * CLHEP::cm);          
                double direction_z = static_cast<double>((*pmtHC)[i]->GetPhotonDirectionZ()[i] * CLHEP::cm);          
                CCMMCPE::PhotonSource source = G4doubletoPhotonSource.at((*pmtHC)[i]->GetPhotonCreationProcess()[i]);

                // now make a CCMMCPE!
                I3Position position(position_x, position_y, position_z);
                I3Direction direction(direction_x, direction_y, direction_z);
                CCMMCPE this_mc_pe = CCMMCPE(time, energy, position, direction, source);
                all_CCMMCPE.push_back(this_mc_pe);
            }

            // now save to I3Map
            CCMMCPEMap->insert(std::make_pair(fKey, all_CCMMCPE));

            // now visualization things
            fHitCount += (*pmtHC)[i]->GetPhotonCount();
            reconPos += (*pmtHC)[i]->GetPMTPos() * (*pmtHC)[i]->GetPhotonCount();
            if((*pmtHC)[i]->GetPhotonCount() >= fPMTThreshold) {
                ++fPMTsAboveThreshold;
            }
            else {  // wasn't above the threshold, turn it back off
                (*pmtHC)[i]->SetDrawit(false);
            }
        }


        if(fHitCount > 0) {  // don't bother unless there were hits
            reconPos /= fHitCount;
            if(fVerbose > 0) {
                G4cout << "\tReconstructed position of hits in CCM : " << reconPos / mm << G4endl;
            }
        fReconPos = reconPos;
        }
        pmtHC->DrawAllHits();
    }


    if(fVerbose > 0) {
        // End of event output. later to be controlled by a verbose level
        G4cout << "\tNumber of photons that hit PMTs in this event : " << fHitCount << G4endl;
        G4cout << "\tNumber of PMTs above threshold(" << fPMTThreshold << ") : " << fPMTsAboveThreshold << G4endl;
        G4cout << "\tNumber of photons produced by scintillation in this event : " << fPhotonCount_Scint << G4endl;
        G4cout << "\tNumber of photons produced by cerenkov in this event : " << fPhotonCount_Ceren << G4endl;
        G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : " << fAbsorptionCount << G4endl;
        G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in " << "this event : " << fBoundaryAbsorptionCount << G4endl;
        G4cout << "Unaccounted for photons in this event : " << (fPhotonCount_Ceren - fAbsorptionCount - fHitCount - fBoundaryAbsorptionCount) << G4endl;
    }

    //// update the run statistics
    //auto run = static_cast<G4CCMRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

    //run->IncHitCount(fHitCount);
    //run->IncPhotonCount_Scint(fPhotonCount_Scint);
    //run->IncPhotonCount_Ceren(fPhotonCount_Ceren);
    //run->IncEDep(fTotE);
    //run->IncAbsorption(fAbsorptionCount);
    //run->IncBoundaryAbsorption(fBoundaryAbsorptionCount);
    //run->IncHitsAboveThreshold(fPMTsAboveThreshold);

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
