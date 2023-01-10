/*----------------------------------------------------------
 *
 *   CCMEventTreeHandle
 *
 *     Based on EDMEventHandle and IoEventHandle of E907/MIPP
 *
 *       Provide an interface between the CCMRootIO object and
 *         the RawData (ROOT streamable) object.
 *           Allows for partial IO of various major compoents of the event
 *
 *             Adapter: J.L. Klay (CalPoly)
 *               Date: 11-Sep-2008
 *
 *-----------------------------------------------------------*/

#include <limits>
#include <TProcessID.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <algorithm>

#include "CCMAnalysis/CCMIO/CCMEventTreeHandle.h"

#include "CCMAnalysis/CCMDataStructures/Events.h"
#include "CCMAnalysis/CCMDataStructures/Pulses.h"
#include "CCMAnalysis/CCMDataStructures/MCTruth.h"
#include "CCMAnalysis/CCMDataStructures/RawData.h"
#include "CCMAnalysis/CCMDataStructures/AccumWaveform.h"
#include "CCMAnalysis/CCMDataStructures/SimplifiedEvent.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

static const char* gsRawData = "rawData";
static const char* gsPulses = "pulses";
static const char* gsEvents= "EventTree";

//_____________________________________________________________
CCMEventTreeHandle::CCMEventTreeHandle() 
: fIndex(0), fNumOfEntries(0), fInputFile(), fOutputFile(), 
  fEventsTree(nullptr), fPulsesTree(nullptr), fOutEventTree(nullptr),
  fObjectCnt(TProcessID::GetObjectCount()),
  fRawData(new RawData()),
  fPulses(new Pulses()),
  fEvents(new Events()),
  fAccumWaveform(new AccumWaveform()),
  fMCTruth(new MCTruth())
{
  //constructor
  for (int i = 0; i < kNEventBranch; i++) {
    fIsLoaded[i] = false;
    fBranch[i] = nullptr;
  }

  Long64_t tsize = 150; //set the limit of the event tree to be large (GB)
  for (int i = 0; i < 3; i++) tsize *= 1000; //Convert to bytes from GB
  TTree::SetMaxTreeSize(tsize);

  this->ClearLoadFlags();
}

//_____________________________________________________________
CCMEventTreeHandle::~CCMEventTreeHandle()
{
  //destructor
}

//_____________________________________________________________
void CCMEventTreeHandle::Clear()
{
  this->ClearLoadFlags();
}

//_____________________________________________________________
int CCMEventTreeHandle::ClearLoadFlags()
{
  //Mark all branches as unloaded.  Next read attempt will have to go
  //back to disk.  Return the set of load flags prior to the clear.

  int isave = 0;
  TProcessID::SetObjectCount(fObjectCnt);

  for (int i = 0; i < kNEventBranch; i++) {
    //mark all branches as unloaded
    if (fIsLoaded[i] == true) isave |= (1 << i);
    fIsLoaded[i] = false;
  }
  return isave;

}

//_____________________________________________________________
void CCMEventTreeHandle::SetLoaded(int branchId) const
{
  if (branchId < 0 || branchId > kNEventBranch) 
    //MsgFatal("Requested branch not defined.  Should not happen.");
    ::abort();

  //Use kNEventBranch to mark that all branches are loaded
  if(branchId == kNEventBranch) {
    for (int i = 0; i < kNEventBranch; i++) fIsLoaded[i] = true;
    return;
  }

  //Use others to mark branches one-by-one
  fIsLoaded[branchId] = true;

}

//_____________________________________________________________
bool CCMEventTreeHandle::IsLoaded(CCMEventBranchID_t id) const
{
  if (id < kNEventBranch && id >= 0) return fIsLoaded[id];
  //MsgFatal("Requested branchId does not exist.  Should not happen.");
  return false;
}

//_____________________________________________________________
bool CCMEventTreeHandle::ReadTree() // Is none not in the config
{
  if (std::find(fReadBranches.begin(),fReadBranches.end(),"none") != fReadBranches.end()) {
    return false;
  }

  return true;
}

//_____________________________________________________________
int CCMEventTreeHandle::SetupInputFile(TFile & f)
{
  // set up an event file to read events
  fInputFile = &f;

  if (MsgLog::GetGlobalDebugLevel() >= 2) {
    MsgDebug(2,"Printing Input File");
    fInputFile->Print();
  }

  /*
  if (!ReadTree()) { // If none is in the config
    for ( auto & branch : fBranch) {
      branch = nullptr;
    }

    fNumOfEntries = 0;

    //Reset to the start of the file
    fIndex = 0;

    //Flag all branches as unloaded
    this->ClearLoadFlags();

    return 1;
  }
  */

  //In with the new...
  fInputFile->GetObject(gsRawData,fEventsTree);

  if (fEventsTree == nullptr) {
    fInputFile->GetObject(gsEvents,fEventsTree);
  } else {
    fInputFile->GetObject(gsPulses,fPulsesTree);
  }
  if(fEventsTree == nullptr) {
    fInputFile->ls();
    MsgError(MsgLog::Form("Failed to find %s or %s in file %s",gsRawData,gsEvents,fInputFile->GetName()));
    return 0;
  }

  if (MsgLog::GetGlobalDebugLevel() >= 2) {
    MsgDebug(2,"Printing Events Tree");
    fEventsTree->Print();
    if (fPulsesTree != nullptr) {
      MsgDebug(2,"Printing Pulses Tree");
      fPulsesTree->Print();
    }
  }

  TBranch * branch = 0;

  if (std::find(fReadBranches.begin(),fReadBranches.end(),"rawData") != fReadBranches.end() || 
      std::find(fReadBranches.begin(),fReadBranches.end(),"all") != fReadBranches.end()) {
    branch = fBranch[kRawDataID] = fEventsTree->GetBranch("rawData");
    if(branch != 0) {
      fBranch[kRawDataID]->SetAddress(&fRawData);
    }
  } else {
    fBranch[kRawDataID] = nullptr;
  }

  if (std::find(fReadBranches.begin(),fReadBranches.end(),"pulses") != fReadBranches.end() || 
      std::find(fReadBranches.begin(),fReadBranches.end(),"all") != fReadBranches.end()) {
    if (fPulsesTree == nullptr) {
      branch = fBranch[kPulsesID] = fEventsTree->GetBranch("pulses");
    } else {
      branch = fBranch[kPulsesID] = fPulsesTree->GetBranch("pulses");
    }
    if (branch != 0) {
      fBranch[kPulsesID]->SetAddress(&fPulses);
    }
  } else {
    fBranch[kPulsesID] = nullptr;
  }

  if (std::find(fReadBranches.begin(),fReadBranches.end(),"events") != fReadBranches.end() || 
      std::find(fReadBranches.begin(),fReadBranches.end(),"all") != fReadBranches.end()) {
    branch = fBranch[kEventsID] = fEventsTree->GetBranch("events");
    if(branch != 0) {
      fBranch[kEventsID]->SetAddress(&fEvents);
    }
  } else {
    fBranch[kEventsID] = nullptr;
  }

  if (std::find(fReadBranches.begin(),fReadBranches.end(),"accumWaveform") != fReadBranches.end() || 
      std::find(fReadBranches.begin(),fReadBranches.end(),"all") != fReadBranches.end()) {
    branch = fBranch[kAccumWaveformID] = fEventsTree->GetBranch("accumWaveform");
    if(branch != 0) {
      fBranch[kAccumWaveformID]->SetAddress(&fAccumWaveform);
    }
  } else {
    fBranch[kAccumWaveformID] = nullptr;
  }

  if (std::find(fReadBranches.begin(),fReadBranches.end(),"mcTruth") != fReadBranches.end() || 
      std::find(fReadBranches.begin(),fReadBranches.end(),"all") != fReadBranches.end()) {
    branch = fBranch[kMCTruthID] = fEventsTree->GetBranch("mcTruth");
    if(branch != 0) {
      fBranch[kMCTruthID]->SetAddress(&fMCTruth);
    }
  } else {
    fBranch[kMCTruthID] = nullptr;
  }

  if (MsgLog::GetGlobalDebugLevel() >= 2) {
    MsgDebug(2,"Printing Branches");
    for (const auto & branch : fBranch) {
      if (branch != nullptr) {
        branch->Print();
      }
    }
  }

  fNumOfEntries = (uint32_t)fEventsTree->GetEntries();

  //Reset to the start of the file
  fIndex = 0;

  //Flag all branches as unloaded
  this->ClearLoadFlags();

  return 1;

}

//_____________________________________________________________
int CCMEventTreeHandle::SetupOutputFile(TFile & f)
{
  //Out with the old
  if (fOutEventTree != nullptr) {
    fOutEventTree = nullptr;
  }

  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"none") != fSaveBranches.end()) {
    fOutputFile = nullptr;
    return 1;
  }

  //Set up an output file for writing events
  fOutputFile = &f;
  
  fOutputFile->cd();
  fOutEventTree = new TTree(gsEvents, gsEvents);
  if(fOutEventTree == nullptr) {
    MsgError("Problem setting up output tree!");
    return 0;
  }

  int sLvl = fgkSplitLevel;

  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"events") != fSaveBranches.end() ||
      std::find(fSaveBranches.begin(),fSaveBranches.end(),"all") != fSaveBranches.end()) {
    if (fEvents == nullptr) {
      fEvents = new Events();
    }
    fOutEventTree->Branch("events", &fEvents, 320000, sLvl);
  }
  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"accumWaveform") != fSaveBranches.end() ||
      std::find(fSaveBranches.begin(),fSaveBranches.end(),"all") != fSaveBranches.end()) {
    if (fAccumWaveform == nullptr) {
      fAccumWaveform = new AccumWaveform();
    }
    fOutEventTree->Branch("accumWaveform", &fAccumWaveform, 320000, sLvl);
  }
  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"mcTruth") != fSaveBranches.end() ||
      std::find(fSaveBranches.begin(),fSaveBranches.end(),"all") != fSaveBranches.end()) {
    if (fMCTruth == nullptr) {
      fMCTruth = new MCTruth();
    }
    fOutEventTree->Branch("mcTruth", &fMCTruth, 320000, sLvl);
  }
  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"rawData") != fSaveBranches.end() ||
      std::find(fSaveBranches.begin(),fSaveBranches.end(),"all") != fSaveBranches.end()) {
    if (fRawData == nullptr) {
      fRawData = new RawData();
    }
    fOutEventTree->Branch("rawData", &fRawData, 320000, sLvl);
  }
  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"pulses") != fSaveBranches.end() ||
      std::find(fSaveBranches.begin(),fSaveBranches.end(),"all") != fSaveBranches.end()) {
    if (fPulses == nullptr) {
      fPulses = new Pulses();
    }
    fOutEventTree->Branch("pulses", &fPulses, 320000, sLvl);
  }

  return 1;
}

//_____________________________________________________________
uint32_t CCMEventTreeHandle::Advance(uint32_t n)
{
  //advance n places in the input stream.  Return the difference
  //between the new position and the old.

  if(fEventsTree == nullptr || n < 1) {
    return 0;
  }

  uint32_t indexSave = fIndex;

  fIndex += n;
  if(fIndex >= fNumOfEntries) {
    fIndex = fNumOfEntries - 1;
  }

  //Since we've moved on, mark all branches as unloaded
  this->ClearLoadFlags();

  return fIndex-indexSave;
}

//_____________________________________________________________
uint32_t CCMEventTreeHandle::Rewind(uint32_t n)
{
  //rewind n places in the input stream.  Return the difference
  //between the new position and the old

  if(fRawData == nullptr || n < 1) {
    return 0;
  }

  uint32_t indexSave = fIndex;

  fIndex -= n;

  //Since we've moved on, mark all branches as unloaded
  this->ClearLoadFlags();

  return indexSave-fIndex;
}


//_____________________________________________________________
int CCMEventTreeHandle::Load(int branchID) const
{
  //Loads one branch if branchID < kNEventBranch or loads entire event if
  //branchId == kNEventBranch

  //Do not load anything if ROOT input file not specified
  if (fInputFile == nullptr) {
    return 0;
  }

  if(branchID < 0 || branchID > kNEventBranch) {
    abort();
  }

  //if branchID is maxed out, load all the branches. Load one-by-one
  //to avoid overwirting branches already loaded

  if(branchID == kNEventBranch){
    for (int i=0; i<kNEventBranch; i++) {
      this->Load(i);
    }
    return 1;
  }

  //check if branch is already loaded
  if (this->IsLoaded((CCMEventBranchID_t) branchID)) {
    return 0;
  }

  fInputFile->cd();
  //Make sure that this branch exists before actually loading it
  if(fBranch[branchID] != nullptr) {
    fBranch[branchID]->GetEntry(fIndex);
    this->SetLoaded(branchID);
    return 1;
  }

  return 0;
}

//_____________________________________________________________
int CCMEventTreeHandle::Write()
{
  if (std::find(fSaveBranches.begin(),fSaveBranches.end(),"none") != fSaveBranches.end()) {
    return 1;
  }

  if (fOutputFile != nullptr && fOutEventTree != nullptr) {
    //Make sure we load any branches which haven't been loaded yet so
    //they get written

    this->Load(kNEventBranch);

    fOutEventTree->Fill();
  }
  return 1;
}

//_____________________________________________________________
void CCMEventTreeHandle::Close()
{
  //MsgDebug(3,"Calling Close method");

  if (fOutputFile != nullptr && fOutEventTree != nullptr) {
    fOutputFile->cd();
    fOutEventTree->Write();
    fOutEventTree = nullptr;
  }

}

//_____________________________________________________________
AccumWaveform& CCMEventTreeHandle::CurrentAccumWaveform()
{
  this->Load(kAccumWaveformID);
  return *fAccumWaveform;
}

//_____________________________________________________________
MCTruth& CCMEventTreeHandle::CurrentMCTruth()
{
  this->Load(kMCTruthID);
  return *fMCTruth;
}

//_____________________________________________________________
Events& CCMEventTreeHandle::CurrentEvents()
{
  this->Load(kEventsID);
  return *fEvents;
}

//_____________________________________________________________
RawData& CCMEventTreeHandle::CurrentRawData()
{
  this->Load(kRawDataID);
  return *fRawData;
}

//_____________________________________________________________
Pulses& CCMEventTreeHandle::CurrentPulses()
{
  this->Load(kPulsesID);
  return *fPulses;
}

//_____________________________________________________________
void CCMEventTreeHandle::SetCurrentAccumWaveform(const AccumWaveform& accumWaveform)
{
  if (fAccumWaveform == nullptr) {
    fAccumWaveform = new AccumWaveform(accumWaveform);
    return;
  }
  fAccumWaveform->operator=(accumWaveform);
}

//_____________________________________________________________
void CCMEventTreeHandle::SetCurrentMCTruth(const MCTruth& mcTruth)
{
  if (fMCTruth == nullptr) {
    fMCTruth = new MCTruth(mcTruth);
    return;
  }
  fMCTruth->operator=(mcTruth);
}

//_____________________________________________________________
void CCMEventTreeHandle::SetCurrentEvents(const Events& events)
{
  if (fEvents == nullptr) {
    fEvents = new Events(events);
    return;
  }
  fEvents->operator=(events);
}

//_____________________________________________________________
void CCMEventTreeHandle::SetCurrentRawData(const RawData& rawData)
{
  if (fRawData == nullptr) {
    fRawData = new RawData(rawData);
    return;
  }
  fRawData->operator=(rawData);
}

//_____________________________________________________________
void CCMEventTreeHandle::SetCurrentPulses(const Pulses& pulses)
{
  if (fPulses == nullptr) {
    fPulses = new Pulses(pulses);
  }
  fPulses->operator=(pulses);
}

//_____________________________________________________________
void CCMEventTreeHandle::Report()
{
  if (fEventsTree != nullptr) {
    MsgInfo("fEventsTree: ");
    fRawData->Print();
  }
  if(fOutEventTree != nullptr) {
    MsgInfo("OutEventTree: ");
    fOutEventTree->Print();
  }
}

