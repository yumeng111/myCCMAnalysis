/*----------------------------------------------------------
 *
 *   CCMRootIO
 *
 *     Adapted from SBRootIO of SciBath
 *     Handles the input/output of ROOT objects
 *
 *     Adapter: R. T. Thornton (LANL)
 *     Date: May 7, 2020
 *
 *-----------------------------------------------------------*/
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <sstream>

#include "MsgLog.h"
#include "CCMRootIO.h"
#include "RawData.h"
#include "Pulses.h"
#include "AccumWaveform.h"
#include "Events.h"
#include "SimplifiedEvent.h"

#include "TFile.h"
#include "TObjString.h"

std::shared_ptr<CCMRootIO> CCMRootIO::fgInstance = nullptr;

//__________________________________________________
std::shared_ptr<CCMRootIO> CCMRootIO::GetInstance()
{
  if (fgInstance != nullptr) {
    return fgInstance;
  }

  fgInstance = std::shared_ptr<CCMRootIO>(new CCMRootIO());
  return fgInstance;
}

//__________________________________________________
CCMRootIO::CCMRootIO()
  : fEventHandle(std::make_unique<CCMEventTreeHandle>()),

  fFileIndex(-1),
  fInFile(0),
  fOwnHandle(true),
  fInFileList(0),
  fOutFileName(""),
  fOutFile(0),
  fFlushFreq(100),
  fOutSizeLimit(0)
{
    fTriggerNumber = 0;
    fReadOK = false;
    fNWrite = 0;
}

//__________________________________________________
CCMRootIO::~CCMRootIO()
{
  Close();

  if (fInFile) {
    fInFile = nullptr;
  }
  if(fOwnHandle && fEventHandle) {
    fEventHandle = nullptr;
  }

}

//__________________________________________________
void CCMRootIO::SetInFileName(const char* infile)
{
  fInFileList.clear();
  std::string ifName(infile);
  fInFileList.push_back(ifName);
  fFileIndex = 0;

}

//__________________________________________________
void CCMRootIO::SetInFileList(std::vector<std::string> infileList)
{
  fInFileList = infileList;
  if (fFileIndex<0 && !fInFileList.empty()) {
    fFileIndex = 0;
  }

}

//__________________________________________________
int CCMRootIO::AddFile(const char* file_regexp, bool hasWildcards)
{
  int nfiles = 0;
  if (hasWildcards) {
    fInFileList = Utility::GetListOfFiles(file_regexp);
  } else {
    fInFileList.push_back(std::string(file_regexp));
    ++nfiles;
  }

  // Setup the first file in the list
  if (fFileIndex<0 && fInFileList.size()>0) {
    fFileIndex = 0;
    fEventHandle->Clear();
    if (this->SetupInputFile()) {
      fReadOK = true;
    } else {
      fReadOK = false;
    }
  }
  return nfiles;

}

//__________________________________________________
int CCMRootIO::RemoveFile(const char* /*file_regexp*/)
{
  //not implemented yet
  return 0;
}

//__________________________________________________
int CCMRootIO::GoToFile(const char* file)
{
  unsigned int i;
  for (i=0; i< fInFileList.size(); ++i) {
    if (fInFileList[i] == file) {
      fFileIndex = i;
      this->SetupInputFile();
      return 1;
    }
  }
  // File not in list
  return 0;
}

//__________________________________________________
int CCMRootIO::AdvanceFile(int n)
{
  //======================================================================
  // Advance n positions in the file list
  //======================================================================
  if (n<=0) return 0;

  int indexMax  = fInFileList.size()-1;
  int indexSave = fFileIndex;

  fFileIndex += n;
  if (fFileIndex > indexMax) {
    fFileIndex = indexMax+1;
    fReadOK = false;
    return 0;
  }

  // Do something to open new file...
  this->SetupInputFile();

  return (fFileIndex-indexSave);
}

//__________________________________________________
int CCMRootIO::RewindFile(int n)
{
  //======================================================================
  // Rewind n positions in the file list
  //======================================================================
  if (n<=0) return 0;

  int indexMin  = 0;
  int indexSave = fFileIndex;

  fFileIndex -= n;
  if (fFileIndex < indexMin) {
    fFileIndex = 0;
    fReadOK = false;
    return 0;
  }

  // Do something to open new file...
  this->SetupInputFile();

  return (indexSave-fFileIndex);
}

//__________________________________________________
uint32_t CCMRootIO::GoTo(uint32_t event) 
{
  //Go to specified event

  //Get an initial guess at how many events to skip
  this->UpdateTriggerNumbers();

  //we're already there
  if (fTriggerNumber == event) return 1;

  //Advance/Rewind in the data stream
  while (event < fTriggerNumber) {
    uint32_t nrew = this->Rewind(1);
    this->UpdateTriggerNumbers();
    if (nrew == 0) break;
  }
  while (event > fTriggerNumber) {
    uint32_t nadv = this->Advance(1);
    this->UpdateTriggerNumbers();
    if (nadv == 0) break;
  }

  //check if we've found the right event number
  if (fTriggerNumber == event) 
    return 1;
  else
    return 0;

}

//__________________________________________________
uint32_t CCMRootIO::Advance(uint32_t n)
{
  //Advance n positions in the event stream, if possible
  if(!fEventHandle) {
    MsgError("No event tree handle!  Shouldn't happen.");
    exit(-1);
  }


  uint32_t ndone = 0;
  while (ndone < n) {
    //Try to advance far enough to match the request
    uint32_t ntry = n-ndone;
    uint32_t ndid = fEventHandle->Advance(ntry);
    ndone += ndid;
    //if not enough events in the file, return the last event in the file
    if (ndid < ntry) {
      fReadOK = false;
      return ndone;
    }
  }

  this->UpdateTriggerNumbers();
  return ndone;
}

//__________________________________________________
uint32_t CCMRootIO::Rewind(uint32_t n)
{
  //Advance n positions in the event stream, if possible
  if(!fEventHandle) {
    MsgError("No event tree handle!  Shouldn't happen.");
    exit(-1);
  }

  uint32_t ndone = 0;
  while (ndone < n) {
    //Try to rewind far enough to match the request
    uint32_t ntry = n-ndone;
    uint32_t ndid = fEventHandle->Rewind(ntry);
    ndone += ndid;
    //if not enough events in the file, return the last event in the file
    if (ndid < ntry) {
      fReadOK = false;
      return ndone;
    }
  }
  this->UpdateTriggerNumbers();
  return ndone;
}

//__________________________________________________
int CCMRootIO::Reload()
{
  //Mark the event handle as unfilled so that requests for data
  //members have to go back to the event file
  fEventHandle->ClearLoadFlags();
  return 1;
}


//__________________________________________________
int CCMRootIO::SetupInputFile()
{
  //set up the input file for reading
  if (fFileIndex<0 || fFileIndex>=(int)fInFileList.size()) {
    fReadOK = false;
    return 0;
  }

  if(fInFile) {
    fInFile->Flush();
    fInFile = nullptr;
    fReadOK = false;
  }

  fInFile = std::make_shared<TFile>(fInFileList[fFileIndex].c_str(),"READ");
  if(fInFile == nullptr) {
    MsgError(MsgLog::Form("Failed to open file: %s for read",fInFileList[fFileIndex].c_str()));
    fReadOK = false;
    return 0;
  }

  //Sync the event handle to the file
  if (fEventHandle->SetupInputFile(*fInFile)) {
    this->UpdateTriggerNumbers();
    fReadOK = true;
  } else {
    fReadOK = false;
  }

  return 1;

}

//__________________________________________________
void CCMRootIO::SetupOutputFile()
{
  //set up the output file for writing
  if(fEventHandle == nullptr) {
    fEventHandle = std::make_unique<CCMEventTreeHandle>();
  }

  if(fOutFile != nullptr) {
    fEventHandle->Close();
    fOutFile->Flush();
    fOutFile->Write();
    fOutFile = nullptr;
  }
  

  fOutFile = std::make_shared<TFile>(fOutFileName.c_str(),"RECREATE","CCM ROOT Event File",3);
  fEventHandle->SetupOutputFile(*fOutFile);
}


//__________________________________________________
const char* CCMRootIO::CurrentFileName() const
{
  if (fInFileList.size()>0) {
    if (fFileIndex>=0 && fFileIndex<(int)fInFileList.size()) {
      return fInFileList[fFileIndex].c_str();
    }
  }
  return "";
}

//__________________________________________________
const char* CCMRootIO::FileName(int i) const
{
  if (i>=0 && i<(int)fInFileList.size()) return fInFileList[i].c_str();
  return 0;
}

//__________________________________________________
CCMEventTreeHandle& CCMRootIO::GetEventTree()
{
  if (fEventHandle == nullptr) {
    fEventHandle = std::make_unique<CCMEventTreeHandle>();
  }
  return *fEventHandle;
}

//__________________________________________________
void CCMRootIO::UpdateTriggerNumbers()
{
  fTriggerNumber = this->GetEventTree().Index();
  fReadOK = true;

}

//__________________________________________________
int CCMRootIO::WriteTrigger()
{
  if (fOutFile == 0) {
    MsgDebug(2,"No output file set.");
    return 0;
  }

  int aok = 0;

  aok = fEventHandle->Write();

  if(aok) {
    ++fNWrite;
    if(fNWrite%fFlushFreq == 0) { 
      fOutFile->Flush(); 
    }
  }

  UpdateTriggerNumbers();

  return aok;

}

//__________________________________________________
void CCMRootIO::Close()
{

  if(fOutFile != nullptr) {
    MsgInfo(MsgLog::Form("Wrote %d event(s) to Event Tree.",fNWrite));
    MsgInfo("Closing File.");
    fEventHandle->Close();
    fOutFile->Flush();
    fOutFile->Write();
    fOutFile = nullptr;
  }

  if (fInFile != nullptr) {
    fInFile->Flush();
    fInFile = nullptr;
  }
}

//__________________________________________________
void CCMRootIO::Clear()
{
  if (fInFile != nullptr) {
    fInFile->Flush();
    fInFile = nullptr;
  }

  if (!fInFileList.empty()) {
    fInFileList.clear();
    fFileIndex = 0;
    fReadOK = false;
  }
}

//__________________________________________________
void CCMRootIO::Dump()
{
  //inspect the values of all variables
  printf("========================CCMRootIO::Dump()=========================\n");

  printf("CCMRootIO::fCurrentFile = %s\n",fInFileList[fFileIndex].c_str());

  printf("CCMRootIO::fTriggerNumber = %d\n",fTriggerNumber);

  printf("CCMRootIO::fOutFileName = %s\n",fOutFileName.c_str());

  printf("=====================================================================\n");

}

//__________________________________________________
AccumWaveform& CCMRootIO::GetAccumWaveform()
{
  return fEventHandle->CurrentAccumWaveform();
}

//__________________________________________________
Events& CCMRootIO::GetEvents()
{
  return fEventHandle->CurrentEvents();
}

//__________________________________________________
RawData& CCMRootIO::GetRawData()
{
  return fEventHandle->CurrentRawData();
}

//__________________________________________________
Pulses& CCMRootIO::GetPulses()
{
  return fEventHandle->CurrentPulses();
}

//__________________________________________________
void CCMRootIO::SetAccumWaveform(const AccumWaveform & accumWaveform)
{
  fEventHandle->SetCurrentAccumWaveform(accumWaveform);
}

//__________________________________________________
void CCMRootIO::SetEvents(const Events & events)
{
  fEventHandle->SetCurrentEvents(events);
}

//__________________________________________________
void CCMRootIO::SetRawData(const RawData & rawData)
{
  fEventHandle->SetCurrentRawData(rawData);
}

//__________________________________________________
void CCMRootIO::SetPulses(const Pulses & pulses)
{
  fEventHandle->SetCurrentPulses(pulses);
}

//__________________________________________________
void CCMRootIO::SetParameter(std::string /*name*/, const int /*value*/)
{
  // nothing yet
}

//--------------------------------------------------------------------
void CCMRootIO::SetParameter(std::string /*name*/, const double /*value*/)
{
  // nothing yet
}

//--------------------------------------------------------------------
void CCMRootIO::SetParameter(std::string name, std::string value)
{
  MsgInfo(MsgLog::Form("Name %s Value %s",name.c_str(),value.c_str()));
  if (name.find("ReadBranches") != std::string::npos) {
    std::stringstream ss(value);
    std::string branchName = "";
    std::vector<std::string> branches;
    while (ss >> branchName) {
      branches.emplace_back(branchName);
    }
    fEventHandle->SetReadBranches(branches);
  } else if (name.find("SaveBranches") != std::string::npos) {
    std::stringstream ss(value);
    std::string branchName = "";
    std::vector<std::string> branches;
    while (ss >> branchName) {
      branches.emplace_back(branchName);
    }
    fEventHandle->SetSaveBranches(branches);
  }
}

