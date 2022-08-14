/*----------------------------------------------------------
 *
 *   CCMRawIO
 *
 *     Adapted from SBRootIO of SciBath
 *     Handles the input/output of ROOT objects
 *
 *     Adapter: R. T. Thornton (LANL)
 *     Date: May 7, 2020
 *
 *-----------------------------------------------------------*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#include "CCMAnalysis/io/CCMRawIO.h"

#include "CCMAnalysis/ds/RawData.h"
#include "CCMAnalysis/utils/Utility.h"
#include "CCMAnalysis/utils/PMTInfoMap.h"
#include "CCMAnalysis/utils/MsgLog.h"

//__________________________________________________
CCMRawIO::CCMRawIO()
  : fFileIndex(-1),
  fInFile(0),
  fInFileList(0),
  fOutFileName(""),
  fOutFile(0),
  fFlushFreq(100),
  fOutSizeLimit(0),
  fRawData(std::make_unique<RawData>())
{
    fTriggerNumber = 0;
    fReadOK = false;
    fNWrite = 0;
}

//__________________________________________________
CCMRawIO::~CCMRawIO()
{
  Close();

  if (fInFile.is_open()) {
    fInFile.close();
  }

}

//__________________________________________________
void CCMRawIO::SetInFileName(const char* infile)
{
  fInFileList.clear();
  std::string ifName(infile);
  fInFileList.push_back(ifName);
  fFileIndex = 0;

}

//__________________________________________________
void CCMRawIO::SetInFileList(std::vector<std::string> infileList)
{
  fInFileList = infileList;
  if (fFileIndex<0 && fInFileList.size()>0) {
    fFileIndex = 0;
  }

}

//__________________________________________________
int CCMRawIO::AddFile(const char* file_regexp, bool hasWildcards)
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
    if (this->SetupInputFile()) {
      fReadOK = true;
    } else {
      fReadOK = false;
    }
  }
  return nfiles;
}

//__________________________________________________
int CCMRawIO::RemoveFile(const char* /*file_regexp*/)
{
  //not implemented yet
  return 0;
}

//__________________________________________________
int CCMRawIO::GoToFile(const char* file)
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
int CCMRawIO::AdvanceFile(int n)
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
int CCMRawIO::RewindFile(int n)
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
uint32_t CCMRawIO::GoTo(uint32_t event) 
{
  //Go to specified event

  //we're already there
  if (fTriggerNumber == event) return 1;

  //Advance/Rewind in the data stream
  while (event < fTriggerNumber) {
    uint32_t nrew = this->Rewind(1);
    if (nrew == 0) break;
  }
  while (event > fTriggerNumber) {
    uint32_t nadv = this->Advance(1);
    if (nadv == 0) break;
  }

  //check if we've found the right event number
  if (fTriggerNumber == event) 
    return 1;
  else
    return 0;

}

//__________________________________________________
uint32_t CCMRawIO::Advance(uint32_t n)
{
  //Advance n positions in the event stream, if possible
  if(!fInFile.is_open()) {
    MsgFatal("No file open! Need to run SetupInputFile first.");
  }

  long startPosition = fInFile.tellg();
  fInFile.seekg(0,fInFile.end);
  long totalLength = fInFile.tellg();

  long length = sizeof(event_t);
  long ndone = n;
  if(startPosition + length*ndone > totalLength) {
      ndone = (totalLength - startPosition) % length;
  }

  if (ndone <= 0) {
    fReadOK = false;
    return ndone;
  }

  // Reverted change that replaced "ndone" with "(ndone-1)". May cause issues
  fInFile.seekg(startPosition+length*(ndone-1));
  fInFile.read((char*)&fData, length);
  if (fInFile.eof()) {
    fReadOK = false;
    return ndone;
  }

  if (fInFile.fail()) {
    fReadOK = false;
    return ndone;
  }

  ConvertEventToRawData();

  fTriggerNumber += ndone;
  return ndone;
}

//__________________________________________________
uint32_t CCMRawIO::Rewind(uint32_t n)
{
  //Rewind n positions in the event stream, if possible
  if(!fInFile.is_open()) {
    MsgFatal("No file open! Need to run SetupInputFile first.");
  }

  int startPosition = fInFile.tellg();

  int length = sizeof(event_t);
  int ndone = n;
  while (startPosition - length*ndone < 0) {
    --ndone;
  }

  if (ndone == 0) {
    fReadOK = false;
    return ndone;
  }

  fInFile.seekg(startPosition-length*ndone);
  fInFile.read((char*)&fData, length);
  if (fInFile.eof()) {
    fReadOK = false;
    return ndone;
  }

  if (fInFile.fail()) {
    fReadOK = false;
    return ndone;
  }

  ConvertEventToRawData();

  fTriggerNumber -= ndone;
  return ndone;
}

//__________________________________________________
int CCMRawIO::Reload()
{
  //Mark the event handle as unfilled so that requests for data
  //members have to go back to the event file
  fInFile.seekg(0,fInFile.beg);
  fInFile.clear();
  return 1;
}


//__________________________________________________
int CCMRawIO::SetupInputFile()
{
  //set up the input file for reading
  if (fFileIndex<0 || fFileIndex>=(int)fInFileList.size()) {
    fReadOK = false;
    return 0;
  }

  if(fInFile) {
    fInFile.close();
    fReadOK = false;
  }

  fInFile.open(fInFileList[fFileIndex].c_str(), std::ios::binary);
  if(!fInFile.is_open()) {
    MsgFatal(MsgLog::Form("Failed to open file: %s for read",fInFileList[fFileIndex].c_str()));
  }
  if(!fInFile.good()) {
    MsgFatal(MsgLog::Form("file: %s is not good",fInFileList[fFileIndex].c_str()));
  }

  //Sync the event handle to the file
  fTriggerNumber = 0;
  fReadOK = true;

  return 1;

}

//__________________________________________________
void CCMRawIO::SetupOutputFile()
{
  //set up the output file for writing
  if(fOutFile.is_open()) {
    fOutFile.close();
  }
  
  fOutFile.open(fOutFileName.c_str(),std::ios::binary);

}


//__________________________________________________
const char* CCMRawIO::CurrentFileName() const
{
  if (fInFileList.size()>0) {
    if (fFileIndex>=0 && fFileIndex<(int)fInFileList.size()) {
      return fInFileList[fFileIndex].c_str();
    }
  }
  return "";
}

//__________________________________________________
const char* CCMRawIO::FileName(int i) const
{
  if (i>=0 && i<(int)fInFileList.size()) return fInFileList[i].c_str();
  return 0;
}

//__________________________________________________
int CCMRawIO::WriteTrigger()
{
  if (!fOutFile.is_open()) {
    MsgDebug(2,"No output file set.");
    return 0;
  }

  ConvertRawDataToEvent();

  fOutFile.write((char*)&fData,sizeof(event_t));

  ++fNWrite;
  if(fNWrite%fFlushFreq == 0) { 
    fOutFile << std::flush;
  }

  return 1;

}

//__________________________________________________
void CCMRawIO::Close()
{

  if (fOutFileName != "") {
    if(fOutFile.is_open()) {
      MsgInfo(MsgLog::Form("Wrote %d event(s) to Raw File.",fNWrite));
      MsgInfo("Closing OutFile.");
      fOutFile.close();
    }
  }

  if (!fInFileList.empty()) {
    if (fInFile.is_open()) {
      MsgInfo("Closing InFile.");
      fInFile.close();
    }
  }

}

//__________________________________________________
void CCMRawIO::Dump()
{
  //inspect the values of all variables
  printf("========================CCMRawIO::Dump()=========================\n");

  printf("CCMRawIO::fCurrentFile = %s\n",fInFileList[fFileIndex].c_str());

  printf("CCMRawIO::fTriggerNumber = %d\n",fTriggerNumber);

  printf("CCMRawIO::fOutFileName = %s\n",fOutFileName.c_str());

  printf("=====================================================================\n");

}

//__________________________________________________
void CCMRawIO::SetParameter(std::string /*name*/, const int /*value*/)
{
  // nothing yet
}

//--------------------------------------------------------------------
void CCMRawIO::SetParameter(std::string /*name*/, const double /*value*/)
{
  // nothing yet
}

//--------------------------------------------------------------------
void CCMRawIO::SetParameter(std::string /*name*/, std::string /*value*/)
{
  // nothing yet
}

//--------------------------------------------------------------------
void CCMRawIO::ConvertEventToRawData()
{
  fRawData->Reset(NDIGITIZERS,
      NCHANNELS,
      NSAMPLES,
      fData.evNum,
      fData.computerTime.tv_sec,
      fData.computerTime.tv_nsec);

  fRawData->SetBoardEventNum(fData.digitizers.evNum);
  fRawData->SetClockTime(fData.digitizers.time);

  for (int board = 0; board < NDIGITIZERS; ++board) {
    fRawData->SetChannelMask(board,fData.digitizers.chMask[board]);
    fRawData->SetSize(board,fData.digitizers.size[board]);
    fRawData->SetTemp(board,fData.digitizers.temperatures[board]);

    // loop through the channels
    for (int channel = 0; channel < NCHANNELS; ++channel) {
      auto key = PMTInfoMap::CreateKey(board, channel);
      //MsgInfo(MsgLog::Form("Board  %d Event Number %d",board,fData.digitizers.samples[board][channel]));
      fRawData->SetWaveform(key,fData.digitizers.samples[board][channel]);
    } // end for channel
  } // end for board

}

//--------------------------------------------------------------------
void CCMRawIO::ConvertRawDataToEvent()
{
  fData.evNum = fRawData->GetEventNumber();
  fData.computerTime.tv_nsec = fRawData->GetGPSNSIntoSec();
  fData.computerTime.tv_sec = fRawData->GetGPSSecIntoDay();

  fRawData->SetBoardEventNum(fData.digitizers.evNum);
  fRawData->SetClockTime(fData.digitizers.time);

  for (int board = 0; board < NDIGITIZERS; ++board) {
    fData.digitizers.time[board] = fRawData->GetClockTime(board);
    fData.digitizers.evNum[board] = fRawData->GetBoardEventNum(board);

    // loop through the channels
    for (int channel = 0; channel < NCHANNELS; ++channel) {
      fData.digitizers.chMask[board][channel] = fRawData->GetChannelMask(board,channel);
      fData.digitizers.size[board][channel] = fRawData->GetSize(board,channel);
      fData.digitizers.temperatures[board][channel] = fRawData->GetTemp(board,channel);

      // loop through the samples
      int key = board*NCHANNELS+channel;
      for (int sample = 0; sample < NSAMPLES; ++sample) {
        fData.digitizers.samples[board][channel][sample] = fRawData->GetSample(key,sample);
      }
    } // end for channel
  } // end for board
}

//__________________________________________________
RawData& CCMRawIO::GetRawData()
{
  return *fRawData;
}

//__________________________________________________
void CCMRawIO::SetRawData(const RawData & rawData)
{
  fRawData->operator=(rawData);
}
