/*!**********************************************
 * \file CCMFindPulses
 * \author R. T. Thornton (LANL)
 * \date February 24, 2020
 * \brief Read in a binary, find pulses, save as a ROOT file
 *
 * Read in a binary file and find all the pulses for each trigger
 * save the output file as a ROOT file with two trees. The first tree
 * is the raw waveform (currently not saved) and the pulse information
 * is in a second tree. The output file is compressed using a
 * ROOT compression scheme.
 *
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMFindPulses.h"
#include "CCMModuleTable.h"

#include "MsgLog.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "Pulses.h"
#include "RawData.h"

#include "data_structures.hh"

#include "TF1.h"
#include "TFile.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TSQLServer.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

//See CCMModuleTable for info
MODULE_DECL(CCMFindPulses);

//_______________________________________________________________________________________
CCMFindPulses::CCMFindPulses(const char* version) 
  : CCMModule("CCMFindPulses"),
    fTriggerType("ALL"),
    fTruncateWaveform(false),
    fFromRootFile(false),
    fResetPulses(true),
    fWriteDBEntry(false),
    fDBHost(""),
    fDBUser(""),
    fDBPwd("")
{
  //Default constructor
  this->SetCfgVersion(version);

  if (!fTriggerTime.empty()) {
    fTriggerTime.clear();
  }
}

//_______________________________________________________________________________________
CCMFindPulses::CCMFindPulses(const CCMFindPulses& clufdr) 
: CCMModule(clufdr),
  fTriggerType(clufdr.fTriggerType),
  fTruncateWaveform(clufdr.fTruncateWaveform),
  fFromRootFile(clufdr.fFromRootFile),
  fResetPulses(clufdr.fResetPulses),
  fWriteDBEntry(clufdr.fWriteDBEntry),
  fDBHost(clufdr.fDBHost),
  fDBUser(clufdr.fDBUser),
  fDBPwd(clufdr.fDBPwd),
  fTriggerTime(clufdr.fTriggerTime),
  fReadData(clufdr.fReadData),
  fPulses(clufdr.fPulses),
  fRawData(clufdr.fRawData)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMFindPulses::~CCMFindPulses()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMFindPulses::ProcessTrigger()
{

  RawData localCopy(*fReadData);
  if (fFromRootFile) {
    localCopy = *fRawData;
  }

  if (!localCopy.IsTriggerPresent(fTriggerType)) {
    return kCCMDoNotWrite;
  }

  // fTriggerTime stores the time each board saw the trigger that is saved in channel 15
  ++fNEventsTotal;

  size_t board = 0;
  size_t channel = 0;
  size_t sample = 0;

  const size_t kNDigitizers = localCopy.GetNumBoards();
  const size_t kNChannels = localCopy.GetNumChannels();
  const size_t kNSamples = localCopy.GetNumSamples();

  const size_t kNPMTs = localCopy.GetNumWaveforms();
  const size_t kNEffDigitizers = (kNPMTs == kNChannels) ? 1 : kNDigitizers;

  fTriggerTime.assign(kNDigitizers,9e9);


  // check to make sure all the board event numbers are the same
  int channelEvtNum0 = localCopy.GetBoardEventNum(0);
  bool evNumCheck = true;
  for (board = 0; board < kNDigitizers && !fFromRootFile; ++board) {
    int channelEvtNum = localCopy.GetBoardEventNum(board);
    if (channelEvtNum != channelEvtNum0) {
      evNumCheck = false;
    }
  }

  // count number of times we have a mismatched trigger
  if (!evNumCheck) {
    MsgError(MsgLog::Form("Mismatched Trigger %ld",fNEventsTotal));
    ++fNEventsSkipped;
    //continue; // do not skip because of information for the database
  }


  /////////////////////////////////////////////////////////////////
  // Calculate the trigger time of each board in the DAQ window
  // Use the earliest trigger time as the time of the trigger
  // and shift all the other boards to line up with the earliest board
  /////////////////////////////////////////////////////////////////

  double minTriggerTime = 9e9;
  double maxTriggerTime = 9e9;
  // used if not saving trigger to channel 15
  if (!fFromRootFile) {
    std::vector<double> times;
    std::vector<double> amps;
    double baseline = 0.0;
    double adjustedADC = 0.0;
    int firstSample = 0;
    for (board = 0; board < kNDigitizers; ++board) {
      channel = kNChannels - 1; // only look at beam trigger channel
      int key = PMTInfoMap::CreateKey(board,channel);
      baseline = 0.0;
      for (sample = 0; sample < 50; ++sample) {
        baseline += localCopy.GetSample(key,sample);
      }
      baseline /= 50.0;

      firstSample = 0;
      for (sample = 0; sample < kNSamples; ++sample) {
        adjustedADC = baseline - static_cast<double>(localCopy.GetSample(key,sample));
        if (adjustedADC > 400) {
          firstSample = sample;
          break;
        }
      } // end for sample < kNSamples

      fTriggerTime.at(board) = firstSample;

    } // end for board < kNDigitizers
  } else {
    for (board = 0; board < kNDigitizers; ++board) {
      fTriggerTime.at(board) = localCopy.GetOffset(board);
    }
  }// end if !fFromRootFile

  // skips event if trigger time did not change (no channel 15 NIM signal)
  minTriggerTime = *std::min_element(fTriggerTime.begin(),fTriggerTime.end());
  maxTriggerTime = *std::max_element(fTriggerTime.begin(),fTriggerTime.end());
  if (minTriggerTime == 9e9 || maxTriggerTime == 9e9) {
    MsgFatal("No min and/or max trigger time found. Something went wrong");
  }

  // Get the computer time of the trigger and keep track of the first and last trigger time
  // this is used to fill the data base for nearline monitoring
  std::chrono::time_point <std::chrono::system_clock,std::chrono::duration<unsigned int>> tp_seconds (
      std::chrono::duration<unsigned int>(localCopy.GetGPSSecIntoDay()));
  std::chrono::system_clock::time_point tp (tp_seconds);
  if (fFirstTriggerTime == 0) {
    fFirstTriggerTime = std::chrono::system_clock::to_time_t(tp);
  }

  fLastTriggerTime = std::chrono::system_clock::to_time_t(tp);

  // If the event had mismatched triggers do not find pulses
  if (!evNumCheck) {
    return kCCMDoNotWrite;
  }

  /////////////////////////////////////////////////////////////////
  // Start finding pulses
  /////////////////////////////////////////////////////////////////

  // clear pulses and raw data for the new trigger window
  if (fResetPulses) {
    fPulses->Reset();
  }
  fRawData->operator=(localCopy);

  fPulses->SetEventNumber(localCopy.GetEventNumber());
  fPulses->SetComputerSecIntoEpoch(localCopy.GetGPSSecIntoDay());
  fPulses->SetComputerNSIntoSec(localCopy.GetGPSNSIntoSec());


  double offset = 0.0;
  double pmtOffset = 0.0;

  if (fTruncateWaveform) {
    fRawData->TruncateWaveform(1);
  }

  const size_t kBoardOffset = (kNEffDigitizers == kNDigitizers)? 0 : 10;

  // Loop through each board
  for (board = 0; board < kNEffDigitizers; ++board) {
    offset = 0.0;

    // determine the board time offset
    offset = -(fTriggerTime[board + kBoardOffset] - minTriggerTime);

    // loop through the channels
    for (channel = 0; channel < kNChannels; ++channel) {

      auto key = PMTInfoMap::CreateKey(board+kBoardOffset,channel);

      //fHighestTemp = std::max(fHighestTemp,localCopy.GetTemp(board+kBoardOffset,channel));

      auto samples = localCopy.GetSamples(PMTInfoMap::CreateKey(board,channel));

      // if fTruncateWaveform = true
      // only save fRawData information for the last board
      // this board contains the trigger decode information
      // and the waveforms from the EJ-301 detectors
      //
      // otherwise all the waveforms have already been copied 
      // tot he fRawData object above
      if (board == kNDigitizers-1 && fTruncateWaveform) {
        fRawData->SetWaveform(channel,&samples.front());
      }

      // Check to see if the board/channel combo is turned off in the
      // data analysis.
      // If so, skip the channel.
      if (!PMTInfoMap::IsActive(key)) {
        // Get the #PMTInformation for the current board,channel
        auto pmtInfo = PMTInfoMap::GetPMTInfo(key);
        // If not pmtInfo exists check to see of the 1in PMT time offset needs to be applied
        if (pmtInfo != nullptr) {
          // Check to see if the current pmt is a 1 in pmt.
          // If so, set pmtOffset to 21.0 becuse the 1 in pmts
          // are 42 ns faster than the 8 in pmts based off the 
          // data sheets provided by Hamamatsu 
          if (pmtInfo->Is1in()) {
            pmtOffset = 21.0;
          } else {
            pmtOffset = 0.0;
          }
        }
      }

      // Tell the pulses object which board and channel is currently being looked at
      fPulses->SetBoardChannel(board+kBoardOffset,channel);

      // Apply the #Pulses::DerivativeFilter to the pmt waveform at hand to find the pulses in the waveform
      // Parameters that are passed
      // - \p &samples.front() - The waveform for the current PMT
      // - \p NSAMPELS - The number of samples in the current waveform
      // - \p minTriggerTime - The minimum trigger time for the event (this gets saved for each waveform)
      // - 3 is the number of points to use in the derivative (I believe this variable is not longer being used)
      // - 0 is the threshold applied to the potential pulse that is found (0 means no threshold)
      // - \p offset is the time offset of the board
      // - \p pmtOffset is the 1" to 8" pmt offset
      // - 1.0 is the ADC to PE calibration to be applied (1.0 means no calibration)
      fPulses->DerivativeFilter(&samples.front(),kNSamples,minTriggerTime,3,0,offset,pmtOffset,1.0);

    } // end for channel
  } // end for board

  // Get the number of pulses found over all PMTs and print some diagnostic information
  size_t numPulses = fPulses->GetNumPulses();
  for (size_t pulse = 0; pulse < numPulses; ++pulse) {
    if (fPulses->GetPulseTime(pulse) < 0) {
      MsgInfo(MsgLog::Form("event %d pulse %zu time %g",fNEventsTotal,pulse,fPulses->GetPulseTime(pulse)));
    }
  }

  // Reset the \p fTriggerTime
  fTriggerTime.assign(kNDigitizers,9e9);

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMFindPulses::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("TriggerType").Get(fTriggerType);
  c("TruncateWaveform").Get(fTruncateWaveform);
  c("FromRootFile").Get(fFromRootFile);
  c("ResetPulses").Get(fResetPulses);
  c("WriteDBEntry").Get(fWriteDBEntry);
  if (fWriteDBEntry) {
    c("DBHost").Get(fDBHost);
    c("DBUser").Get(fDBUser);
    c("DBPwd").Get(fDBPwd);
  }

  fIsInit = true;

}

//_______________________________________________________________________________________
CCMResult_t CCMFindPulses::EndOfJob() 
{ 
  WriteDB();

  return kCCMSuccess;
}

//_______________________________________________________________________________________
CCMResult_t CCMFindPulses::NewRun(uint32_t run, uint32_t subRun)
{
  if(fCurrentRun == run && fCurrentSubRun == subRun) {
    return kCCMSuccess;
  }

  if (fNEventsTotal == 0) {
    fNEventsSkipped = 0;
    fHighestTemp = 0;

    fCurrentRun = run;
    fCurrentSubRun = subRun;
    return kCCMSuccess;
  }

  WriteDB();

  fNEventsTotal = 0;
  fNEventsSkipped = 0;
  fHighestTemp = 0;
  
  fCurrentRun = run;
  fCurrentSubRun = subRun;

  return kCCMSuccess;

}

//_______________________________________________________________________________________
void CCMFindPulses::WriteDB()
{
  if (!fWriteDBEntry) {
    return;
  }

  // Write information to data base (code will crash if no connection to the database is made)
  MsgInfo(MsgLog::Form("Number of triggers looped over = %ld, Number of triggers skipped = %ld",fNEventsTotal,fNEventsSkipped));
  struct tm * timeInfo = localtime(&fFirstTriggerTime);
  char bufferFirstTriggerTime[80];
  strftime(bufferFirstTriggerTime,80,"%F %T",timeInfo);
  timeInfo = localtime(&fLastTriggerTime);
  char bufferLastTriggerTime[80];
  strftime(bufferLastTriggerTime,80,"%F %T",timeInfo);

  // defaults
  // - Host = mysql://neutrondaq.lanl.gov/ccmdb
  // - User = ccmdaq
  // - pwd = ask someone
  TSQLServer *db = TSQLServer::Connect(fDBHost.c_str(),fDBUser.c_str(),fDBPwd.c_str());
  std::string insertCommand = "insert into ccmdb.runinfo (start_time,end_time,run_number,subrun_number,num_beam_triggers,num_misaligned_triggers,fHighestTemp) value (";
  MsgInfo(MsgLog::Form("insert command = %s",insertCommand.c_str()));
  insertCommand  += "'" + std::string(bufferFirstTriggerTime) + "','" + 
    std::string(bufferLastTriggerTime) + "','" + 
    fCurrentRun + "','" + 
    fCurrentSubRun + "','"  + 
    std::to_string(fNEventsTotal) + "','" + 
    std::to_string(fNEventsSkipped) + "','" + 
    std::to_string(fHighestTemp) + "');";

  TSQLResult *res = db->Query(insertCommand.c_str());
  res->Print("v");


  delete res;
  delete db;

}

