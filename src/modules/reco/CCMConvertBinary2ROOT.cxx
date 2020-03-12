/*!**********************************************
 * \file CCMConvertBinary2ROOT
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
#include "CCMConvertBinary2ROOT.h"
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
MODULE_DECL(CCMConvertBinary2ROOT);

//_______________________________________________________________________________________
CCMConvertBinary2ROOT::CCMConvertBinary2ROOT(const char* version) 
  : CCMModule("CCMConvertBinary2ROOT"),
    fInFileName(""),
    fOutFileName(""),
    fWriteDBEntry(false),
    fDBHost(""),
    fDBUser(""),
    fDBPwd("")
{
  //Default constructor
  this->SetCfgVersion(version);

  fTriggerTime.resize(NDIGITIZERS,9e9),
  fPulses = new Pulses(0,0);
  fRawData = new RawData(NDIGITIZERS,NCHANNELS,NSAMPLES,0,0,0);
}

//_______________________________________________________________________________________
CCMConvertBinary2ROOT::CCMConvertBinary2ROOT(const CCMConvertBinary2ROOT& clufdr) 
: CCMModule(clufdr),
  fInFileName(clufdr.fInFileName),
  fOutFileName(clufdr.fOutFileName),
  fWriteDBEntry(clufdr.fWriteDBEntry),
  fDBHost(clufdr.fDBHost),
  fDBUser(clufdr.fDBUser),
  fDBPwd(clufdr.fDBPwd),
  fTriggerTime(clufdr.fTriggerTime),
  fPulses(clufdr.fPulses),
  fRawData(clufdr.fRawData)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMConvertBinary2ROOT::~CCMConvertBinary2ROOT()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMConvertBinary2ROOT::ProcessEvent()
{

  // Define various variables used in the loops below
  unsigned long totalEventNum = 0;
  unsigned long totalSkippedEvents = 0;

  int board               = 0;
  int channel             = 0;
  int sample              = 0;

  // \p fTriggerTime stores the time each board saw the trigger that is saved in channel 15
  fTriggerTime.assign(NDIGITIZERS,9e9);

  unsigned short highestTemp = 0;

  std::time_t firstTriggerTime = 0;//std::chrono::system_clock::to_time_t(tp);
  std::time_t lastTriggerTime = 0;//std::chrono::system_clock::to_time_t(tp);
  //unsigned int firstTriggerNS = 0;
  //unsigned int lastTriggerNS = 0;

  // Parse the file name to get the \p runNumber, and the \p subRunNumber
  std::string fileNameForParse = fInFileName;
  size_t locOfRun = fileNameForParse.find("run");
  std::string runNumber = "";
  std::string subRunNumber = "";
  runNumber.assign(fileNameForParse,locOfRun+3,6);
  subRunNumber.assign(fileNameForParse,locOfRun+3+6+1,6);
  MsgInfo(MsgLog::Form("locOfRun = %zu RunNumber is %s SubRunNumber is %s",locOfRun,runNumber.c_str(),subRunNumber.c_str()));


  // Load the default PMT map from csv files
  // TODO: Need to switch to reading from sql map
  // but may not be doable with running on various
  // networks
  PMTInfoMap::DefaultFillPMTMap();

  // Declare the TFile and TTrees that are saved
  TFile * outROOTFile = 0;//= TFile::Open(outfileName.c_str(),"RECREATE","CCM_calibration",1);
  TTree * pulsesTree = 0;//= new TTree("pulses","pulses");
  TTree * rawTree = 0;

  MsgInfo(MsgLog::Form("[%s] Looking at file : %s",Utility::tstamp(), fInFileName.c_str()));

  // Open the output file and define the trees to be saved
  outROOTFile = new TFile(fOutFileName.c_str(),"RECREATE","CCM_calibration",1);
  outROOTFile->cd();

  if (gROOT->FindObject("pulses")) {
    delete gROOT->FindObject("pulses");
  }
  pulsesTree = new TTree("pulses","pulses");
  pulsesTree->Branch("pulses",&fPulses);
  if (gROOT->FindObject("rawData")) {
    delete gROOT->FindObject("rawData");
  }
  rawTree = new TTree("rawData","rawData");
  rawTree->Branch("rawData",&fRawData);

  // The object that reads in the binary file
  std::ifstream binInFile(fInFileName.c_str(), std::ios::binary);

  // Loop through the input binary file one event at a time
  while (binInFile.read((char*)&fReadData, sizeof(event_t))) {

    if (binInFile.eof()) {
      MsgError("Reached end of binary file");
      break;
    }

    if (binInFile.fail()) {
      MsgFatal(MsgLog::Form("%s did not read successfully",fInFileName.c_str()));
    }

    ++totalEventNum;

    // check to make sure all the board event numbers are the same
    int channelEvtNum0 = fReadData.digitizers.evNum[0];
    bool evNumCheck = true;
    for (board = 0; board < NDIGITIZERS; ++board) {
      int channelEvtNum = fReadData.digitizers.evNum[board];
      if (channelEvtNum != channelEvtNum0) {
        evNumCheck = false;
      }
    }

    // count number of times we have a mismatched trigger
    if (!evNumCheck) {
      MsgError(MsgLog::Form("Mismatched Trigger %ld",totalEventNum));
      ++totalSkippedEvents;
      //continue; // do not skip because of information for the database
    }


    // monitor output
    if (totalEventNum % 100 == 0) {
      MsgInfo(MsgLog::Form("[%s] At trigger %ld",Utility::tstamp(), totalEventNum));
    }

    /////////////////////////////////////////////////////////////////
    // Calculate the trigger time of each board in the DAQ window
    // Use the earliest trigger time as the time of the trigger
    // and shift all the other boards to line up with the earliest board
    /////////////////////////////////////////////////////////////////

    // used if not saving trigger to channel 15
    double channelClockTime = fReadData.digitizers.time[0]*8e-9;
    for (board = 1; board < NDIGITIZERS; ++board) {
      //MsgInfo(MsgLog::Form("board %d channelClockTime %g",board,fReadData.digitizers.time[board]*8e-9));
      channelClockTime = std::min(channelClockTime,static_cast<double>(fReadData.digitizers.time[board]*8e-9));
    }

    double minTriggerTime = 9e9;
    double maxTriggerTime = 9e9;
    std::vector<double> times;
    std::vector<double> amps;
    double baseline = 0.0;
    double adjustedADC = 0.0;
    int firstSample = 0;
    for (board = 0; board < NDIGITIZERS; ++board) {
      channel = NCHANNELS - 1; // only look at beam trigger channel
      baseline = 0.0;
      for (sample = 0; sample < 50; ++sample) {
        baseline += fReadData.digitizers.samples[board][channel][sample];
      }
      baseline /= 50.0;

      firstSample = 0;
      for (sample = 0; sample < NSAMPLES; ++sample) {
        adjustedADC = baseline - static_cast<double>(fReadData.digitizers.samples[board][channel][sample]);
        if (adjustedADC > 400) {
          firstSample = sample;
          break;
        }
      } // end for sample < NSAMPLES

      fTriggerTime.at(board) = firstSample;

    } // end for board < NDIGITIZERS

    // skips event if trigger time did not change (no channel 15 NIM signal)
    minTriggerTime = *std::min_element(fTriggerTime.begin(),fTriggerTime.end());
    maxTriggerTime = *std::max_element(fTriggerTime.begin(),fTriggerTime.end());
    if (minTriggerTime == 9e9 || maxTriggerTime == 9e9) {
      MsgWarning("Should not be here");
      continue;
    }

    // Get the computer time of the trigger and keep track of the first and last trigger time
    // this is used to fill the data base for nearline monitoring
    std::chrono::time_point <std::chrono::system_clock,std::chrono::duration<unsigned int>> tp_seconds (std::chrono::duration<unsigned int>(fReadData.computerTime.tv_sec));
    std::chrono::system_clock::time_point tp (tp_seconds);
    if (firstTriggerTime == 0) {
      firstTriggerTime = std::chrono::system_clock::to_time_t(tp);
      //firstTriggerNS = fReadData.computerTime.tv_nsec;
    }

    lastTriggerTime = std::chrono::system_clock::to_time_t(tp);
    //lastTriggerNS = fReadData.computerTime.tv_nsec;

    // If the event had mismatched triggers do not find pulses
    if (!evNumCheck) {
      continue;
    }

    /////////////////////////////////////////////////////////////////
    // Start finding pulses
    /////////////////////////////////////////////////////////////////

    // clear pulses and raw data for the new trigger window
    fPulses->Reset();
    fRawData->Reset(NDIGITIZERS,NCHANNELS,NSAMPLES,fReadData.evNum,fReadData.computerTime.tv_sec,fReadData.computerTime.tv_nsec);

    // Set meta data information for both pulses and raw data
    fRawData->SetOffset(fTriggerTime);
    fRawData->SetBoardEventNum(fReadData.digitizers.evNum);
    fRawData->SetClockTime(fReadData.digitizers.time);

    fPulses->SetEventNumber(fReadData.evNum);
    fPulses->SetComputerSecIntoEpoch(fReadData.computerTime.tv_sec);
    fPulses->SetComputerNSIntoSec(fReadData.computerTime.tv_nsec);


    double offset = 0.0;
    double pmtOffset = 0.0;
    
    // Loop through each board
    for (board = 0; board < NDIGITIZERS; ++board) {
      offset = 0.0;
      // fill fRawData information and calculate the highest temperature
      fRawData->SetChannelMask(board,fReadData.digitizers.chMask[board]);
      fRawData->SetSize(board,fReadData.digitizers.size[board]);
      highestTemp = std::max(highestTemp,fRawData->SetTemp(board,fReadData.digitizers.temperatures[board]));

      // determine the board time offset
      offset = -(fTriggerTime[board] - minTriggerTime);

      // loop through the channels
      for (channel = 0; channel < NCHANNELS; ++channel) {

        // only save fRawData information for the last board
        // this board contains the trigger decode information
        // and the waveforms from the EJ-301 detectors
        if (board == NDIGITIZERS-1) {
          fRawData->SetWaveform(channel,fReadData.digitizers.samples[board][channel]);
        }

        // Check to see if the board/channel combo is turned off in the
        // data analysis.
        // If so, skip the channel.
        if (!PMTInfoMap::IsActive(board,channel)) {
          continue;
        }

        // Get the #PMTInformation for the current board,channel
        const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(board,channel);
        // If not \p pmtInfo exists, move to the next channel
        if (pmtInfo == nullptr) {
          continue;
        }

        // Check to see if the current pmt is a 1 in pmt.
        // If so, set \p pmtOffset to 21.0 becuse the 1 in pmts
        // are 42 ns faster than the 8 in pmts based off the 
        // data sheets provided by Hamamatsu 
        if (pmtInfo->Is1in()) {
          pmtOffset = 21.0;
        } else {
          pmtOffset = 0.0;
        }

        // Tell the pulses object which board and channel is currently being looked at
        fPulses->SetBoardChannel(board,channel);

        // Apply the #Pulses::DerivativeFilter to the pmt waveform at hand to find the pulses in the waveform
        // Parameters that are passed
        // - \p fReadData.digitizers.samples - The waveform for the current PMT
        // - \p NSAMPELS - The number of samples in the current waveform
        // - \p minTriggerTime - The minimum trigger time for the event (this gets saved for each waveform)
        // - 3 is the number of points to use in the derivative (I believe this variable is not longer being used)
        // - 0 is the threshold applied to the potential pulse that is found (0 means no threshold)
        // - \p offset is the time offset of the board
        // - \p pmtOffset is the 1" to 8" pmt offset
        // - 1.0 is the ADC to PE calibration to be applied (1.0 means no calibration)
        fPulses->DerivativeFilter(fReadData.digitizers.samples[board][channel],NSAMPLES,minTriggerTime,3,0,offset,pmtOffset,1.0);

      } // end for channel
    } // end for board

    // Get the number of pulses found over all PMTs and print some diagnostic information
    size_t numPulses = fPulses->GetNumPulses();
    for (size_t pulse = 0; pulse < numPulses; ++pulse) {
      if (fPulses->GetPulseTime(pulse) < 0) {
        MsgInfo(MsgLog::Form("event %d pulse %zu time %g",totalEventNum,pulse,fPulses->GetPulseTime(pulse)));
      }
    }

    // Save the \p rawTree and \p pulses information to file
    outROOTFile->cd();
    rawTree->Fill();
    if (fPulses->GetNumPulses()) {
      pulsesTree->Fill();
    }

    // Reset the \p fTriggerTime
    fTriggerTime.assign(NDIGITIZERS,9e9);

  } // end while binInFile.good

  // Close input and output files
  binInFile.close();

  outROOTFile->cd();
  //outROOTFile->Flush();
  rawTree->Write();
  pulsesTree->Write();
  delete outROOTFile;

  if (fWriteDBEntry) {

    // Write information to data base (code will crash if no connection to the database is made)
    MsgInfo(MsgLog::Form("[%s] Number of triggers looped over = %ld, Number of triggers skipped = %ld",Utility::tstamp(), totalEventNum,totalSkippedEvents));
    struct tm * timeInfo = localtime(&firstTriggerTime);
    char bufferFirstTriggerTime[80];
    strftime(bufferFirstTriggerTime,80,"%F %T",timeInfo);
    timeInfo = localtime(&lastTriggerTime);
    char bufferLastTriggerTime[80];
    strftime(bufferLastTriggerTime,80,"%F %T",timeInfo);

    // defaults
    // - Host = mysql://neutrondaq.lanl.gov/ccmdb
    // - User = ccmdaq
    // - pwd = ask someone
    TSQLServer *db = TSQLServer::Connect(fDBHost.c_str(),fDBUser.c_str(),fDBPwd.c_str());
    std::string insertCommand = "insert into ccmdb.runinfo (start_time,end_time,run_number,subrun_number,num_beam_triggers,num_misaligned_triggers,highestTemp) value (";
    MsgInfo(MsgLog::Form("insert command = %s",insertCommand.c_str()));
    insertCommand  += "'" + std::string(bufferFirstTriggerTime) + "','" + 
      std::string(bufferLastTriggerTime) + "','" + 
      runNumber + "','" + 
      subRunNumber + "','"  + 
      std::to_string(totalEventNum) + "','" + 
      std::to_string(totalSkippedEvents) + "','" + 
      std::to_string(highestTemp) + "');";
    TSQLResult *res = db->Query(insertCommand.c_str());
    res->Print("v");

    delete res;
    delete db;
  } // end if fWriteDBEntry

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMConvertBinary2ROOT::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  if (&c != 0)
  {
    c("InFileName").Get(fInFileName);
    c("OutFileName").Get(fOutFileName);
    c("WriteDBEntry").Get(fWriteDBEntry);
    c("DBHost").Get(fDBHost);
    c("DBUser").Get(fDBUser);
    c("DBPwd").Get(fDBPwd);

    fIsInit = true;
  } else {
    fIsInit = false;
  }

}

