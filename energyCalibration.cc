/*!**********************************************
 * \file energyCalibration.cc
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
#include "MsgLog.h"
#include "Pulses.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "Utility.h"
#include "RawData.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"

#include "TSystem.h"
#include <cstdlib>
#include <cstdio>


#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>

#include "data_structures.hh"

extern "C" {
#include <unistd.h>
#include <sys/time.h>
#include "getopt.h"
}

/*!**********************************************
 * \fn void Usage()
 * Prints the Usage statement for the current executable
 ***********************************************/
void Usage() 
{

  std::stringstream* usage = new std::stringstream;
  (*usage) << "Usage: energyCalibration [options] \n"
    << "options are:\n"
    << "  -i file.bin           : Add input data file\n"
    << "  -o file.root          : Set output data file\n"
    << "  -t runType            : Set run type name (beam or strobe)\n";

  MsgError((char*)(usage->str().c_str()));
  delete usage;
  return;
}

/*!**********************************************
 * \fn int main(int argc, char ** argv)
 * \brief The main function for energyCalibration
 * \param[in] argc The number of words passed in the command line
 * \param[in] argv An array of the words that are passed
 ***********************************************/
int main(int argc, char ** argv) 
{

  // Used in order to save a vector type into a TTree
  // may not be needed for newer versions of ROOT and g++
  gROOT->ProcessLine("#include <vector>");

  // Declare an #event_t struct to save the data that is read in
  event_t read_data;

  // Objects to save input variables to
  std::string inFileName = "test.bin";
  std::string outFileName = "test.root";
  std::string runType = "default";

  // Declare input parameter options
  static const int kInputOpt    = 'i';
  static const int kOutputOpt   = 'o';
  static const int kRunTypeOpt  = 't';
  static const int kHelpOpt     = 'h';
  static struct option long_options[] = {
    {"input",        required_argument, 0, kInputOpt},
    {"output",       required_argument, 0, kOutputOpt},
    {"runType",      required_argument, 0, kRunTypeOpt},
    {"help",         no_argument,       0, kHelpOpt},
    { NULL, 0, 0, 0} // This is a filler for -1
  };

  // Loop through the input words and pick out the options
  // Save the option to the corresponding input variable
  while (1) {
    int c;
    int optindx = 0;
    c = getopt_long(argc, argv, "c:i:o:t:h", long_options, &optindx);

    if (c==-1) break;
    std::string fname;
    switch (c) {
      case kInputOpt:     inFileName = std::string(optarg); break;
      case kOutputOpt:    outFileName = std::string(optarg); break;
      case kRunTypeOpt:   runType = std::string(optarg); break;
      case kHelpOpt:      Usage(); exit(0); break;
      default:            MsgError(MsgLog::Form("Unknown option %d %s",optind,argv[optind]));
                          Usage();
                          exit(1);
    }
  }

  // Define various variables used in the loops below
  unsigned long totalEventNum = 0;
  unsigned long totalSkippedEvents = 0;

  UInt_t channelEvtNum    = 0;
  UShort_t adcValue       = 0;
  bool above              = false;
  double adjustedADC       = 0.0;
  double baseline          = 0.0;
  double threshold         = 0.0;//static_cast<double>(0.5)*1e-3*std::pow(2.0,14.0)/2.0;
  int board               = 0;
  int channel             = 0;
  int startSample         = 0;
  int endSample           = 0;
  int lastSample          = 0;
  int sample              = 0;
  int sample2             = 0;

  double pulseTime = 0.0;
  double pulseIntegral = 0.0;
  int pulseLength = 0;

  // \p triggerTime stores the time each board saw the trigger that is saved in channel 15
  std::vector<float> triggerTime(NDIGITIZERS,9e9);
  Pulses* pulses = new Pulses(0,0);
  RawData * rawData = new RawData(NDIGITIZERS,NCHANNELS,NSAMPLES,0,0,0);

  unsigned short highestTemp = 0;

  std::time_t firstTriggerTime = 0;//std::chrono::system_clock::to_time_t(tp);
  std::time_t lastTriggerTime = 0;//std::chrono::system_clock::to_time_t(tp);
  unsigned int firstTriggerNS = 0;
  unsigned int lastTriggerNS = 0;

  // Parse the file name to get the \p runNumber, and the \p subRunNumber
  std::string fileNameForParse = inFileName;
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

  MsgInfo(MsgLog::Form("[%s] Looking at file : %s",Utility::tstamp(), inFileName.c_str()));

  // Open the output file and define the trees to be saved
  outROOTFile = new TFile(outFileName.c_str(),"RECREATE","CCM_calibration",1);
  outROOTFile->cd();

  if (gROOT->FindObject("pulses")) {
    delete gROOT->FindObject("pulses");
  }
  pulsesTree = new TTree("pulses","pulses");
  pulsesTree->Branch("pulses",&pulses);
  if (gROOT->FindObject("rawData")) {
    delete gROOT->FindObject("rawData");
  }
  rawTree = new TTree("rawData","rawData");
  rawTree->Branch("rawData",&rawData);

  // The object that reads in the binary file
  std::ifstream binInFile(inFileName.c_str(), std::ios::binary);

  // Loop through the input binary file one event at a time
  while (binInFile.read((char*)&read_data, sizeof(event_t))) {

    if (binInFile.eof()) {
      MsgError("Reached end of binary file");
      break;
    }

    if (binInFile.fail()) {
      MsgFatal(MsgLog::Form("%s did not read successfully",inFileName.c_str()));
    }

    ++totalEventNum;

    // check to make sure all the board event numbers are the same
    int channelEvtNum0 = read_data.digitizers.evNum[0];
    bool evNumCheck = true;
    for (board = 0; board < NDIGITIZERS; ++board) {
      int channelEvtNum = read_data.digitizers.evNum[board];
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
    double channelClockTime = read_data.digitizers.time[0]*8e-9;
    for (board = 1; board < NDIGITIZERS; ++board) {
      //MsgInfo(MsgLog::Form("board %d channelClockTime %g",board,read_data.digitizers.time[board]*8e-9));
      channelClockTime = std::min(channelClockTime,static_cast<double>(read_data.digitizers.time[board]*8e-9));
    }

    double minTriggerTime = 9e9;
    double maxTriggerTime = 9e9;
    if (runType.find("strobe") == std::string::npos) {
      // strobe does not contain the trigger waveform, so have to align based of the reported clock time
      // finds when channel 15 NIM single occurs for each board ("board offset")
      std::vector<double> times;
      std::vector<double> amps;
      double endBaseline = 0.0;
      double baseline = 0.0;
      double adjustedADC = 0.0;
      int firstSample = 0;
      double par0 = 0.0;
      double par1 = 0.0;
      for (board = 0; board < NDIGITIZERS; ++board) {
        channel = NCHANNELS - 1; // only look at beam trigger channel
        baseline = 0.0;
        for (sample = 0; sample < 50; ++sample) {
          baseline += read_data.digitizers.samples[board][channel][sample];
        }
        baseline /= 50.0;

        firstSample = 0;
        for (sample = 0; sample < NSAMPLES; ++sample) {
          adjustedADC = baseline - static_cast<double>(read_data.digitizers.samples[board][channel][sample]);
          if (adjustedADC > 400) {
            firstSample = sample;
            break;
          }
        } // end for sample < NSAMPLES

        triggerTime.at(board) = firstSample;

      } // end for board < NDIGITIZERS

      // skips event if trigger time did not change (no channel 15 NIM signal)
      minTriggerTime = *std::min_element(triggerTime.begin(),triggerTime.end());
      maxTriggerTime = *std::max_element(triggerTime.begin(),triggerTime.end());
      if (minTriggerTime == 9e9 || maxTriggerTime == 9e9) {
        MsgWarning("Should not be here");
        continue;
      }

    } // end if runType does not contain strobe

    // Get the computer time of the trigger and keep track of the first and last trigger time
    // this is used to fill the data base for nearline monitoring
    std::chrono::time_point <std::chrono::system_clock,std::chrono::duration<unsigned int>> tp_seconds (std::chrono::duration<unsigned int>(read_data.computerTime.tv_sec));
    std::chrono::system_clock::time_point tp (tp_seconds);
    if (firstTriggerTime == 0) {
      firstTriggerTime = std::chrono::system_clock::to_time_t(tp);
      firstTriggerNS = read_data.computerTime.tv_nsec;
    }

    lastTriggerTime = std::chrono::system_clock::to_time_t(tp);
    lastTriggerNS = read_data.computerTime.tv_nsec;

    // If the event had mismatched triggers do not find pulses
    if (!evNumCheck) {
      continue;
    }

    /////////////////////////////////////////////////////////////////
    // Start finding pulses
    /////////////////////////////////////////////////////////////////

    // clear pulses and raw data for the new trigger window
    pulses->Reset();
    rawData->Reset(NDIGITIZERS,NCHANNELS,NSAMPLES,read_data.evNum,read_data.computerTime.tv_sec,read_data.computerTime.tv_nsec);

    // Set meta data information for both pulses and raw data
    rawData->SetOffset(triggerTime);
    rawData->SetBoardEventNum(read_data.digitizers.evNum);
    rawData->SetClockTime(read_data.digitizers.time);

    pulses->SetEventNumber(read_data.evNum);
    pulses->SetComputerSecIntoEpoch(read_data.computerTime.tv_sec);
    pulses->SetComputerNSIntoSec(read_data.computerTime.tv_nsec);


    double offset = 0.0;
    double pmtOffset = 0.0;
    
    // Loop through each board
    for (board = 0; board < NDIGITIZERS; ++board) {
      offset = 0.0;
      // fill rawData information and calculate the highest temperature
      rawData->SetChannelMask(board,read_data.digitizers.chMask[board]);
      rawData->SetSize(board,read_data.digitizers.size[board]);
      highestTemp = std::max(highestTemp,rawData->SetTemp(board,read_data.digitizers.temperatures[board]));

      // determine the board time offset
      if (runType.find("strobe") == std::string::npos) {
        offset = -(triggerTime[board] - minTriggerTime);
      } else {
        offset = -std::floor(static_cast<double>(read_data.digitizers.time[board]*8e-9) - channelClockTime/2.0);
        minTriggerTime = 0;
      }

      // loop through the channels
      for (channel = 0; channel < NCHANNELS; ++channel) {

        // only save rawData information for the last board
        // this board contains the trigger decode information
        // and the waveforms from the EJ-301 detectors
        if (board == NDIGITIZERS-1) {
          rawData->SetWaveform(channel,read_data.digitizers.samples[board][channel]);
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
        pulses->SetBoardChannel(board,channel);

        // Apply the #Pulses::DerivativeFilter to the pmt waveform at hand to find the pulses in the waveform
        // Parameters that are passed
        // - \p read_data.digitizers.samples - The waveform for the current PMT
        // - \p NSAMPELS - The number of samples in the current waveform
        // - \p minTriggerTime - The minimum trigger time for the event (this gets saved for each waveform)
        // - 3 is the number of points to use in the derivative (I believe this variable is not longer being used)
        // - 0 is the threshold applied to the potential pulse that is found (0 means no threshold)
        // - \p offset is the time offset of the board
        // - \p pmtOffset is the 1" to 8" pmt offset
        // - 1.0 is the ADC to PE calibration to be applied (1.0 means no calibration)
        pulses->DerivativeFilter(read_data.digitizers.samples[board][channel],NSAMPLES,minTriggerTime,3,0,offset,pmtOffset,1.0);

      } // end for channel
    } // end for board

    // Get the number of pulses found over all PMTs and print some diagnostic information
    size_t numPulses = pulses->GetNumPulses();
    for (size_t pulse = 0; pulse < numPulses; ++pulse) {
      if (pulses->GetPulseTime(pulse) < 0) {
        MsgInfo(MsgLog::Form("event %d pulse %zu time %g",totalEventNum,pulse,pulses->GetPulseTime(pulse)));
      }
    }

    // Save the \p rawTree and \p pulses information to file
    outROOTFile->cd();
    rawTree->Fill();
    if (pulses->GetNumPulses()) {
      pulsesTree->Fill();
    }

    // Reset the \p triggerTime
    triggerTime.assign(NDIGITIZERS,9e9);

  } // end while binInFile.good

  // Close input and output files
  binInFile.close();

  outROOTFile->cd();
  //outROOTFile->Flush();
  rawTree->Write();
  pulsesTree->Write();
  delete outROOTFile;

  // Write information to data base (code will crash if no connection to the database is made)
  MsgInfo(MsgLog::Form("[%s] Number of triggers looped over = %ld, Number of triggers skipped = %ld",Utility::tstamp(), totalEventNum,totalSkippedEvents));
  struct tm * timeInfo = localtime(&firstTriggerTime);
  char bufferFirstTriggerTime[80];
  strftime(bufferFirstTriggerTime,80,"%F %T",timeInfo);
  timeInfo = localtime(&lastTriggerTime);
  char bufferLastTriggerTime[80];
  strftime(bufferLastTriggerTime,80,"%F %T",timeInfo);

  TSQLServer *db = TSQLServer::Connect("mysql://neutrondaq.lanl.gov/ccmdb","ccmdaq", "ccmdbpwd");
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

  return EXIT_SUCCESS;

}

