/*!**********************************************
 * \file nearlineDiag.cc
 * \brief Main function for the NearlineDiag executable
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date February 24, 2020
 * 
 * Main code for the NearlineDiag executable.
 ***********************************************/
#include "MsgLog.h"
#include "Pulses.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "Utility.h"
#include "NearlineSPEDiag.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

/*!**********************************************
 * \fn int main(int argc, char ** argv)
 * \brief The main function for NearlineDiag
 * \param[in] argc The number of words passed in the command line
 * \param[in] argv An array of the words that are passed
 ***********************************************/
int main(int argv, char ** argc) 
{

  gROOT->ProcessLine("#include <vector>");

  std::ifstream fout;

  // list of file names to loop over
  std::string fileListName = argc[1];

  std::ostream_iterator<std::string> out_string(std::cout,"\n");
  std::vector<std::string> fileNames;

  // save the file names to a vector
  fout.open(fileListName.c_str());
  std::string fileName = "";
  while (fout >> fileName) {
    fileNames.push_back(fileName);
  }
  fout.close();
  std::sort(fileNames.begin(),fileNames.end());
  std::copy(fileNames.begin(),fileNames.end(),out_string);

  // load the default pmt map
  PMTInfoMap::DefaultFillPMTMap();

  // load the calibration file (note the calibration used was an old calibration
  // it allowed us to see a consistent trend throughout the run)
  TFile * calibrationFile = TFile::Open("root_out_2019ledData_run118_.root","READ");
  std::map<int,std::pair<float,float>> calibrationValues;
  TTreeReader calibrationTree("spe",calibrationFile);
  TTreeReaderValue<float> speValue(calibrationTree,"speValue");
  TTreeReaderValue<float> speThreshold(calibrationTree,"endNoiseWallFitRangeStart");
  TTreeReaderValue<int> pmtID(calibrationTree,"pmtID");
  while (calibrationTree.Next()) {
    calibrationValues.insert(std::make_pair(*pmtID,std::make_pair(*speValue,*speThreshold)));
  }
  delete calibrationFile;

  // define various values needed for the analysis
  Pulses* pulses = 0;

  size_t overlapTriggers = 0;
  size_t totalTriggers = 0;


  size_t fileCount = 0;

  TFile * inROOTFile = 0;

  TTree * tree = 0;

  std::array<int,800> totalNumPMTs;
  std::array<double,800> totalIntegral;
  std::array<int,8000> numPMTs;
  std::array<double,8000> integral;

  std::array<int,800> totalNumPMTsVeto;
  std::array<double,800> totalIntegralVeto;
  std::array<int,8000> numPMTsVeto;
  std::array<double,8000> integralVeto;

  std::array<int,160> specount; //spe count for each of the 160 pmts in the first 10 boards. use as specount[key] in the pulse loop.       
  std::array<int,160> pmttype;// pmt type (0=1inveto,1=uncoated,2=coated) for the 160 pmts. use as pmttype[key] in pulse loop.              

  totalNumPMTs.fill(0);
  totalIntegral.fill(0);

  double preBeamTriggers = 0;
  double inBeamIntegral = 0.0;
  double preBeamPrior = 0;

  std::time_t firstTriggerTime = 0;//std::chrono::system_clock::to_time_t(tp);
  std::time_t lastTriggerTime = 0;//std::chrono::system_clock::to_time_t(tp);
  unsigned int firstTriggerNS = 0;
  unsigned int lastTriggerNS = 0;

  std::ofstream tfile;
  tfile.open("testratecounts.txt");

  // loop over each file name
  for (const auto & file : fileNames) {

    inROOTFile = TFile::Open(file.c_str(),"READ");
    if (!inROOTFile) {
      MsgFatal(MsgLog::Form("Could not open %s",file.c_str()));
    }

    inROOTFile->GetObject("pulses",tree);
    if (!tree) {
      MsgWarning(MsgLog::Form("Could not find tree %s in file %s","pulses",file.c_str()));
      continue;
    }

    tree->SetBranchAddress("pulses",&pulses);

    MsgInfo(MsgLog::Form("Looking at file : %s",file.c_str()));

    const long kNTriggers = tree->GetEntries();
    totalTriggers += kNTriggers;

    double testnPMTS = 0;
    double testninteg = 0;
    double testnPMTsVeto = 0;
    double testnintegVeto = 0;

    // loop over each trigger in the file
    for (long trigger = 0; trigger < kNTriggers; ++trigger) {

      tree->GetEntry(trigger);

      if (!pulses) {
        MsgWarning(MsgLog::Form("Trigger %ld does not contain Pulses*",trigger));
        continue;
      }

      if (trigger % 100 == 0) {
        MsgInfo(MsgLog::Form("[%s] At event %ld out of %ld (%.4f%%)",Utility::tstamp(),
              trigger,kNTriggers,static_cast<double>(trigger)/static_cast<double>(kNTriggers)*100.0));
      }


      if (pulses->GetTriggerTime()*2e-3 - 9.92 <= -0.5) {
        ++overlapTriggers;
        continue;
      }

      pulses->Sort();

      numPMTs.fill(0);
      integral.fill(0.0);
      numPMTsVeto.fill(0);
      integralVeto.fill(0.0);

      int tempStart = -100;
      int tempEnd = -100;
      double longestPulse = 0;
      double longPulseInt = 0;

      // grab the time of the trigger
      std::chrono::time_point <std::chrono::system_clock,std::chrono::duration<unsigned int>> tp_seconds (std::chrono::duration<unsigned int>(pulses->GetComputerSecIntoEpoch()));
      std::chrono::system_clock::time_point tp (tp_seconds);
      if (firstTriggerTime == 0) {
        firstTriggerTime = std::chrono::system_clock::to_time_t(tp);
        firstTriggerNS = pulses->GetComputerNSIntoSec();
      }

      lastTriggerTime = std::chrono::system_clock::to_time_t(tp);
      lastTriggerNS = pulses->GetComputerNSIntoSec();


      // loop over each pulse in the DAQ window
      const size_t kNumPulses = pulses->GetNumPulses();
      for (size_t loc = 0; loc < kNumPulses; ++loc) {
        int key = pulses->GetKey(loc);
        const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(key);
        if (!pmtInfo) {
          continue;
        }

        double time = pulses->GetPulseTime(loc);

        if (std::isnan(time) || std::isinf(time)) {
          continue;
        }

        double threshold = calibrationValues[key].second;
        double spe = calibrationValues[key].first;

        double pulseIntegral = pulses->GetPulseIntegral(loc);
        if (std::isnan(pulseIntegral) || std::isinf(pulseIntegral)) {
          continue;
        }

        if (pmtInfo->IsVeto()) {
          if (pmtInfo->Is1in()) {
            if (pulseIntegral > 10.0) {
              ++specount[key];
              pmttype[key] = 1;
            }
          } else {
            if (pulseIntegral > 10.0) {
              ++specount[key];
              pmttype[key] = 0;
            }
          }
        } else {
          if (pulseIntegral < threshold || pulseIntegral < 0 || spe <= 0) {
            continue;
          } else if (time*2e-3 - 9.92 <= -1) {
            ++specount[key];
            if (pmtInfo->IsUncoated()){
              pmttype[key]=2;
            }else{
              pmttype[key]=3;
            }
            // if PulseIntegral is above threshold and 0 and the spe if not less than 0, add 1 to the specount for the [key] pmt              
            //then use the pmtInfo functions Is1in and IsUncoated to determine the type of pmt and define that.                               
            //redefined for every pulse, but this should not be an issue. Will not define a type if there are no spes (deemed acceptable)     
          }
        }

        if (pmtInfo->IsVeto()){
          pulseIntegral /= 10.0;
        }else{
          pulseIntegral /= spe;
        }

        double length = pulses->GetPulseLength(loc);

        //if (longPulseInt < integral/ADCToPEDer(pulses->GetBoard(loc),pulses->GetChannel(loc))) {
        //  longestPulse = length;
        //  longPulseInt = integral/ADCToPEDer(pulses->GetBoard(loc),pulses->GetChannel(loc));
        //}
        //longPulseHist->Fill(length,integral/ADCToPEDer(pulses->GetBoard(loc),pulses->GetChannel(loc)));


        int firstBin = time;
        int lastBin  = time+length;
        double middle = (lastBin + firstBin)/2.0;
        if (pmtInfo->IsVeto()) {
          for (int bin = firstBin; bin < lastBin && bin < 8000; ++bin) {
            //if (bin != firstBin) {
            ++numPMTsVeto[bin];
            ++totalNumPMTsVeto[bin/10];
            //}
            if (bin < middle) {
              integralVeto[bin] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(bin - firstBin)/(middle - firstBin));
              totalIntegralVeto[bin/10] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(bin - firstBin)/(middle - firstBin));
            } else if (bin != middle) {
              integralVeto[bin] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(lastBin - bin)/(lastBin - middle));
              totalIntegralVeto[bin/10] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(lastBin - bin)/(lastBin - middle));
            } else {
              integralVeto[bin] += 2.0*pulseIntegral/(lastBin-firstBin);
              totalIntegralVeto[bin/10] += 2.0*pulseIntegral/(lastBin-firstBin);
            }
          }  
          continue;
        }

        for (int bin = firstBin; bin < lastBin && bin < 8000; ++bin) {
          //if (bin != firstBin) {
          ++numPMTs[bin];
          ++totalNumPMTs[bin/10];
          //}
          if (bin < middle) {
            integral[bin] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(bin - firstBin)/(middle - firstBin));
            totalIntegral[bin/10] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(bin - firstBin)/(middle - firstBin));
          } else if (bin != middle) {
            integral[bin] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(lastBin - bin)/(lastBin - middle));
            totalIntegral[bin/10] += 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(lastBin - bin)/(lastBin - middle));
          } else {
            integral[bin] += 2.0*pulseIntegral/(lastBin-firstBin);
            totalIntegral[bin/10] += 2.0*pulseIntegral/(lastBin-firstBin);
          }
        } // end for bin = firstBin
      } // end for size_t loc = 0

      bool below = false;
      bool belowhigher = false;

      for (int bin = 0; bin < 8000; ++bin) {
        //MsgInfo(MsgLog::Form("bin %d hit %d int %.2f",bin,fPassedThresholdCount[bin],fPassedThresholdInt[bin]));
        if (bin < 4460) {
          if (!below && numPMTs[bin] >= 5 && integral[bin] >= 0.2) {
            below = true;
            //            ++preBeamTriggers;
            testnPMTS = numPMTs[bin];
            testninteg = integral[bin];
            testnPMTsVeto = static_cast<double>(numPMTsVeto[bin]);
            testnintegVeto = static_cast<double>(integralVeto[bin]);
            //add one to preBeamTriggers and store the number of PMTs and integral for each event.                                          
          } else if (below) {
            if (integral[bin] > testninteg) {
              testnPMTS = numPMTs[bin];
              testninteg = integral[bin];
            }
            //if the integral is higher ffor a later part of the waveform in the event, replace the stored npmts and integral with the higher value.
            if (static_cast<double>(integralVeto[bin]) > testnintegVeto){
              testnPMTsVeto = static_cast<double>(numPMTsVeto[bin]);
              testnintegVeto = static_cast<double>(integralVeto[bin]);
            }
            // within a given event, find the highest veto integral and change the stored npmts veto and integral veto with that value 
            if (numPMTs[bin] < 5 || integral[bin] < 0.2) {
              //MsgInfo("End of \"Event\"");                             
              below = false;
	      //        lfile << lfilen <<  '\t' << testnPMTS <<  '\t' << testninteg <<   '\t' << testnPMTsVeto <<  '\t' << testnintegVeto << '\n';
            }
          } // end if (numPMTs and integral above threshold

          if (!belowhigher && numPMTs[bin] >=20 && integral[bin] >= 5) {
            belowhigher = true;
            ++preBeamTriggers;
          } else if (belowhigher && numPMTs[bin] < 20 || integral[bin] < 5){
            belowhigher = false;
          }
        } else if (bin > 4460 && bin < 5460) {
          if (!std::isnan(integral[bin]) && !std::isinf(integral[bin])) {
            inBeamIntegral += integral[bin];
          }
        } // end if-else on  bin
      } // end for bin
    } // end for trigger

    //    lfile.close();

    MsgInfo(MsgLog::Form("File Count %d Number of pre beam events %g inBeamIntegral %g",fileCount,preBeamTriggers,inBeamIntegral));

    delete inROOTFile;
    inROOTFile = 0;

    MsgInfo(MsgLog::Form("Increment file count to %d out of %zu",++fileCount,fileNames.size()));

    tfile << file.c_str() << '\t' <<(preBeamTriggers-preBeamPrior) <<'\t' << inBeamIntegral <<'\t' << totalTriggers <<'\t' << overlapTriggers<<'\n';
    preBeamPrior = preBeamTriggers;

  } // end loop over name vector

  tfile.close();

  ///////////////////////////////////////////////
  // Save various information to text files
  // and to the database if connection exists
  ///////////////////////////////////////////////

  std::string specountoutput = "Board\tChannel\tSPEcount\tSPErate\tType\n";
  int digit = 0;
  int channel = 0;
  double sperate = 0.0;

  std::array<int,3> pmtsoftype;
  std::array<int,3> spespertype;
  std::array<std::string,3> pmtstypen;

  for (int pmtc = 0; pmtc<160; ++pmtc){
    //PMTInfoMap::DecodeKey(pmtc,digit,channel);
    //This was giving highly inaccurate values (channels up to 60 for board 9, which did not make any sense.

    digit = pmtc/16;
    channel = pmtc%16;
    if (channel == 15){
      continue;
    }
    if (specount[pmtc] == 0) {
      continue;
    }

    sperate = specount[pmtc]/(8.92*static_cast<double>(totalTriggers))*1000.0;
    spespertype[pmttype[pmtc]-1] += sperate;
    //MsgInfo(MsgLog::Form("Board %d Channel %d SPE count %d SPE rate %g",digit,channel,specount[pmtc],sperate));                           
    specountoutput = specountoutput+Form("%d\t%d\t%d\t%g\t%d\n",digit,channel,specount[pmtc],sperate,pmttype[pmtc]);
    if (pmttype[pmtc]==1){
      ++pmtsoftype[0];
      pmtstypen[0] = "Veto";
    }else if (pmttype[pmtc]==2) {// && sperate > 100.){
      ++pmtsoftype[1];
      pmtstypen[1] = "Uncoated";
    }else if (pmttype[pmtc]==3) {// && sperate > 100.){
      ++pmtsoftype[2];
      pmtstypen[2] = "Coated";
    }
    }

    std::array<float,3> sperates;

    for (int type=0;type<3;++type){
      sperates[type] = spespertype[type]/pmtsoftype[type];
      specountoutput = specountoutput+pmtstypen[type]+Form("\t%d\t%g\t%d\n",spespertype[type],sperates[type],pmtsoftype[type]);
      MsgInfo(MsgLog::Form("%s\t%d\t%g kHz\n",pmtstypen[type].c_str(),pmtsoftype[type],sperates[type]));
    }

    std::ofstream sfile;
    sfile.open("speratestest.txt");
    sfile << specountoutput;
    sfile.close();

    totalTriggers -= overlapTriggers;

    MsgInfo(MsgLog::Form("Number of pre beam events %g inBeamIntegral %g total triggers = %zu overlap triggers = %zu",
          preBeamTriggers,inBeamIntegral,totalTriggers,overlapTriggers));

    //  return EXIT_SUCCESS;

    struct tm * timeInfo = localtime(&firstTriggerTime);
    char bufferFirstTriggerTime[80];
    strftime(bufferFirstTriggerTime,80,"%F %T",timeInfo);
    timeInfo = localtime(&lastTriggerTime);
    char bufferLastTriggerTime[80];
    strftime(bufferLastTriggerTime,80,"%F %T",timeInfo);

    MsgInfo(MsgLog::Form("First Time %s Last Time %s",bufferFirstTriggerTime,bufferLastTriggerTime));

    float prebeamrate = preBeamTriggers/static_cast<double>(totalTriggers)/8.92*1000.;

    std::string insertCommand = "insert into ccmdb.nearline_diag (start_time, end_time, prebeam_event_rate, spe_rate_1in, spe_rate_uncoated, spe_rate_coated) value (";
    insertCommand = Form("%s'%s','%s','%f','%f','%f','%f');",
        insertCommand.c_str(),
        bufferFirstTriggerTime,bufferLastTriggerTime,
        prebeamrate,sperates[0],sperates[1],sperates[2]);

    TSQLServer *db = TSQLServer::Connect("mysql://neutrondaq.lanl.gov/ccmdb","ccmdaq", "ccmdbpwd");
    MsgInfo(MsgLog::Form("Insert Command = %s",insertCommand.c_str()));
    TSQLResult *res = db->Query(insertCommand.c_str());
    delete res;
    /*
       std::string insertCommand = "insert into ccmdb.nearline_event(start_time,end_time,num_pre_beam_events,in_beam_integral) value (";
       insertCommand = Form("%s'%s','%s','%f','%f');",insertCommand.c_str(), bufferFirstTriggerTime,bufferLastTriggerTime,
       preBeamTriggers/static_cast<double>(totalTriggers)/8.92,inBeamIntegral/static_cast<double>(totalTriggers));
    //insertCommand += "'"+std::string(bufferFirstTriggerTime)+"','"+std::string(bufferLastTriggerTime)+"','"+std::to_string(preBeamTriggers)+"','"+std::to_string(inBeamIntegral)+"');";
    MsgInfo(MsgLog::Form("Insert Command = %s",insertCommand.c_str()));
    TSQLResult *res = db->Query(insertCommand.c_str());
    delete res;

    insertCommand = "insert into ccmdb.nearline_ped(insert_date,start_date,end_date";
    for (int i=0; i < 176; ++i) {
    insertCommand += ",pmt"+std::to_string(i);
    }
    insertCommand += ") value (NOW(),";
    insertCommand += "'"+std::string(bufferFirstTriggerTime)+"','"+std::string(bufferLastTriggerTime)+"'";
    for (int i=0; i < 176; ++i) {
    double value = speFinder->GetThreshold(i);
    if (value == 9999) {
    insertCommand += ",''";
    } else {
    insertCommand += ",'"+std::to_string(value)+"'";
    }
    }
    insertCommand += ");";
    res = db->Query(insertCommand.c_str());
    delete res;

    insertCommand = "insert into ccmdb.nearline_rates(insert_date,start_date,end_date";
    for (int i=0; i < 176; ++i) {
    insertCommand += ",pmt"+std::to_string(i);
    }
    insertCommand += ") value (NOW(),";
    insertCommand += "'"+std::string(bufferFirstTriggerTime)+"','"+std::string(bufferLastTriggerTime)+"'";
    for (int i=0; i < 176; ++i) {
    double value = speFinder->GetRate(i);
    if (value == 9999) {
    insertCommand += ",''";
    } else {
    insertCommand += ",'"+std::to_string(value)+"'";
    }
    }
    insertCommand += ");";
    res = db->Query(insertCommand.c_str());
    delete res;

    insertCommand = "insert into ccmdb.nearline_num_pmt_waveform(start_time,end_time";
    for (int i=0; i < 800; ++i) {
    insertCommand += ",bin"+std::to_string(i);
    }
    insertCommand += ") value ('"+std::string(bufferFirstTriggerTime)+"','"+std::string(bufferLastTriggerTime)+"'";
    for (const auto & value : totalNumPMTs) {
    insertCommand += ",'"+std::to_string(value)+"'";
    }
    insertCommand += ");";
    res = db->Query(insertCommand.c_str());
    delete res;

    insertCommand = "insert into ccmdb.nearline_integral_waveform(start_time,end_time";
    for (int i=0; i < 800; ++i) {
    insertCommand += ",bin"+std::to_string(i);
    }
    insertCommand += ") value ('"+std::string(bufferFirstTriggerTime)+"','"+std::string(bufferLastTriggerTime)+"'";
    for (const auto & value : totalIntegral) {
    insertCommand += ",'"+std::to_string(value)+"'";
    }
    insertCommand += ");";
    res = db->Query(insertCommand.c_str());

    delete res;
    */
    delete db;

    MsgInfo("Finished running over NearlineDiag");


    return EXIT_SUCCESS;

    }


