/*!**********************************************
 * \file CCMFindEvents.cxx
 * \author R.T. Thornton (LANL)
 * \date February 24,2020
 *
 * Main code to find events in the detector that
 * that are for a specific trigger
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMFindEvents.h"
#include "CCMModuleTable.h"

#include "RawData.h"
#include "Pulses.h"
#include "PMTInfoMap.h"
#include "SinglePulse.h"
#include "PMTInformation.h"
#include "SimplifiedEvent.h"
#include "MsgLog.h"
#include "Utility.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TChain.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <map>
#include <locale>

//See CCMModuleTable for info
MODULE_DECL(CCMFindEvents);

//_______________________________________________________________________________________
CCMFindEvents::CCMFindEvents(const char* version) 
  : CCMModule("CCMFindEvents"),
    fTriggerType("BEAM"),
    fHVOffList(""),
    fCalibrationFile(""),
    fOutFileName(""),
    fInFileName("")
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMFindEvents::CCMFindEvents(const CCMFindEvents& clufdr) 
: CCMModule(clufdr),
  fTriggerType(clufdr.fTriggerType),
  fHVOffList(clufdr.fHVOffList),
  fCalibrationFile(clufdr.fCalibrationFile),
  fOutFileName(clufdr.fOutFileName),
  fInFileName(clufdr.fInFileName)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMFindEvents::~CCMFindEvents()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMFindEvents::ProcessEvent()
{

  TH1::AddDirectory(0);

  // Create \p pulsesBeam as a smart pointer.
  // Used if \p triggerType == "beam"
  std::shared_ptr<Pulses> pulsesBeam = std::make_shared<Pulses>(0,0);

  // Fill #PMTInformation with default map
  PMTInfoMap::DefaultFillPMTMap();

  // Get the list of which HV are turned off.
  // Could move this to #PMTInfoMap::IsActive() function
  std::map<int,HighVoltage_t> highVoltage;
  PMTInfoMap::LoadHVOffList(fHVOffList);
  PMTInfoMap::LoadCalibrationFile(fCalibrationFile);

  // Setup TFile to read in input root file and find the
  // rawData and pulses TTrees
  TFile * file = TFile::Open(fInFileName.c_str(),"READ");
  TTree * tree = 0;
  TTree * treePulses = 0;

  file->GetObject("rawData",tree);
  file->GetObject("pulses",treePulses);

  if (!tree) {
    MsgError("rawData tree was not loaded");
  }
  if (!treePulses) {
    MsgError("pulses tree was not loaded");
  }

  // Set where to point the \p tree and \p treePulses branch pointers
  RawData * rawData = 0;
  Pulses * ccmPulses = 0;

  tree->SetBranchAddress("rawData",&rawData);
  treePulses->SetBranchAddress("pulses",&ccmPulses);

  long nEntries = tree->GetEntries();
  long skipped = 0;
  //long nEntries = 10;

  long numTriggers = 0;

  // Setup output files
  std::shared_ptr<TFile> outfile = std::make_shared<TFile>(Form("%s",fOutFileName.c_str()),"RECREATE");
  TTree * eventsTree = new TTree("events","events");
  SimplifiedEvent* events = new SimplifiedEvent();
  eventsTree->Branch("events",&events);

  // Set the number of bins and bin width
  // This is hard coded and should be taken 
  // from either data_structures.hh or a data base
  // since they (in principle) could change
  const int kNumBins = 8000;
  const double kBinWidth = 2e-3;

  // Define Histograms and vectors that will be filled with
  // various information below.
  TH1D * twoNSDist = new TH1D("twoNSDist","",10000,0,100);
  auto twoNSDistU = dynamic_cast<TH1D*>(twoNSDist->Clone("twoNSDistU"));
  auto twoNSDistC = dynamic_cast<TH1D*>(twoNSDist->Clone("twoNSDistC"));

  std::vector<float> pulsesTime(kNumBins,0);
  auto integralTime = pulsesTime;
  auto integralDer(integralTime);


  auto vetoBottomTime(pulsesTime);
  auto vetoTopTime(pulsesTime);
  auto vetoCTTime(pulsesTime);
  auto vetoCBTime(pulsesTime);

  auto itIntegralTimeBegin = integralTime.begin();
  auto itIntegralTimeEnd = integralTime.end();
  auto itIntegralDerBegin = integralDer.begin();
  auto itIntegralDerEnd = integralDer.end();
  auto itPulsesTimeBegin = pulsesTime.begin();
  auto itPulsesTimeEnd = pulsesTime.end();
  auto itVetoBottomTimeBegin = vetoBottomTime.begin();
  auto itVetoBottomTimeEnd = vetoBottomTime.end();
  auto itVetoCBTimeBegin = vetoCBTime.begin();
  auto itVetoCBTimeEnd = vetoCBTime.end();
  auto itVetoCTTimeBegin = vetoCTTime.begin();
  auto itVetoCTTimeEnd = vetoCTTime.end();
  auto itVetoTopTimeBegin = vetoTopTime.begin();
  auto itVetoTopTimeEnd = vetoTopTime.end();

  std::vector<std::vector<float>> pmtWaveform(160);
  for (int i=0; i < 160; ++i) {
    pmtWaveform[i].resize(pulsesTime.size(),0.0);
  }

  ///////////////////////////////////////////////
  // Time to loop through all the triggers in the file.
  // This loop builds the accumulated waveforms for the
  // various triggers that equate to the \p fTriggerType
  // that is passed from command line.
  ///////////////////////////////////////////////
  long e = 0;
  for (e = 0; e < nEntries; ++e) {
    if (e % 100 == 0) {
      MsgInfo(MsgLog::Form("At trigger %ld out of %ld",e,nEntries));
    }
    int nBytes = tree->GetEntry(e);
    if (nBytes <= 0) {
      continue;
    }
    if (!rawData) {
      MsgError("CCM RawData is set to 0");
      continue;
    }

    if (e % 1000 == 0) {
      outfile->Flush();
    }

    // Reset histograms and counters
    for (int i=0; i< 160; ++i) {
      std::fill(pmtWaveform[i].begin(),pmtWaveform[i].end(),0.0);
    }

    std::fill(itPulsesTimeBegin,itPulsesTimeEnd,0);
    std::fill(itIntegralTimeBegin,itIntegralTimeEnd,0.0);
    std::fill(itIntegralDerBegin,itIntegralDerEnd,0.0);
    std::fill(itVetoBottomTimeBegin,itVetoBottomTimeEnd,0);
    std::fill(itVetoTopTimeBegin,itVetoTopTimeEnd,0);
    std::fill(itVetoCTTimeBegin,itVetoCTTimeEnd,0);
    std::fill(itVetoCBTimeBegin,itVetoCBTimeEnd,0);

    // Check which trigger occured in DAQ window
    // make sure it is strobe
    size_t numWaveforms = rawData->GetNumWaveforms();
    int boardOffset = 0;
    if (numWaveforms != rawData->GetNumChannels()) {
      boardOffset = (rawData->GetNumBoards()-1)*rawData->GetNumChannels();
    }

    int firstSample = 0;
    if (fTriggerType.find("STROBE") != std::string::npos) {
      firstSample = Utility::FindFirstSample(boardOffset+1,rawData);
      if (!firstSample) {
        ++skipped;
        continue;
      }
    } else if (fTriggerType.find("BEAM") != std::string::npos) {
      firstSample = Utility::FindFirstSample(boardOffset+2,rawData);
      if (!firstSample) {
        ++skipped;
        continue;
      }
    } else if (fTriggerType.find("LED") != std::string::npos) {
      firstSample = Utility::FindFirstSample(boardOffset+3,rawData);
      if (!firstSample) {
        ++skipped;
        continue;
      }
    } else {
      MsgWarning(MsgLog::Form("Trigger type is %s and does not equal strobe, beam, or led. Skipping Event",fTriggerType.c_str()));
      continue;
    }

    // Load pulses Tree
    // find board to board adjustments for the 11th board
    nBytes = treePulses->GetEntry(e);
    if (nBytes <= 0) {
      continue;
    }
    if (!ccmPulses) {
      MsgError("CCM Pulses is set to 0");
      continue;
    }

    int beamTime = 8001;

    // If \p fTriggerType == beam
    // Find when the BCM triggered in the DAQ window 
    if (fTriggerType.find("BEAM") != std::string::npos) {
      double adjust11  = rawData->GetOffset(10) - ccmPulses->GetTriggerTime();
      std::shared_ptr<Pulses> pulsesBeam = std::make_shared<Pulses>(0,0);
      float offset = adjust11;
      pulsesBeam->DerivativeFilter(&rawData->GetSamples(0).front(),8000,0,3,2800,offset,1.0);

      for (size_t pulse = 0; pulse < pulsesBeam->GetNumPulses(); ++pulse) {
        double time = pulsesBeam->GetPulseTime(pulse);
        if (beamTime == 8001) {
          beamTime = time;
          break;
        }
      }
      //beamTimeHist->AddBinContent(beamTime,100.0);
      pulsesBeam->Reset();
    } else if (beamTime == 8001) {
      beamTime = 0;
    }

    // count trigger
    ++numTriggers;

    // loop through the pulses
    const size_t kNumPulses = ccmPulses->GetNumPulses();
    for (size_t loc = 0; loc < kNumPulses; ++loc) {

      // make sure the time is within the range of the DAQ window for analysis
      double time = ccmPulses->GetPulseTime(loc);
      if (std::isnan(time) || std::isinf(time) || time > 7960.0) {
        continue;
      }

      // make sure the PMT is in the database
      int key = ccmPulses->GetKey(loc);
      if (!PMTInfoMap::IsActive(key)) {
        continue;
      }

      const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(key);
      if (!pmtInfo) {
        continue;
      }

      // get the integral of the pulse and make sure it is a real number
      double pulseIntegral = ccmPulses->GetPulseIntegral(loc);
      if (std::isnan(pulseIntegral) || std::isinf(pulseIntegral)) {
        continue;
      }

      // get threshold and ADCtoPE for the PMT
      double threshold = pmtInfo->GetADCThreshold();
      double pe = pmtInfo->GetADCToPE();

      //time = time*kBinWidth - 9.92; // convert time from bin number to us

      // if the PMT was a veto pmt see if the integral is above 10
      // keep track the number of times that pmt fired in the DAQ window
      if (pmtInfo->IsVeto()) {
        auto name = pmtInfo->GetLocName();
        if (pulseIntegral  > 10) {
          int firstBin = std::max(time,0.0);
          int lastBin  = std::min(time + ccmPulses->GetPulseLength(loc),static_cast<double>(kNumBins));
          for (int bin = firstBin; bin <lastBin; ++bin) {
            if (bin < 0 || bin >= static_cast<int>(kNumBins)) {
              MsgFatal(MsgLog::Form("for veto pulse\tevent %ld pulse %zu key %d bin %d firstBin %d lastBin %d",e,loc,key,bin,firstBin,lastBin));
            }
            // veto top tubes
            if (name.find("VT") != std::string::npos) {
              ++vetoTopTime.at(bin);
              // veto column top tubes
            } else if (name.find("VCT") != std::string::npos) {
              ++vetoCTTime.at(bin);
              // veto column bottom tubes
            } else if (name.find("VCB") != std::string::npos) {
              ++vetoCBTime.at(bin);
            } else {
              ++vetoBottomTime.at(bin);
            } // end if-else over where the veto is located
          } // end for each bin the pulse is over
        } // end if pulseIntegral is greater than 10

        // move to the next event
        continue;
      } // end if pmtInfo->IsVeto()

      // check if pulse integral is below threshold or
      // the ADCtoPE calibration is negative or zero
      // if either is true do not analyze the pulse
      if (pulseIntegral < threshold || pe <= 0) {
        continue;
      }

      // check length of pulse
      // make sure it is greater than or equal to 20 ns
      double length = ccmPulses->GetPulseLength(loc);
      if (length*2.0 < 20.0) {
        continue;
      }

      // apply the ADCtoPE to the pulse
      pulseIntegral /= pe;

      // get the first and last bin in the DAQ window the pulse spanned
      int firstBin = std::max(time,0.0);
      int lastBin  = std::min(time + length,static_cast<double>(kNumBins));
      double middle = (lastBin + firstBin)/2.0;


      // add the pulse to the integral and pmt histograms
      // reconstruct the pulse for the integral histogram
      // as a triangle where the area of the triangle is
      // equal to to the pulseIntegral
      for (int bin = firstBin; bin < lastBin; ++bin) {
        double integral = 0.0;
        double hitFraction = 0.0;
        if (bin < middle) {
          integral = 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(bin - firstBin)/(middle - firstBin));
          hitFraction = 2.0/(lastBin-firstBin)*(static_cast<double>(bin - firstBin)/(middle - firstBin));
        } else if (bin != middle) {
          integral = 2.0*pulseIntegral/(lastBin-firstBin)*(static_cast<double>(lastBin - bin)/(lastBin - middle));
          hitFraction = 2.0/(lastBin-firstBin)*(static_cast<double>(lastBin - bin)/(lastBin - middle));
        } else {
          integral = 2.0*pulseIntegral/(lastBin-firstBin);
          hitFraction = 2.0/(lastBin-firstBin);
        }

        if (bin < 0 || bin >= static_cast<int>(kNumBins)) {
          MsgFatal(MsgLog::Form("event %ld pulse %zu key %d bin %d firstBin %d lastBin %d middle %f",e,loc,key,bin,firstBin,lastBin,middle));
        }
        pmtWaveform.at(key).at(bin) += integral;
        integralTime.at(bin) += integral;
        pulsesTime.at(bin) += hitFraction;

      } // end for bin = firstBin
    } // for over each pulse

    for (int sampleLoc=2; sampleLoc<static_cast<int>(kNumBins)-1; ++sampleLoc) {
      integralDer.at(sampleLoc) = (integralTime.at(sampleLoc+1) - integralTime.at(sampleLoc-1))/3.0/0.5;
    }

    /////////////////////////////////////////////
    // time to find events in the DAQ window
    /////////////////////////////////////////////
    int numBins20ns = 0.02/kBinWidth;
    int numBins90ns = 0.090/kBinWidth;
    //int numBins1p6us = 1.6/kBinWidth;
    //int prevStart = -1000;

    for (size_t bin = 0; bin < kNumBins; ++bin) {
      if (-9.92 + static_cast<double>(bin)*kBinWidth >= static_cast<double>(7960-1600)*2e-3 - 9.92) {
        break;
      }
      if (integralTime.at(bin) >= 0.4) {
        TVector3 pos;
        int startBin = bin;

        if (startBin >= static_cast<int>(kNumBins)) {
          MsgFatal(MsgLog::Form("event %ld bin %zu startBin %d",e,bin,startBin));
        }

        // find end of event: defined as two consecutive empty bins
        int endBin = kNumBins-1;
        for (size_t bin2 = bin; bin2 < kNumBins-10; ++bin2) {
          if (std::accumulate(itPulsesTimeBegin+bin2,itPulsesTimeBegin+bin2+10,0.f) == 0) {
            endBin = bin2;
            break;
          } // end if no hits for 20 ns
        } // end for size_t bin2 = bin

        if (endBin >= static_cast<int>(kNumBins)) {
          MsgFatal(MsgLog::Form("event %ld bin %zu endBin %d",e,bin,endBin));
        }

        int start = startBin - numBins90ns;
        if (start < 0) {
          start = 0;
        }

        int vetoActivityTop = std::accumulate(itVetoTopTimeBegin+start,itVetoTopTimeBegin+endBin,0);
        int vetoActivityCT = std::accumulate(itVetoCTTimeBegin+start,itVetoCTTimeBegin+endBin,0);
        int vetoActivityCB = std::accumulate(itVetoCBTimeBegin+start,itVetoCBTimeBegin+endBin,0);
        int vetoActivityBottom = std::accumulate(itVetoBottomTimeBegin+start,itVetoBottomTimeBegin+endBin,0);

        //int maxBin = std::distance(itIntegralTimeBegin,std::max_element(itIntegralTimeBegin+startBin,itIntegralTimeBegin+endBin));

        //prevStart = startBin;
        // + 15 to calibrate to EJ detectors
        // subtract the beam time to remove the jitter of the beam timing
        double startTime = static_cast<double>(startBin-beamTime+15-1)*kBinWidth;
        double endTime = static_cast<double>(endBin-beamTime+15-1)*kBinWidth;

        double promptIntegral = 0;
        if (startBin+numBins90ns >= kNumBins) {
          promptIntegral = std::accumulate(itIntegralTimeBegin+startBin,itIntegralTimeEnd,0.0);
        } else {
          promptIntegral = std::accumulate(itIntegralTimeBegin+startBin,itIntegralTimeBegin+startBin+numBins90ns,0.0);
        }

        //double integral = std::accumulate(itIntegralTimeBegin+startBin,itIntegralTimeBegin+endBin,0.0);

        double promptCoated = 0.0;
        double promptUncoated = 0.0;
        double totalCoated = 0.0;
        double totalUncoated = 0.0;
        double promptFit = 0;
        std::vector<float> percentOfPMT;
        double largestPMTFraction = 0.0;
        //double largestPMTCharge = 0.0;
        //int largestPMT = 0.0;
        int numPMTs = 0;

        for (int i=0; i < 160; ++i) {
          const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(i);
          if (!pmtInfo) {
            continue;
          }
          if (pmtInfo->IsVeto()) {
            continue;
          }
          double prompt = 0;
          double prompt90 = 0;
          double total = 0; 
          // find integral for the first 20 ns
          if (startBin+numBins20ns < kNumBins) {
            prompt = std::accumulate(pmtWaveform.at(i).begin()+startBin,pmtWaveform.at(i).begin()+startBin+numBins20ns,0.0);
          } else {
            prompt = std::accumulate(pmtWaveform.at(i).begin()+startBin,pmtWaveform.at(i).end(),0.0);
          }

          // find integral for the first 90 ns
          if (startBin+numBins90ns < kNumBins) {
            prompt90 = std::accumulate(pmtWaveform.at(i).begin()+startBin,pmtWaveform.at(i).begin()+startBin+numBins90ns,0.0);
          } else {
            prompt90= std::accumulate(pmtWaveform.at(i).begin()+startBin,pmtWaveform.at(i).end(),0.0);
          }

          // find integral for the entire event window
          if (endBin < kNumBins) {
            total = std::accumulate(pmtWaveform.at(i).begin()+startBin,pmtWaveform.at(i).begin()+endBin,0.0);
          }  else {
            total = std::accumulate(pmtWaveform.at(i).begin()+startBin,pmtWaveform.at(i).end(),0.0);
          }

          if (prompt90 != 0) {
            ++numPMTs;
            percentOfPMT.push_back(prompt90/promptIntegral);
          }
          if (prompt90/promptIntegral > largestPMTFraction) {
            largestPMTFraction = prompt90/promptIntegral;
            //largestPMT = i;
            //largestPMTCharge = prompt90;
          }

          if (pmtInfo->IsUncoated()) {
            promptUncoated += prompt90;
            totalUncoated += total;
          } else {
            promptCoated += prompt90;
            totalCoated += total;
          }
          if (prompt) {
            pos += *(pmtInfo->GetPosition())*prompt*prompt;
            promptFit += prompt*prompt;
          }
        } // end for over all digitizer channels
        pos *= 1.0/promptFit;
        events->SetTriggerNumber(rawData->GetEventNumber());
        events->SetComputerSecIntoEpoch(rawData->GetGPSSecIntoDay());
        events->SetComputerNSIntoSec(rawData->GetGPSNSIntoSec());
        events->SetStartTime(startTime);
        if (startBin+numBins90ns >= kNumBins) {
          events->SetNumCoated(std::accumulate(itPulsesTimeBegin+startBin,itPulsesTimeEnd,0.f),true);
        } else {
          events->SetNumCoated(std::accumulate(itPulsesTimeBegin+startBin,itPulsesTimeBegin+startBin+numBins90ns,0.f),true);
        }
        if (endBin < kNumBins) {
          events->SetNumCoated(std::accumulate(itPulsesTimeBegin+startBin,itPulsesTimeBegin+endBin,0.f),false);
          std::vector<float> hits(itPulsesTimeBegin+startBin,itPulsesTimeBegin+endBin);
          std::vector<float> energy(itIntegralTimeBegin+startBin,itIntegralTimeBegin+endBin);
          events->AddWaveforms(hits,energy);
        } else {
          events->SetNumCoated(std::accumulate(itPulsesTimeBegin+startBin,itPulsesTimeEnd,0.f),false);
          std::vector<float> hits(itPulsesTimeBegin+startBin,itPulsesTimeEnd);
          std::vector<float> energy(itIntegralTimeBegin+startBin,itIntegralTimeEnd);
          events->AddWaveforms(hits,energy);
        }
        events->SetLength(endTime-startTime);
        events->SetLargestPMTFraction(largestPMTFraction);
        events->SetIntegralCoated(promptCoated,true);
        events->SetIntegralCoated(totalCoated,false);
        events->SetIntegralUncoated(promptUncoated,true);
        events->SetIntegralUncoated(totalUncoated,false);
        events->SetXPosition(pos.X());
        events->SetYPosition(pos.Y());
        events->SetZPosition(pos.Z());
        events->SetNumVetoTop(vetoActivityTop);
        events->SetNumVetoBottom(vetoActivityBottom);
        events->SetNumVetoBack(vetoActivityCB);
        events->SetNumVetoFront(vetoActivityCT);
        outfile->cd();
        eventsTree->Fill();
        events->Reset();
        double avgPercent = std::accumulate(percentOfPMT.begin(),percentOfPMT.end(),0.0);
        avgPercent /= static_cast<double>(percentOfPMT.size());
        double sigmaPercent = 0.0;
        double skewPercent = 0.0;
        for (const auto & pmt : percentOfPMT) {
          sigmaPercent += std::pow(pmt-avgPercent,2.0);
          skewPercent += std::pow(pmt-avgPercent,3.0);
        }
        skewPercent = skewPercent/static_cast<double>(percentOfPMT.size())/std::pow(1.0/(static_cast<double>(percentOfPMT.size())-1.0)*sigmaPercent,3.0/2.0);
        sigmaPercent = std::sqrt(1.0/(static_cast<double>(percentOfPMT.size())-1.0)*sigmaPercent);

        bin = endBin;
      } else {
        twoNSDist->Fill(integralTime.at(bin));
        if (bin >= 0 && bin < pmtWaveform.front().size()) {
          for (int i=0; i < 160; ++i) {
            const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(i);
            if (!pmtInfo) {
              continue;
            }
            if (pmtInfo->IsVeto()) {
              continue;
            }
            double pe = pmtInfo->GetADCToPE();
            if (pe <= 0) {
              continue;
            }
            if (pmtInfo->IsUncoated()) {
              twoNSDistU->Fill(pmtWaveform.at(i).at(bin));
            } else {
              twoNSDistC->Fill(pmtWaveform.at(i).at(bin));
            } //end if-else IsUncoated
          } // end for over all digitizer channels
        } // end if bin is in vector range
      } // end threshold if-else
    } // end for int bin
  } // end for e < nEntries

  // Save tree and histograms to file
  outfile->cd();
  eventsTree->Write();
  twoNSDist->Write();
  twoNSDistC->Write();
  twoNSDistU->Write();
  outfile->Close();

  MsgInfo(MsgLog::Form("Num Triggers %ld",numTriggers));

  if (file) {
    delete file;
  }

  /////////////////////////////////////////////
  // delete histograms
  /////////////////////////////////////////////
  if (twoNSDist) {
    delete twoNSDist;
  }
  if (twoNSDistC) {
    delete twoNSDistC;
  }
  if (twoNSDistU) {
    delete twoNSDistU;
  }

  PMTInfoMap::ClearMap();

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMFindEvents::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  if (&c != 0)
  {
    c("TriggerType").Get(fTriggerType);

    std::locale loc;
    for (auto & c : fTriggerType) {
      std::toupper(c,loc);
    }

    c("HVOffList").Get(fHVOffList);
    c("CalibrationFile").Get(fCalibrationFile);
    c("OutFileName").Get(fOutFileName);
    c("InFileName").Get(fInFileName);

    fIsInit = true;
  } else {
    fIsInit = false;
  }

}

