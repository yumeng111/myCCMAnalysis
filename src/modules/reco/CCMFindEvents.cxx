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
#include "PMTInformation.h"
#include "SimplifiedEvent.h"
#include "Events.h"
#include "MsgLog.h"
#include "Utility.h"

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
    fThreshold(0.4),
    fEvents(nullptr),
    fRawData(nullptr),
    fPulses(nullptr),
    fNumTriggers(0)
{
  //Default constructor
  this->SetCfgVersion(version);
  DefineVectors();
}

//_______________________________________________________________________________________
CCMFindEvents::CCMFindEvents(const CCMFindEvents& clufdr) 
: CCMModule(clufdr),
  fTriggerType(clufdr.fTriggerType),
  fThreshold(clufdr.fThreshold),
  fEvents(clufdr.fEvents),
  fRawData(clufdr.fRawData),
  fPulses(clufdr.fPulses),
  fNumTriggers(clufdr.fNumTriggers)
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

  // create new Events object. This will delete
  // if a previous object was created
  fEvents->Reset();

  fEvents->SetEventNumber(fRawData->GetEventNumber());
  fEvents->SetComputerSecIntoEpoch(fRawData->GetGPSSecIntoDay());
  fEvents->SetComputerNSIntoSec(fRawData->GetGPSNSIntoSec());

  std::unique_ptr<SimplifiedEvent> event = std::make_unique<SimplifiedEvent>();

  auto itIntegralTimeBegin = fIntegralTime.begin();
  auto itIntegralTimeEnd = fIntegralTime.end();
  auto itIntegralDerBegin = fIntegralDer.begin();
  auto itIntegralDerEnd = fIntegralDer.end();
  auto itPulsesTimeBegin = fPulsesTime.begin();
  auto itPulsesTimeEnd = fPulsesTime.end();
  auto itVetoBottomTimeBegin = fVetoBottomTime.begin();
  auto itVetoBottomTimeEnd = fVetoBottomTime.end();
  auto itVetoCLeftTimeBegin = fVetoCLeftTime.begin();
  auto itVetoCLeftTimeEnd = fVetoCLeftTime.end();
  auto itVetoCRightTimeBegin = fVetoCRightTime.begin();
  auto itVetoCRightTimeEnd = fVetoCRightTime.end();
  auto itVetoCFrontTimeBegin = fVetoCFrontTime.begin();
  auto itVetoCFrontTimeEnd = fVetoCFrontTime.end();
  auto itVetoCBackTimeBegin = fVetoCBackTime.begin();
  auto itVetoCBackTimeEnd = fVetoCBackTime.end();
  auto itVetoTopTimeBegin = fVetoTopTime.begin();
  auto itVetoTopTimeEnd = fVetoTopTime.end();
  auto itVetoTotalTimeBegin = fVetoTotalTime.begin();
  auto itVetoTotalTimeEnd = fVetoTotalTime.end();

  // Reset histograms and counters
  for (auto & waveform : fPMTWaveform) {
    std::fill(waveform.begin(),waveform.end(),0.f);
  }
  for (auto & waveform : fPMTWaveformCount) {
    std::fill(waveform.begin(),waveform.end(),0);
  }

  std::fill(itPulsesTimeBegin,itPulsesTimeEnd,0);
  std::fill(itIntegralTimeBegin,itIntegralTimeEnd,0.f);
  std::fill(itIntegralDerBegin,itIntegralDerEnd,0.f);

  std::fill(itVetoBottomTimeBegin,itVetoBottomTimeEnd,0);
  std::fill(itVetoTopTimeBegin,itVetoTopTimeEnd,0);
  std::fill(itVetoTotalTimeBegin,itVetoTotalTimeEnd,0);
  std::fill(itVetoCRightTimeBegin,itVetoCRightTimeEnd,0);
  std::fill(itVetoCFrontTimeBegin,itVetoCFrontTimeEnd,0);
  std::fill(itVetoCBackTimeBegin,itVetoCBackTimeEnd,0);
  std::fill(itVetoCLeftTimeBegin,itVetoCLeftTimeEnd,0);

  // Check which trigger occured in DAQ window
  // make sure it is strobe
  if (!fRawData->IsTriggerPresent(fTriggerType)) {
    return kCCMDoNotWrite;
  }


  // If fTriggerType == beam
  // Find when the BCM triggered in the DAQ window 
  int beamTime = fgkNumBins+1;
  if (fTriggerType.find("BEAM") != std::string::npos) {
    beamTime = fRawData->GetBCMTime();
  } else if (beamTime == fgkNumBins+1) {
    beamTime = 0;
  }

  // count trigger
  ++fNumTriggers;

  // declear variables used in the for loops
  TVector3 pos;
  auto pmtInfo = PMTInfoMap::GetPMTInfo(0);
  auto pmtWaveformBeginIter = fPMTWaveform.front().begin();

  double c0 = 0;
  double c1 = 0;
  double hitFraction = 0;
  double integral = 0;
  double largestPMTFraction = 0.0;
  double lastBin = 0;
  double length = 0;
  double middle = 0;
  double pe = 0;
  double prompt = 0;
  double prompt90 = 0;
  double promptCoated = 0.0;
  double promptFit = 0;
  double promptIntegral = 0;
  double promptUncoated = 0.0;
  double pulseIntegral = 0;
  double threshold = 0;
  double time = 0;
  double timeCoated = 0.0;
  double timeUncoated = 0.0;
  double total = 0;
  double totalCoated = 0.0;
  double totalUncoated = 0.0;
  double vetoTimeBottom = 0;
  double vetoTimeCBack = 0;
  double vetoTimeCFront = 0;
  double vetoTimeCLeft = 0;
  double vetoTimeCRight = 0;
  double vetoTimeTop = 0;

  int bin = 0;
  int bin2 = 0;
  int current = 0;
  int end20nsBin = 0;
  int end90nsBin = 0;
  int endBin = 0;
  int endTotalBin = 0;
  int firstBin = 0;
  int key = 0;
  int maxVeto = 0;
  int maxVetoEnd = 0;
  int maxVetoPrompt = 0;
  int maxVetoPromptEnd = 0;
  int maxVetoPromptStart = 0;
  int maxVetoStart = 0;
  int numPMTs = 0;
  int peakLoc = 0;
  int pmt = 0;
  int start = 0;
  int startBin = 0;
  int vetoActivityBottom = 0;
  int vetoActivityCBack = 0;
  int vetoActivityCFront = 0;
  int vetoActivityCLeft = 0;
  int vetoActivityCRight = 0;
  int vetoActivityPromptBottom = 0;
  int vetoActivityPromptCBack = 0;
  int vetoActivityPromptCFront = 0;
  int vetoActivityPromptCLeft = 0;
  int vetoActivityPromptCRight = 0;
  int vetoActivityPromptTop = 0;
  int vetoActivityTop = 0;
  int vetoBin = 0;
  int vetoIntLength = 0;

  std::vector<double> xPoints;
  std::vector<double> yPoints;
  std::vector<float> energy;
  std::vector<float> hits;
  std::vector<float> percentOfPMT;
  std::vector<int> count;



  // loop through the pulses
  const size_t kNumPulses = fPulses->GetNumPulses();
  for (size_t loc = 0; loc < kNumPulses; ++loc) {

    // make sure the time is within the range of the DAQ window for analysis
    time = fPulses->GetPulseTime(loc);
    if (std::isnan(time) || std::isinf(time) || time > 7960.0) {
      continue;
    }

    // make sure the PMT is in the database
    key = fPulses->GetKey(loc);
    if (!PMTInfoMap::IsActive(key)) {
      continue;
    }

    pmtInfo = PMTInfoMap::GetPMTInfo(key);
    if (!pmtInfo) {
      continue;
    }

    // get the integral of the pulse and make sure it is a real number
    pulseIntegral = fPulses->GetPulseIntegral(loc);
    if (std::isnan(pulseIntegral) || std::isinf(pulseIntegral)) {
      continue;
    }

    // get threshold and ADCtoPE for the PMT
    threshold = pmtInfo->GetADCThreshold();
    pe = pmtInfo->GetADCToPE();

    //time = time*fgkBinWidth - 9.92; // convert time from bin number to us

    // if the PMT was a veto pmt see if the integral is above 10
    // keep track the number of times that pmt fired in the DAQ window
    if (pmtInfo->IsVeto()) {
      auto name = pmtInfo->GetLocName();
      if (pulseIntegral  > 5) {
        firstBin = std::max(time,0.0);
        ++fPMTWaveformCount.at(key).at(firstBin);
        ++fPMTWaveform.at(key).at(firstBin);
        // veto top tubes
        if (name.find("VT") != std::string::npos) {
          ++fVetoTopTime.at(firstBin);
          // veto column top tubes
        } else if (name.find("VC") != std::string::npos) {
          // Back (columes 22 through 3)
          // Left (columes 4 through 9)
          // Right (columes 16 through 21)
          // Front (columes 10 through 15)
          int column = pmtInfo->GetColumn();
          if (column >= 4 && column <= 9) {
            ++fVetoCLeftTime.at(firstBin);
          } else if (column >= 10 && column <= 15) {
            ++fVetoCFrontTime.at(firstBin);
          } else if (column >= 16 && column <= 21) {
            ++fVetoCRightTime.at(firstBin);
          } else if (column > 0) {
            ++fVetoCBackTime.at(firstBin);
          }
        } else {
          ++fVetoBottomTime.at(firstBin);
        } // end if-else over where the veto is located
        ++fVetoTotalTime.at(firstBin);
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
    length = fPulses->GetPulseLength(loc);
    if (length*2.0 < 20.0) {
      continue;
    }

    // apply the ADCtoPE to the pulse
    pulseIntegral /= pe;

    // get the first and last bin in the DAQ window the pulse spanned
    firstBin = std::max(time,0.0);
    lastBin  = std::min(time + length,static_cast<double>(fgkNumBins));
    middle = (lastBin + firstBin)/2.0;


    // add the pulse to the integral and pmt histograms
    // reconstruct the pulse for the integral histogram
    // as a triangle where the area of the triangle is
    // equal to to the pulseIntegral
    ++fPMTWaveformCount.at(key).at(firstBin);
    for (bin = firstBin; bin < lastBin; ++bin) {
      integral = 0.0;
      hitFraction = 0.0;
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

      if (bin < 0 || bin >= static_cast<int>(fgkNumBins)) {
        MsgFatal(MsgLog::Form("event %ld pulse %zu key %d bin %d firstBin %d lastBin %d middle %f",
              fPulses->GetEventNumber(),loc,key,bin,firstBin,lastBin,middle));
      }
      fPMTWaveform.at(key).at(bin) += integral;
      fIntegralTime.at(bin) += integral;
      fPulsesTime.at(bin) += hitFraction;

    } // end for bin = firstBin
  } // for over each pulse

  for (int sampleLoc=2; sampleLoc<static_cast<int>(fgkNumBins)-1; ++sampleLoc) {
    fIntegralDer.at(sampleLoc) = (fIntegralTime.at(sampleLoc+1) - fIntegralTime.at(sampleLoc-1))/3.0/0.5;
  }

  /////////////////////////////////////////////
  // time to find events in the DAQ window
  /////////////////////////////////////////////
  int numBins20ns = 0.02/fgkBinWidth;
  int numBins90ns = 0.090/fgkBinWidth;
  //int numBins30ns = 0.030/fgkBinWidth;
  //int numBins1p6us = 1.6/fgkBinWidth;
  //int prevStart = -1000;

  for (size_t timeBin = 0; timeBin < fgkNumBins; ++timeBin) {
    MsgDebug(3,MsgLog::Form("New Time bin = %zu time = %.3f limit = %.3f",timeBin,
          -9.92 + static_cast<double>(timeBin)*fgkBinWidth,
          static_cast<double>(7960-1600)*2e-3 - 9.92));
    if (-9.92 + static_cast<double>(timeBin)*fgkBinWidth >= static_cast<double>(7960-1600)*2e-3 - 9.92) {
      break;
    }
    MsgDebug(3,"Passed DAQ time cut");
    if (fIntegralTime.at(timeBin) >= fThreshold) {
      MsgDebug(2,"Above threshold");

      startBin = timeBin;

      if (startBin >= static_cast<int>(fgkNumBins)) {
        MsgFatal(MsgLog::Form("Start Bin is greater than or equal to total number of bins: event %ld bin %zu startBin %d",
              fPulses->GetEventNumber(),timeBin,startBin));
      } else {
        MsgDebug(2,MsgLog::Form("event %ld bin %zu start bin = %d",fPulses->GetEventNumber(),timeBin,startBin));
      }

      // find where the first peak is in the event
      peakLoc = startBin;
      for (bin2 = startBin; bin2 < fgkNumBins-1; ++bin2) {
        if (fIntegralTime.at(bin2+1) < fIntegralTime.at(bin2)) {
          peakLoc = bin2;
          break;
        } // end if fIntegralTime.at(bin2+1) < ...
      } // end for size_t bin2 = startBin
      MsgDebug(2,MsgLog::Form("Peak Location = %d",peakLoc));

      if ((peakLoc+startBin)/2-startBin > 3) {
        for (bin2 = startBin; bin2 < (peakLoc+startBin)/2; ++bin2) {
          MsgDebug(4,MsgLog::Form("Peak location loop bin %d",bin2));
          xPoints.push_back(bin2);
          yPoints.push_back(fIntegralTime.at(bin2));
        }
        Utility::LinearUnweightedLS(xPoints.size(),&xPoints.front(),&yPoints.front(),c0,c1);
        MsgDebug(4,"Did LinearUnweightedLS");

        // find the first non-empty start bin
        // important for calculating position
        // and integrals
        MsgDebug(4,MsgLog::Form("Fit results: -c0 = %.3f c1 = %.3f start = %d",-c0,c1,static_cast<int>(-c0/c1)));
        if (-c0/c1 > 0) {
          startBin = static_cast<int>(-c0/c1);
          while (fPulsesTime.at(startBin) == 0) {
            ++startBin;
          }
          MsgDebug(4,MsgLog::Form("Final startBin = %d",startBin));
        }

        xPoints.clear();
        yPoints.clear();
      }

      // find end of event: defined as two consecutive empty bins
      endBin = fgkNumBins-1;
      for (bin2 = timeBin; bin2 < fgkNumBins-10; ++bin2) {
        if (std::accumulate(itPulsesTimeBegin+bin2,itPulsesTimeBegin+bin2+10,0.f) == 0) {
          endBin = bin2;
          break;
        } // end if no hits for 20 ns
      } // end for size_t bin2 = bin

      if (endBin >= static_cast<int>(fgkNumBins)) {
        MsgFatal(MsgLog::Form("event %ld bin %zu endBin %d",fPulses->GetEventNumber(),timeBin,endBin));
      } else {
        MsgDebug(2,MsgLog::Form("event %ld bin %zu start bin = %d and end bin = %d",fPulses->GetEventNumber(),timeBin,startBin,endBin));
      }

      start = startBin - numBins90ns;
      if (start < 0) {
        start = 0;
      }
      end90nsBin = (startBin+numBins90ns < fgkNumBins) ? startBin+numBins90ns : fgkNumBins-1;
      end20nsBin = (startBin+numBins20ns < fgkNumBins) ? startBin+numBins20ns : fgkNumBins-1;
      endTotalBin = (endBin < fgkNumBins) ? endBin : fgkNumBins-1;

      MsgDebug(2,MsgLog::Form("end90nsBin = %d end20nsBin = %d endTotalBin = %d",end90nsBin,end20nsBin,endTotalBin));


      maxVeto = 0;
      maxVetoStart = 0;
      maxVetoEnd = 0;
      maxVetoPrompt = 0;
      maxVetoPromptStart = 0;
      maxVetoPromptEnd = 0;
      vetoIntLength = std::min(numBins90ns,endBin-start);
      MsgDebug(2,MsgLog::Form("Length of veto integral = %d",vetoIntLength));
      for (vetoBin = start; vetoBin+vetoIntLength <= endBin; ++vetoBin) {
        current = std::accumulate(itVetoTotalTimeBegin+vetoBin,itVetoTotalTimeBegin+vetoBin+vetoIntLength,0);
        if (current > maxVeto) {
          maxVeto = current;
          maxVetoStart = vetoBin;
          maxVetoEnd = vetoBin+vetoIntLength;
        }
        if (vetoBin <= startBin) {
          current = std::accumulate(itVetoTotalTimeBegin+vetoBin,itVetoTotalTimeBegin+vetoBin+vetoIntLength,0);
          if (current > maxVetoPrompt) {
            maxVetoPrompt = current;
            maxVetoPromptStart = vetoBin;
            maxVetoPromptEnd = vetoBin+vetoIntLength;
          }
        }
      }

      MsgDebug(1,MsgLog::Form("maxVetoStart %d maxVetoEnd %d",maxVetoStart,maxVetoEnd));
      vetoActivityTop = std::accumulate(itVetoTopTimeBegin+maxVetoStart,itVetoTopTimeBegin+maxVetoEnd,0);
      vetoActivityCRight = std::accumulate(itVetoCRightTimeBegin+maxVetoStart,itVetoCRightTimeBegin+maxVetoEnd,0);
      vetoActivityCLeft = std::accumulate(itVetoCLeftTimeBegin+maxVetoStart,itVetoCLeftTimeBegin+maxVetoEnd,0);
      vetoActivityCFront = std::accumulate(itVetoCFrontTimeBegin+maxVetoStart,itVetoCFrontTimeBegin+maxVetoEnd,0);
      vetoActivityCBack = std::accumulate(itVetoCBackTimeBegin+maxVetoStart,itVetoCBackTimeBegin+maxVetoEnd,0);
      vetoActivityBottom = std::accumulate(itVetoBottomTimeBegin+maxVetoStart,itVetoBottomTimeBegin+maxVetoEnd,0);

      MsgDebug(1,"Veto prompt integral");
      vetoActivityPromptTop = std::accumulate(itVetoTopTimeBegin+maxVetoPromptStart,itVetoTopTimeBegin+maxVetoPromptEnd,0);
      vetoActivityPromptCRight = std::accumulate(itVetoCRightTimeBegin+maxVetoPromptStart,itVetoCRightTimeBegin+maxVetoPromptEnd,0);
      vetoActivityPromptCLeft = std::accumulate(itVetoCLeftTimeBegin+maxVetoPromptStart,itVetoCLeftTimeBegin+maxVetoPromptEnd,0);
      vetoActivityPromptCFront = std::accumulate(itVetoCFrontTimeBegin+maxVetoPromptStart,itVetoCFrontTimeBegin+maxVetoPromptEnd,0);
      vetoActivityPromptCBack = std::accumulate(itVetoCBackTimeBegin+maxVetoPromptStart,itVetoCBackTimeBegin+maxVetoPromptEnd,0);
      vetoActivityPromptBottom = std::accumulate(itVetoBottomTimeBegin+maxVetoPromptStart,itVetoBottomTimeBegin+maxVetoPromptEnd,0);

      MsgDebug(1,"Get first none empty bin for the veto tubes");
      vetoTimeTop = Utility::FindFirstNoneEmptyBin<int>(itVetoTopTimeBegin,itVetoTopTimeBegin+start,itVetoTopTimeBegin+endBin);
      vetoTimeCRight = Utility::FindFirstNoneEmptyBin<int>(itVetoCRightTimeBegin,itVetoCRightTimeBegin+start,itVetoCRightTimeBegin+endBin);
      vetoTimeCLeft = Utility::FindFirstNoneEmptyBin<int>(itVetoCLeftTimeBegin,itVetoCLeftTimeBegin+start,itVetoCLeftTimeBegin+endBin);
      vetoTimeCFront = Utility::FindFirstNoneEmptyBin<int>(itVetoCFrontTimeBegin,itVetoCFrontTimeBegin+start,itVetoCFrontTimeBegin+endBin);
      vetoTimeCBack = Utility::FindFirstNoneEmptyBin<int>(itVetoCBackTimeBegin,itVetoCBackTimeBegin+start,itVetoCBackTimeBegin+endBin);
      vetoTimeBottom = Utility::FindFirstNoneEmptyBin<int>(itVetoBottomTimeBegin,itVetoBottomTimeBegin+start,itVetoBottomTimeBegin+endBin);

      // + 15 to calibrate to EJ detectors
      // subtract the beam time to remove the jitter of the beam timing
      if (vetoTimeTop >= 0) {
        vetoTimeTop = static_cast<double>(vetoTimeTop - beamTime + 15)*fgkBinWidth;
      }
      if (vetoTimeCRight >= 0) {
        vetoTimeCRight = static_cast<double>(vetoTimeCRight - beamTime + 15)*fgkBinWidth;
      }
      if (vetoTimeCLeft >= 0) {
        vetoTimeCLeft = static_cast<double>(vetoTimeCLeft - beamTime + 15)*fgkBinWidth;
      }
      if (vetoTimeCFront >= 0) {
        vetoTimeCFront = static_cast<double>(vetoTimeCRight - beamTime + 15)*fgkBinWidth;
      }
      if (vetoTimeCBack >= 0) {
        vetoTimeCBack = static_cast<double>(vetoTimeCBack - beamTime + 15)*fgkBinWidth;
      }
      if (vetoTimeBottom >= 0) {
        vetoTimeBottom = static_cast<double>(vetoTimeBottom - beamTime + 15)*fgkBinWidth;
      }
      MsgDebug(2,"Got first none empty bin for the veto tubes");

      double startTime = static_cast<double>(startBin-beamTime+15)*fgkBinWidth;
      double endTime = static_cast<double>(endBin-beamTime+15)*fgkBinWidth;

      //MsgInfo(MsgLog::Form("Max Veto Location(Prompt) %d(%d) Amount(Prompt) %d(%d)",
      //      maxVetoStart,maxVetoPromptStart,maxVeto,maxVetoPrompt));

      promptIntegral = std::accumulate(itIntegralTimeBegin+startBin,itIntegralTimeBegin+end90nsBin,0.f);

      timeCoated = 0.0;
      timeUncoated = 0.0;
      promptCoated = 0.0;
      promptUncoated = 0.0;
      totalCoated = 0.0;
      totalUncoated = 0.0;
      promptFit = 0;
      if (!percentOfPMT.empty()) {
        percentOfPMT.clear();
      }
      largestPMTFraction = 0.0;
      numPMTs = 0;

      MsgDebug(3,"Loop over pmts");
      for (pmt=0; pmt < fgkNumPMTs; ++pmt) {
        pmtInfo = PMTInfoMap::GetPMTInfo(pmt);
        if (!pmtInfo) {
          continue;
        }
        if (pmtInfo->IsVeto()) {
          continue;
        }

        pmtWaveformBeginIter = fPMTWaveform.at(pmt).begin();

        // find integral for the first 20 ns
        MsgDebug(4,Form("PMT %d first 20ns",pmt));
        prompt = std::accumulate(pmtWaveformBeginIter+startBin,pmtWaveformBeginIter+end20nsBin,0.f);

        // find integral for the first 90 ns
        MsgDebug(4,Form("PMT %d first 90ns",pmt));
        prompt90 = std::accumulate(pmtWaveformBeginIter+startBin,pmtWaveformBeginIter+end90nsBin,0.f);

        // find integral for the entire event window
        MsgDebug(4,Form("PMT %d total",pmt));
        total = std::accumulate(pmtWaveformBeginIter+startBin,pmtWaveformBeginIter+endTotalBin,0.f);
        time = Utility::FindFirstNoneEmptyBin<float>(pmtWaveformBeginIter,pmtWaveformBeginIter+startBin,pmtWaveformBeginIter+endTotalBin);

        if (prompt90 != 0) {
          ++numPMTs;
          percentOfPMT.push_back(prompt90/promptIntegral);
        }
        if (prompt90 > largestPMTFraction) {
          largestPMTFraction = prompt90;
          //largestPMT = pmt;
          //largestPMTCharge = prompt90;
        }

        if (pmtInfo->IsUncoated()) {
          promptUncoated += prompt90;
          totalUncoated += total;
          if (time != -1) {
            timeUncoated = std::min(time,timeUncoated);
          }
        } else {
          promptCoated += prompt90;
          totalCoated += total;
          if (time != -1) {
            timeCoated = std::min(time,timeCoated);
          }
        }
        if (prompt) {
          pos += *(pmtInfo->GetPosition())*prompt*prompt;
          promptFit += prompt*prompt;
        }
      } // end for over all digitizer channels
      pos *= 1.0/promptFit;

      timeUncoated = static_cast<double>(timeUncoated - beamTime + 15)*fgkBinWidth;
      timeCoated = static_cast<double>(timeCoated - beamTime + 15)*fgkBinWidth;

      MsgDebug(2,"Set events parameters");
      event->SetStartTime(startTime);
      event->SetStartTimeUncoated(timeUncoated);
      event->SetStartTimeCoated(timeCoated);
      event->SetStartTimeVetoTop(vetoTimeTop);
      event->SetStartTimeVetoRight(vetoTimeCRight);
      event->SetStartTimeVetoLeft(vetoTimeCLeft);
      event->SetStartTimeVetoFront(vetoTimeCFront);
      event->SetStartTimeVetoBack(vetoTimeCBack);
      event->SetStartTimeVetoBottom(vetoTimeBottom);

      event->SetNumCoated(std::accumulate(itPulsesTimeBegin+startBin,itPulsesTimeBegin+end90nsBin,0.f),true);
      event->SetNumCoated(std::accumulate(itPulsesTimeBegin+startBin,itPulsesTimeBegin+endTotalBin,0.f),false);

      hits.resize(endTotalBin-startBin);
      energy.resize(endTotalBin-startBin);
      count.resize(endTotalBin-startBin);

      hits.assign(itPulsesTimeBegin+startBin,itPulsesTimeBegin+endTotalBin);
      energy.assign(itIntegralTimeBegin+startBin,itIntegralTimeBegin+endTotalBin);
      event->AddWaveforms(hits,energy);
      for (pmt=0; pmt < fgkNumPMTs; ++pmt) {
        count.assign(fPMTWaveformCount.at(pmt).begin()+startBin,fPMTWaveformCount.at(pmt).begin()+endTotalBin);
        if (std::accumulate(count.begin(),count.end(),0) == 0) {
          continue;
        }
        energy.assign(fPMTWaveform.at(pmt).begin()+startBin,fPMTWaveform.at(pmt).begin()+endTotalBin);
        event->AddPMTWaveforms(pmt,count,energy);
      }

      event->SetLength(endTime-startTime);
      event->SetLargestPMTFraction(largestPMTFraction);

      event->SetIntegralCoated(promptCoated,true);
      event->SetIntegralCoated(totalCoated,false);
      event->SetIntegralUncoated(promptUncoated,true);
      event->SetIntegralUncoated(totalUncoated,false);

      event->SetXPosition(pos.X());
      event->SetYPosition(pos.Y());
      event->SetZPosition(pos.Z());

      event->SetNumVetoTop(vetoActivityTop,false);
      event->SetNumVetoBottom(vetoActivityBottom,false);
      event->SetNumVetoLeft(vetoActivityCLeft,false);
      event->SetNumVetoRight(vetoActivityCRight,false);
      event->SetNumVetoFront(vetoActivityCFront,false);
      event->SetNumVetoBack(vetoActivityCBack,false);

      event->SetNumVetoTop(vetoActivityPromptTop,true);
      event->SetNumVetoBottom(vetoActivityPromptBottom,true);
      event->SetNumVetoLeft(vetoActivityPromptCLeft,true);
      event->SetNumVetoRight(vetoActivityPromptCRight,true);
      event->SetNumVetoFront(vetoActivityPromptCFront,true);
      event->SetNumVetoBack(vetoActivityPromptCBack,true);

      event->SetPMTHits(numPMTs,percentOfPMT);

      MsgDebug(2,"Add Event");
      fEvents->AddSimplifiedEvent(SimplifiedEvent(*event));

      MsgDebug(2,"Reset");
      event->Reset();

      // currently not using the commented out code
      //double avgPercent = std::accumulate(percentOfPMT.begin(),percentOfPMT.end(),0.0);
      //avgPercent /= static_cast<double>(percentOfPMT.size());
      //double sigmaPercent = 0.0;
      //double skewPercent = 0.0;
      //for (const auto & pmt : percentOfPMT) {
      //  sigmaPercent += std::pow(pmt-avgPercent,2.0);
      //  skewPercent += std::pow(pmt-avgPercent,3.0);
      //}
      //skewPercent = skewPercent /
      //              static_cast<double>(percentOfPMT.size()) /
      //              std::pow(1.0/(static_cast<double>(percentOfPMT.size())-1.0)*sigmaPercent,3.0/2.0);
      //sigmaPercent = std::sqrt(1.0/(static_cast<double>(percentOfPMT.size())-1.0)*sigmaPercent);

      MsgDebug(4,Form("End of bin loop, bin = %zu endBin = %d",timeBin,endBin));
      timeBin = endBin;
    } // end threshold if-else
  } // end for size_t timeBin

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMFindEvents::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("TriggerType").Get(fTriggerType);
  c("Threshold").Get(fThreshold);


  MsgInfo("Input parameter values");
  MsgInfo(MsgLog::Form("-TriggerType: %s",fTriggerType.c_str()));
  MsgInfo(MsgLog::Form("-Threshold: %f",fThreshold));

  fIsInit = true;

}

//_______________________________________________________________________________________
void CCMFindEvents::DefineVectors()
{
  fPulsesTime.resize(fgkNumBins,0.f);
  fIntegralTime = fPulsesTime;
  fIntegralDer = fPulsesTime;

  fVetoBottomTime.resize(fgkNumBins,0);
  fVetoTopTime.resize(fgkNumBins,0);
  fVetoCRightTime.resize(fgkNumBins,0);
  fVetoCLeftTime.resize(fgkNumBins,0);
  fVetoCFrontTime.resize(fgkNumBins,0);
  fVetoCBackTime.resize(fgkNumBins,0);
  fVetoTotalTime.resize(fgkNumBins,0);

  fPMTWaveform.resize(fgkNumPMTs,std::vector<float>());
  for (auto & waveform : fPMTWaveform) {
    waveform = fPulsesTime;
  }

  fPMTWaveformCount.resize(fgkNumPMTs,std::vector<int>());
  for (auto & waveform : fPMTWaveformCount) {
    waveform.resize(fgkNumBins,0);
  }
}

//_______________________________________________________________________________________
CCMResult_t CCMFindEvents::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Num Triggers %ld passed trigger type cut",fNumTriggers));

  return kCCMSuccess;
}

