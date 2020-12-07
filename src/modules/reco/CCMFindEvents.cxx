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

#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "SimplifiedEvent.h"
#include "Events.h"
#include "AccumWaveform.h"
#include "MsgLog.h"
#include "Utility.h"

#include "TROOT.h"
#include "TFile.h"

#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <map>
#include <locale>

//See CCMModuleTable for info
MODULE_DECL(CCMFindEvents);

const int gkMaxBinLoc = 7960;

//_______________________________________________________________________________________
CCMFindEvents::CCMFindEvents(const char* version) 
  : CCMModule("CCMFindEvents"),
    fThreshold(0.4),
    fAccumWaveform(nullptr),
    fEvents(nullptr),
    fEventFinderID(kCCMDynamicLengthEventID),
    fAccumWaveformMethodID(kCCMAccumWaveformTriangleID),
    fNumTriggers(0),
    fNumEvents(0),
    fResetEvents(false),
    fFixedLength(0)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMFindEvents::CCMFindEvents(const CCMFindEvents& clufdr) 
: CCMModule(clufdr),
  fThreshold(clufdr.fThreshold),
  fAccumWaveform(clufdr.fAccumWaveform),
  fEvents(clufdr.fEvents),
  fEventFinderID(clufdr.fEventFinderID),
  fAccumWaveformMethodID(clufdr.fAccumWaveformMethodID),
  fNumTriggers(clufdr.fNumTriggers),
  fNumEvents(clufdr.fNumEvents),
  fResetEvents(clufdr.fResetEvents),
  fFixedLength(clufdr.fFixedLength)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMFindEvents::~CCMFindEvents()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMFindEvents::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting FindEvents for Trigger");
  }

  // create new Events object. This will delete
  // if a previous object was created
  if (fResetEvents) {
    fEvents->Reset();
  }

  fEvents->SetEventNumber(fAccumWaveform->GetEventNumber());
  fEvents->SetComputerSecIntoEpoch(fAccumWaveform->GetComputerSecIntoEpoch());
  fEvents->SetComputerNSIntoSec(fAccumWaveform->GetComputerNSIntoSec());

  fEvents->SetBeamTime(fAccumWaveform->GetBeamOffset());
  fEvents->SetBeamIntegral(fAccumWaveform->GetBeamIntegral());
  fEvents->SetBeamLength(fAccumWaveform->GetBeamLength());

  // count trigger
  ++fNumTriggers;

  int bin2 = 0;
  int endBin = 0;
  int startBin = 0;

  for (size_t timeBin = 0; timeBin < gkMaxBinLoc && timeBin < Utility::fgkNumBins; ++timeBin) {
    if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("Passed DAQ time cut %zu",timeBin));
    }
    //if (MsgLog::GetGlobalDebugLevel() >= 2 && fAccumWaveformMethodID == kCCMAccumWaveformTrianglePulseCutID) {
    //  MsgDebug(2,MsgLog::Form("timeBin %zu value %g threshold %g",timeBin,
    //        fAccumWaveform->GetIndex(timeBin,fAccumWaveformMethodID,kCCMIntegralTimeID),fThreshold));
    //}
    if (fAccumWaveform->GetIndex(timeBin,fAccumWaveformMethodID,kCCMIntegralTimeID) < fThreshold) {
      continue;
    }

    if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("Above threshold %zu",timeBin));
    }

    startBin = timeBin;

    // extrapolate the start time of the event with a linear
    // fit to the rising slope of the accumulated pulse
    startBin = ExtrapolateStartTime(startBin);

    if (startBin >= static_cast<int>(Utility::fgkNumBins)) {
      MsgFatal(MsgLog::Form("Start Bin is greater than or equal to total number of bins: event %ld bin %zu startBin %d",
            fEvents->GetEventNumber(),timeBin,startBin));
    } else if (MsgLog::GetGlobalDebugLevel() >= 3) {
        MsgDebug(3,MsgLog::Form("event %ld bin %zu start bin = %d",fEvents->GetEventNumber(),timeBin,startBin));
    }

    // find end of event: defined as 10 consecutive empty bins
    endBin = Utility::fgkNumBins-1;
    if (fEventFinderID == kCCMDynamicLengthEventID) {
      for (bin2 = timeBin; bin2 < Utility::fgkNumBins-10; ++bin2) {
        if (fAccumWaveform->Integrate(bin2,bin2+10,fAccumWaveformMethodID,kCCMPulsesTimeID) == 0.f) {
          endBin = bin2;
          break;
        } // end if no hits for 20 ns
      } // end for size_t bin2 = bin
    } else {
      endBin = std::min(endBin,startBin+fFixedLength);
    }

    // just a check
    endBin = std::min(Utility::fgkNumBins-1,endBin);

    if (endBin >= static_cast<int>(Utility::fgkNumBins)) {
      MsgFatal(MsgLog::Form("event %ld bin %zu endBin %d",fEvents->GetEventNumber(),timeBin,endBin));
    } else if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("event %ld bin %zu start bin = %d and end bin = %d",fEvents->GetEventNumber(),timeBin,startBin,endBin));
    }

    if (endBin <= startBin) {
      MsgFatal(MsgLog::Form("event %ld bin %zu startBin %d endBin %d",fEvents->GetEventNumber(),timeBin,startBin,endBin));
    } else if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("event %ld bin %zu start bin = %d and end bin = %d",fEvents->GetEventNumber(),timeBin,startBin,endBin));
    }

    // now that we have the start and end of the event grab
    // the necessary information and save it in the Events
    // class container
    SaveEvent(startBin,endBin);

    if (MsgLog::GetGlobalDebugLevel()) {
      MsgDebug(4,Form("End of bin loop, bin = %zu endBin = %d",timeBin,endBin));
    }
    timeBin = endBin;
  } // end for size_t timeBin

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMFindEvents::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("Threshold").Get(fThreshold);
  c("ResetEvents").Get(fResetEvents);

  std::string tempString = "";
  c("EventFinderID").Get(tempString);
  fEventFinderID = Utility::ConvertStringToCCMEventFinderID(tempString);
  c("AccumWaveformMethodID").Get(tempString);
  fAccumWaveformMethodID = Utility::ConvertStringToCCMAccumWaveformMethod(tempString);

  MsgInfo("Input parameter values");
  MsgInfo(MsgLog::Form("-Threshold: %f",fThreshold));
  MsgInfo(MsgLog::Form("-ResetEvents: %d",fResetEvents));
  MsgInfo(MsgLog::Form("-EventFinderID: %s",Utility::ConvertCCMEventFinderIDToString(fEventFinderID).c_str()));
  MsgInfo(MsgLog::Form("-AccumWaveformMethodID: %s",Utility::ConvertCCMAccumWaveformMethodToString(fAccumWaveformMethodID).c_str()));

  if (fEventFinderID == kCCMFixedLengthEventID) {
    c("FixedLength").Get(fFixedLength);
    MsgInfo(MsgLog::Form("-FixedLength: %d",fFixedLength));
  }

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMFindEvents::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Num of Events %ld in %ld triggers with EventFinder %s and AccumWFMethod %s",
        fNumEvents,fNumTriggers, Utility::ConvertCCMEventFinderIDToString(fEventFinderID).c_str(),
        Utility::ConvertCCMAccumWaveformMethodToString(fAccumWaveformMethodID).c_str()));

  return kCCMSuccess;
}

//_______________________________________________________________________________________
int CCMFindEvents::ExtrapolateStartTime(int startBin) 
{
  // find where the first peak is in the event
  int peakLoc = startBin;
  for (int bin2 = startBin; bin2 < Utility::fgkNumBins-1; ++bin2) {
    if (fAccumWaveform->GetIndex(bin2+1,fAccumWaveformMethodID,kCCMIntegralTimeID) < 
        fAccumWaveform->GetIndex(bin2,fAccumWaveformMethodID,kCCMIntegralTimeID)) {
      peakLoc = bin2;
      break;
    } // end if integralTime->at(bin2+1) < ...
  } // end for size_t bin2 = startBin
  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,MsgLog::Form("event %ld Peak Location = %d",fEvents->GetEventNumber(),peakLoc));
  }

  // recalculte the start time based on the rise slope of the event
  if ((peakLoc+startBin)/2-startBin <= 3) {
    return startBin;
  }

  std::vector<float> xPoints;
  std::vector<float> yPoints;

  float c0 = 0;
  float c1 = 0;

  for (int bin2 = startBin; bin2 < (peakLoc+startBin)/2; ++bin2) {
    if (MsgLog::GetGlobalDebugLevel() >= 6) {
      MsgDebug(6,MsgLog::Form("event %ld Peak location loop bin %d",fEvents->GetEventNumber(),bin2));
    }
    xPoints.push_back(bin2);
    yPoints.push_back(fAccumWaveform->GetIndex(bin2,fAccumWaveformMethodID,kCCMIntegralTimeID));
  }
  Utility::LinearUnweightedLS(xPoints.size(),&xPoints.front(),&yPoints.front(),c0,c1);
  if (MsgLog::GetGlobalDebugLevel() >= 6) {
    MsgDebug(6,"Did LinearUnweightedLS");
  }

  // find the first non-empty start bin
  // important for calculating position
  // and integrals
  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,MsgLog::Form("Event %ld Fit results: -c0 = %.3f c1 = %.3f start = %d",fEvents->GetEventNumber(),-c0,c1,static_cast<int>(-c0/c1)));
  }
  if (-c0/c1 > 0 && -c0/c1 < startBin) {
    startBin = static_cast<int>(-c0/c1);
    while (fAccumWaveform->GetIndex(startBin,fAccumWaveformMethodID,kCCMIntegralTimeID) == 0) {
      ++startBin;
    }
    if (MsgLog::GetGlobalDebugLevel() >= 6) {
      MsgDebug(6,MsgLog::Form("Final startBin = %d",startBin));
    }
  }

  xPoints.clear();
  yPoints.clear();

  return startBin;
} // end int CCMFindEvents::ExtrapolateStartTime

//-------------------------------------------------------------------------------------------------
void CCMFindEvents::SaveEvent(int startBin, int endBin)
{
  std::unique_ptr<SimplifiedEvent> event = std::make_unique<SimplifiedEvent>();

  // declear variables used in the for loops
  TVector3 pos;
  auto pmtInfo = PMTInfoMap::GetPMTInfo(0);

  float largestPMTFraction = 0.0;
  float prompt90 = 0;
  float promptCoated = 0.0;
  float promptFit = 0;
  float promptUncoated = 0.0;
  float time = 0;
  float timeCoated = 0.0;
  float timeUncoated = 0.0;
  float total = 0;
  float totalCoated = 0.0;
  float totalUncoated = 0.0;
  float vetoTimeBottom = 0;
  float vetoTimeCBack = 0;
  float vetoTimeCFront = 0;
  float vetoTimeCLeft = 0;
  float vetoTimeCRight = 0;
  float vetoTimeTop = 0;

  int current = 0;
  int end90nsBin = 0;
  int endTotalBin = 0;
  int maxVeto = 0;
  int maxVetoEnd = 0;
  int maxVetoPrompt = 0;
  int maxVetoPromptEnd = 0;
  int maxVetoPromptStart = 0;
  int maxVetoStart = 0;
  int numPMTs = 0;
  int pmt = 0;
  int start = 0;
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

  size_t maxLoc = 0;
  float maxValue = 0.f;

  std::vector<float> energy;
  std::vector<float> hits;
  std::vector<std::pair<int,float>> percentOfPMT;
  std::vector<int> count;
  //int numBins20ns = 20.0/Utility::fgkBinWidth;
  int numBins90ns = 90.0/Utility::fgkBinWidth;
  //int numBins30ns = 30.0/Utility::fgkBinWidth;
  //int numBins1p6us = 1.6e3/Utility::fgkBinWidth;
  //int prevStart = -1000;

  // calculate start and end time of event
  float startTime = Utility::ShiftTime(startBin,fAccumWaveform->GetBeamOffset());
  float endTime = Utility::ShiftTime(endBin,fAccumWaveform->GetBeamOffset());

  start = std::max(startBin - numBins90ns,0);
  end90nsBin = std::min(startBin+numBins90ns,Utility::fgkNumBins-1);;
  endTotalBin = std::min(endBin,Utility::fgkNumBins-1);

  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,MsgLog::Form("end90nsBin = %d endTotalBin = %d",end90nsBin,endTotalBin));
  }


  maxVeto = 0;
  maxVetoStart = 0;
  maxVetoEnd = 0;
  maxVetoPrompt = 0;
  maxVetoPromptStart = 0;
  maxVetoPromptEnd = 0;
  vetoIntLength = std::min(numBins90ns,endBin-start);
  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,MsgLog::Form("Length of veto integral = %d",vetoIntLength));
  }
  for (vetoBin = start; vetoBin+vetoIntLength <= endBin; ++vetoBin) {
    current = fAccumWaveform->Integrate(vetoBin,vetoBin+vetoIntLength,fAccumWaveformMethodID,kCCMVetoTotalTimeID);
    if (current > maxVeto) {
      maxVeto = current;
      maxVetoStart = vetoBin;
      maxVetoEnd = vetoBin+vetoIntLength;
    }
    if (vetoBin <= startBin) {
      current = fAccumWaveform->Integrate(vetoBin,vetoBin+vetoIntLength,fAccumWaveformMethodID,kCCMVetoTotalTimeID);
      if (current > maxVetoPrompt) {
        maxVetoPrompt = current;
        maxVetoPromptStart = vetoBin;
        maxVetoPromptEnd = vetoBin+vetoIntLength;
      }
    }
  }

  if (MsgLog::GetGlobalDebugLevel() >= 3) {
    MsgDebug(3,MsgLog::Form("maxVetoStart %d maxVetoEnd %d",maxVetoStart,maxVetoEnd));
  }
  vetoActivityTop = fAccumWaveform->Integrate(maxVetoStart,maxVetoEnd,fAccumWaveformMethodID,kCCMVetoTopTimeID);
  vetoActivityCRight = fAccumWaveform->Integrate(maxVetoStart,maxVetoEnd,fAccumWaveformMethodID,kCCMVetoCRightTimeID);
  vetoActivityCLeft = fAccumWaveform->Integrate(maxVetoStart,maxVetoEnd,fAccumWaveformMethodID,kCCMVetoCLeftTimeID);
  vetoActivityCFront = fAccumWaveform->Integrate(maxVetoStart,maxVetoEnd,fAccumWaveformMethodID,kCCMVetoCFrontTimeID);
  vetoActivityCBack = fAccumWaveform->Integrate(maxVetoStart,maxVetoEnd,fAccumWaveformMethodID,kCCMVetoCBackTimeID);
  vetoActivityBottom = fAccumWaveform->Integrate(maxVetoStart,maxVetoEnd,fAccumWaveformMethodID,kCCMVetoBottomTimeID);

  if (MsgLog::GetGlobalDebugLevel() >= 3) {
    MsgDebug(3,"Veto prompt integral");
  }

  vetoActivityPromptTop = fAccumWaveform->Integrate(maxVetoPromptStart,maxVetoPromptEnd,fAccumWaveformMethodID,
      kCCMVetoTopTimeID);
  vetoActivityPromptCRight = fAccumWaveform->Integrate(maxVetoPromptStart,maxVetoPromptEnd,fAccumWaveformMethodID,
      kCCMVetoCRightTimeID);
  vetoActivityPromptCLeft = fAccumWaveform->Integrate(maxVetoPromptStart,maxVetoPromptEnd,fAccumWaveformMethodID,
      kCCMVetoCLeftTimeID);
  vetoActivityPromptCFront = fAccumWaveform->Integrate(maxVetoPromptStart,maxVetoPromptEnd,fAccumWaveformMethodID,
      kCCMVetoCFrontTimeID);
  vetoActivityPromptCBack = fAccumWaveform->Integrate(maxVetoPromptStart,maxVetoPromptEnd,fAccumWaveformMethodID,
      kCCMVetoCBackTimeID);
  vetoActivityPromptBottom = fAccumWaveform->Integrate(maxVetoPromptStart,maxVetoPromptEnd,fAccumWaveformMethodID,
      kCCMVetoBottomTimeID);

  if (MsgLog::GetGlobalDebugLevel() >= 3) {
    MsgDebug(3,"Get first none empty bin for the veto tubes");
  }
  vetoTimeTop = fAccumWaveform->FindFirstNoneEmptyBin(start,endBin,fAccumWaveformMethodID,kCCMVetoTopTimeID);
  vetoTimeCRight =  fAccumWaveform->FindFirstNoneEmptyBin(start,endBin,fAccumWaveformMethodID,kCCMVetoCRightTimeID);
  vetoTimeCLeft =  fAccumWaveform->FindFirstNoneEmptyBin(start,endBin,fAccumWaveformMethodID,kCCMVetoCLeftTimeID);
  vetoTimeCFront =  fAccumWaveform->FindFirstNoneEmptyBin(start,endBin,fAccumWaveformMethodID,kCCMVetoCFrontTimeID);
  vetoTimeCBack =  fAccumWaveform->FindFirstNoneEmptyBin(start,endBin,fAccumWaveformMethodID,kCCMVetoCBackTimeID);
  vetoTimeBottom =  fAccumWaveform->FindFirstNoneEmptyBin(start,endBin,fAccumWaveformMethodID,kCCMVetoBottomTimeID);

  // subtract the beam time to remove the jitter of the beam timing
  if (vetoTimeTop >= 0) {
    vetoTimeTop = Utility::ShiftTime(vetoTimeTop,fAccumWaveform->GetBeamOffset());
  }
  if (vetoTimeCRight >= 0) {
    vetoTimeCRight = Utility::ShiftTime(vetoTimeCRight,fAccumWaveform->GetBeamOffset());
  }
  if (vetoTimeCLeft >= 0) {
    vetoTimeCLeft = Utility::ShiftTime(vetoTimeCLeft,fAccumWaveform->GetBeamOffset());
  }
  if (vetoTimeCFront >= 0) {
    vetoTimeCFront = Utility::ShiftTime(vetoTimeCRight,fAccumWaveform->GetBeamOffset());
  }
  if (vetoTimeCBack >= 0) {
    vetoTimeCBack = Utility::ShiftTime(vetoTimeCBack,fAccumWaveform->GetBeamOffset());
  }
  if (vetoTimeBottom >= 0) {
    vetoTimeBottom = Utility::ShiftTime(vetoTimeBottom,fAccumWaveform->GetBeamOffset());
  }
  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,"Got first none empty bin for the veto tubes");
  }

  //MsgInfo(MsgLog::Form("Max Veto Location(Prompt) %d(%d) Amount(Prompt) %d(%d)",
  //      maxVetoStart,maxVetoPromptStart,maxVeto,maxVetoPrompt));

  fAccumWaveform->Max(maxLoc,maxValue,startBin,endBin,fAccumWaveformMethodID,kCCMIntegralTimeID);

  hits.resize(endTotalBin-startBin);
  energy.resize(endTotalBin-startBin);
  count.resize(endTotalBin-startBin);

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

  if (MsgLog::GetGlobalDebugLevel() >= 6) {
    MsgDebug(6,"Loop over pmts");
  }
  for (pmt=0; pmt < Utility::fgkNumPMTs; ++pmt) {
    pmtInfo = PMTInfoMap::GetPMTInfo(pmt);
    if (!pmtInfo) {
      continue;
    }
    if (pmtInfo->IsVeto()) {
      continue;
    }

    // find integral for the entire event window
    if (MsgLog::GetGlobalDebugLevel() >= 7) {
      MsgDebug(7,Form("PMT %d total",pmt));
    }
    total = fAccumWaveform->Integrate(startBin,endTotalBin,fAccumWaveformMethodID,kCCMPMTWaveformID,pmt);
    //time = Utility::FindFirstNoneEmptyBin<float>(itPMTWaveformBegin,itPMTWaveformBegin+startBin,itPMTWaveformBegin+endTotalBin);
    time = fAccumWaveform->FindFirstNoneEmptyBin(startBin,endTotalBin,fAccumWaveformMethodID,kCCMPMTWaveformID,pmt);

    if (total == 0) {
      continue;
    }

    // find integral for the first 90 ns
    if (MsgLog::GetGlobalDebugLevel() >= 7) {
      MsgDebug(7,Form("PMT %d first 90ns",pmt));
    }
    prompt90= fAccumWaveform->Integrate(startBin,end90nsBin,fAccumWaveformMethodID,kCCMPMTWaveformID,pmt);

    fAccumWaveform->CopyVec<int>(count,startBin,endTotalBin,fAccumWaveformMethodID,kCCMPMTWaveformCountID,pmt);
    fAccumWaveform->CopyVec<float>(energy,startBin,endTotalBin,fAccumWaveformMethodID,kCCMPMTWaveformID,pmt);
    event->AddPMTWaveforms(pmt,count,energy);

    if (prompt90 != 0) {
      ++numPMTs;
      percentOfPMT.push_back(std::make_pair(pmt,prompt90));
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
    if (prompt90) {
      pos += *(pmtInfo->GetPosition())*prompt90*prompt90;
      promptFit += prompt90*prompt90;
    }
  } // end for over all digitizer channels
  pos *= 1.0/promptFit;

  timeUncoated = Utility::ShiftTime(timeUncoated,fAccumWaveform->GetBeamOffset());
  timeCoated = Utility::ShiftTime(timeCoated,fAccumWaveform->GetBeamOffset());

  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,"Set events parameters");
  }

  // keep track which methods were used to find the event since it all gets put into one big vector
  event->SetEventFinderMethod(fEventFinderID);
  event->SetAccumWaveformMethod(fAccumWaveformMethodID);

  event->SetStartTime(startTime);
  event->SetStartTimeUncoated(timeUncoated);
  event->SetStartTimeCoated(timeCoated);
  event->SetStartTimeVetoTop(vetoTimeTop);
  event->SetStartTimeVetoRight(vetoTimeCRight);
  event->SetStartTimeVetoLeft(vetoTimeCLeft);
  event->SetStartTimeVetoFront(vetoTimeCFront);
  event->SetStartTimeVetoBack(vetoTimeCBack);
  event->SetStartTimeVetoBottom(vetoTimeBottom);

  event->SetMaxAccumWaveformTime(Utility::ShiftTime(maxLoc,fAccumWaveform->GetBeamOffset())-startTime);
  event->SetMaxAccumWaveformValue(maxValue);

  event->SetNumCoated(fAccumWaveform->Integrate(startBin,end90nsBin,fAccumWaveformMethodID,kCCMPulsesTimeID),true);
  event->SetNumCoated(fAccumWaveform->Integrate(startBin,endTotalBin,fAccumWaveformMethodID,kCCMPulsesTimeID),false);

  fAccumWaveform->CopyVec<float>(hits,startBin,endTotalBin,fAccumWaveformMethodID,kCCMPulsesTimeID);
  fAccumWaveform->CopyVec<float>(energy,startBin,endTotalBin,fAccumWaveformMethodID,kCCMIntegralTimeID);
  event->AddWaveforms(hits,energy);

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

  event->SetPMTHits(percentOfPMT);

  if (MsgLog::GetGlobalDebugLevel() >= 4) {
    MsgDebug(4,"Add Event");
  }
  fEvents->AddSimplifiedEvent(SimplifiedEvent(*event));
  ++fNumEvents;

  if (MsgLog::GetGlobalDebugLevel()>= 4) {
    MsgDebug(4,"Reset");
  }
  event->Reset();

  // currently not using the commented out code
  //float avgPercent = std::accumulate(percentOfPMT.begin(),percentOfPMT.end(),0.0);
  //avgPercent /= static_cast<float>(percentOfPMT.size());
  //float sigmaPercent = 0.0;
  //float skewPercent = 0.0;
  //for (const auto & pmt : percentOfPMT) {
  //  sigmaPercent += std::pow(pmt-avgPercent,2.0);
  //  skewPercent += std::pow(pmt-avgPercent,3.0);
  //}
  //skewPercent = skewPercent /
  //              static_cast<float>(percentOfPMT.size()) /
  //              std::pow(1.0/(static_cast<float>(percentOfPMT.size())-1.0)*sigmaPercent,3.0/2.0);
  //sigmaPercent = std::sqrt(1.0/(static_cast<float>(percentOfPMT.size())-1.0)*sigmaPercent);

  if (MsgLog::GetGlobalDebugLevel() >= 3) {
    MsgDebug(3,"End of SaveEvent");
  }
}

