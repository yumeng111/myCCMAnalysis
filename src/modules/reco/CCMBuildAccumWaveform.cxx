/*!**********************************************
 * \file CCMBuildAccumWaveform.cxx
 * \author R.T. Thornton (LANL)
 * \date February 24,2020
 *
 * Main code to find events in the detector that
 * that are for a specific trigger
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMBuildAccumWaveform.h"
#include "CCMModuleTable.h"

#include "RawData.h"
#include "Pulses.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "SimplifiedEvent.h"
#include "AccumWaveform.h"
#include "MsgLog.h"
#include "Utility.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"

#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <map>
#include <locale>

//See CCMModuleTable for info
MODULE_DECL(CCMBuildAccumWaveform);

const int gkMaxBinLoc = 5960;//for CCM200.
//const int gkMaxBinLoc = 7960;//for CCM120.

//_______________________________________________________________________________________
CCMBuildAccumWaveform::CCMBuildAccumWaveform(const char* version) 
  : CCMModule("CCMBuildAccumWaveform"),
    fTriggerType("BEAM"),
    fAccumWaveform(nullptr),
    fRawData(nullptr),
    fPulses(nullptr),
    fNumTriggers(0),
    fPulseTimeLowValue(Utility::fgkWindowStartTime),
    fPulseTimeHighValue(Utility::fgkWindowEndTime),
    fBuildMethods(),
    fBeamTime(0)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMBuildAccumWaveform::CCMBuildAccumWaveform(const CCMBuildAccumWaveform& clufdr) 
: CCMModule(clufdr),
  fTriggerType(clufdr.fTriggerType),
  fAccumWaveform(clufdr.fAccumWaveform),
  fRawData(clufdr.fRawData),
  fPulses(clufdr.fPulses),
  fNumTriggers(clufdr.fNumTriggers),
  fPulseTimeLowValue(clufdr.fPulseTimeLowValue),
  fPulseTimeHighValue(clufdr.fPulseTimeHighValue),
  fBuildMethods(clufdr.fBuildMethods),
  fBeamTime(clufdr.fBeamTime)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMBuildAccumWaveform::~CCMBuildAccumWaveform()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMBuildAccumWaveform::ProcessTrigger()
{
  std::cout<<"inside baw"<<"\n";
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting BuildAccumWaveform for Trigger");
  }

  // create new Events object. This will delete if a previous object was created.
  for (auto & method : fBuildMethods) {
    fAccumWaveform->AddMethod(method);
  }
  fAccumWaveform->Reset();

  // Check which trigger occured in DAQ window
  // make sure it is strobe
  bool mcEvent = fTriggerType.find("MC") != std::string::npos;
  if (!mcEvent) {
    if (!fRawData->IsTriggerPresent(fTriggerType)) {
      return kCCMFailure;
    }

    fAccumWaveform->SetEventNumber(fRawData->GetEventNumber());
    fAccumWaveform->SetComputerSecIntoEpoch(fRawData->GetGPSSecIntoDay());
    fAccumWaveform->SetComputerNSIntoSec(fRawData->GetGPSNSIntoSec());
  }


  fBeamTime = Utility::fgkNumBins+1; 

  // If fTriggerType == beam
  // Find when the BCM triggered in the DAQ window 
  fBeamTime = Utility::fgkNumBins+1;
  double beamIntegral = 0;
  double beamLength = 0;
  if (fTriggerType.find("BEAM") != std::string::npos) {
    fBeamTime = fRawData->GetBCMTime(&beamIntegral,&beamLength);
    if (fBeamTime == Utility::fgkNumBins+1) {
      MsgWarning(MsgLog::Form("Event %d No BCM trigger found skipping trigger",fRawData->GetEventNumber()));
      return kCCMFailure;
    }
  } else if (fBeamTime == Utility::fgkNumBins+1) {
    fBeamTime = 0;
  }


  fAccumWaveform->SetBeamOffset(fBeamTime);
  fAccumWaveform->SetBeamIntegral(beamIntegral);
  fAccumWaveform->SetBeamLength(beamLength*Utility::fgkBinWidth);
  fAccumWaveform->SetTriggerTime(fPulses->GetTriggerTime());
  

  // count trigger
  ++fNumTriggers;

  const PMTInformation * pmtInfo;
  double hitFraction = 0;
  double integral = 0;
  double length = 0;
  double pe = 0;
  double pulseIntegral = 0;
  double shiftTime = 0;
  double threshold = 0;
  double time = 0;
  int bin = 0;
  int column = 0;
  int firstBin = 0;
  int key = 0;
  int lastBin = 0;
  int middle = 0;
  std::string name = "";

  const bool kBuildTriangle = 
    std::find(fBuildMethods.begin(),fBuildMethods.end(),kCCMAccumWaveformTriangleID) != fBuildMethods.end();

  const bool kBuildTrianglePulseCut = 
    std::find(fBuildMethods.begin(),fBuildMethods.end(),kCCMAccumWaveformTrianglePulseCutID) != fBuildMethods.end();

  const bool kBuildStart = 
    std::find(fBuildMethods.begin(),fBuildMethods.end(),kCCMAccumWaveformStartID) != fBuildMethods.end();

  const bool kBuildStartPulseCut = 
    std::find(fBuildMethods.begin(),fBuildMethods.end(),kCCMAccumWaveformStartPulseCutID) != fBuildMethods.end();


  // loop through the pulses
  const size_t kNumPulses = fPulses->GetNumPulses();
  for (size_t loc = 0; loc < kNumPulses; ++loc) {

    // make sure the PMT is to be used for the analysis
    key = fPulses->GetKey(loc);
    if (!PMTInfoMap::IsActive(key)) {
      continue;
    }
    pmtInfo = PMTInfoMap::GetPMTInfo(key);
    if (!pmtInfo) {
      continue;
    }
    if((key+1)%16 == 0){
      continue;
    }
    //high rate from Ed.
    if(key == 162){
      continue;
    }
    //dont see any light: LED analysis.
    if(key == 98 || key == 213){
      continue;
    }
    if(key > 254){
      continue;
    }
    //very high threshold(61 ADC)
    if(key == 199){
      continue;
    }


    //make sure the time is within the range of the DAQ window for analysis only include pulses that start before the beam-related background happens
    time = fPulses->GetPulseTime(loc);
    if (std::isnan(time) || std::isinf(time) || time > gkMaxBinLoc) {
      continue;
    }


    //check to see if there is a time cut on which pulses to include only care about the start time of the pulse
    bool passPulseTimeCut = true;
    shiftTime = Utility::ShiftTime(time,fBeamTime);
    if (shiftTime < fPulseTimeLowValue || shiftTime > fPulseTimeHighValue) {
      passPulseTimeCut = false;
    }
    
    //get the integral of the pulse and make sure it is a real number
    pulseIntegral = fPulses->GetPulseIntegral(loc);
    if (std::isnan(pulseIntegral) || std::isinf(pulseIntegral)) {
      continue;
    }

    //get threshold and ADCtoPE for the PMT
    threshold = pmtInfo->GetADCThreshold();
    pe = pmtInfo->GetADCToPE();

    //accounting for high threshold of some PMTs.
    if(key == 86) threshold = 20.0;
    if(key == 88) threshold = 27.0;
    if(key == 90) threshold = 22.0;
    if(key == 93) threshold = 24.0;
    if(key == 97) threshold = 27.0;
    if(key == 102) threshold = 36.0;
    if(key == 113) threshold = 35.0;
    if(key == 196) threshold = 33.0;
    if(key == 197) threshold = 23.0;
    if(key == 203) threshold = 36.0;
    if(key == 206) threshold = 28.0;
    
 
    time = std::max(std::min(time,static_cast<double>(Utility::fgkNumBins-1)),0.0);

    // if the PMT was a veto pmt see if the integral is above 9ADC.
    // keep track the number of times that pmt fired in the DAQ window
    if (pmtInfo->IsVeto()) {
      name = pmtInfo->GetLocName();
      if (pulseIntegral > 9.0) {
	//accounting for higher noise levels of some pmts.
	if(key == 0 && pulseIntegral < 30.0){
	  continue;
	}
	if(key == 1 && pulseIntegral < 37.0){
          continue;
	}
	if(key == 3 && pulseIntegral < 25.0){
          continue;
	}
	if(key == 5 && pulseIntegral < 33.0){
          continue;
	}
	if(key == 10 && pulseIntegral < 36.0){
          continue;
	}
	if(key == 12 && pulseIntegral < 33.0){
          continue;
	}
	if(key == 20 && pulseIntegral < 18.0){
          continue;
	}
	if(key == 22 && pulseIntegral < 34.0){
          continue;
	}
	if(key == 23 && pulseIntegral < 18.0){
          continue;
	}
	if(key == 25 && pulseIntegral < 25.0){
          continue;
	}
	if(key == 30 && pulseIntegral < 18.0){
          continue;
	}

 
        CCMAccWaveform_t waveform = kCCMVetoTotalTimeID;
        // veto top tubes
        if (name.find("VT") != std::string::npos) {
	  waveform = kCCMVetoTopTimeID;
          // veto column top tubes
        } else if (name.find("VC") != std::string::npos) {
          // Back (columes 22 through 3)
          // Left (columes 4 through 9)
          // Right (columes 16 through 21)
          // Front (columes 10 through 15)
          column = pmtInfo->GetColumn();
          if (column >= 4 && column <= 9) {
	     waveform = kCCMVetoCLeftTimeID;
	  }else if (column >= 10 && column <= 15) {
	     waveform = kCCMVetoCFrontTimeID;
	  } else if (column >= 16 && column <= 21) {
	     waveform = kCCMVetoCRightTimeID;
	  } else if (column > 0) {
	     waveform = kCCMVetoCBackTimeID;
	  }
	} else {
	  waveform = kCCMVetoBottomTimeID;
        } // end if-else over where the veto is located

        // there is only one veto PMT accumulated build method
        // save the veto pmt waveform for all build method types
        for (size_t method = 0; method < kCCMAccumWaveformTotalID; ++method) {
          CCMAccumWaveformMethod_t methodID = static_cast<CCMAccumWaveformMethod_t>(method);
          if (std::find(fBuildMethods.begin(),fBuildMethods.end(),methodID) == fBuildMethods.end()) {
            continue;
          }

          if ((methodID == kCCMAccumWaveformStartPulseCutID ||
                methodID == kCCMAccumWaveformTrianglePulseCutID) &&
              !passPulseTimeCut) {
            continue;
          }


          fAccumWaveform->FillIndex(static_cast<int>(time),1.0,methodID,kCCMPMTWaveformCountID,key);
          fAccumWaveform->FillIndex(static_cast<int>(time),1.0,methodID,kCCMPMTWaveformID,key);

          fAccumWaveform->FillIndex(static_cast<int>(time),1.0,methodID,waveform);
          fAccumWaveform->FillIndex(static_cast<int>(time),1.0,methodID,kCCMVetoTotalTimeID);
        }
      } // end if pulseIntegral is greater than 9

      // move to the next event
      continue;
    } // end if pmtInfo->IsVeto()
    

    // check if pulse integral is below threshold or the ADCtoPE calibration is negative or zero
    // if either is true do not analyze the pulse
    if (pulseIntegral < threshold || pe <= 0) {
      continue;
    }


    // check length of pulse make sure it is greater than or equal to 20 ns
    length = fPulses->GetPulseLength(loc);
    if (length*Utility::fgkBinWidth < 20.0) {
      continue;
    }


    // apply the ADCtoPE to the pulse
    pulseIntegral /= pe;


    /////////////////////////////////////////////
    // Fill the accumulated waveform based on
    // the start representation
    // of the pulse
    /////////////////////////////////////////////
    if (kBuildStart) {
      fAccumWaveform->FillIndex(static_cast<int>(time),1,kCCMAccumWaveformStartID,kCCMPMTWaveformCountID,key);
      fAccumWaveform->FillIndex(static_cast<int>(time),pulseIntegral,kCCMAccumWaveformStartID,kCCMPMTWaveformID,key);
      fAccumWaveform->FillIndex(static_cast<int>(time),1.0,kCCMAccumWaveformStartID,kCCMPulsesTimeID);
      fAccumWaveform->FillIndex(static_cast<int>(time),pulseIntegral,kCCMAccumWaveformStartID,kCCMIntegralTimeID);
    }

    if (passPulseTimeCut && kBuildStartPulseCut) {
      fAccumWaveform->FillIndex(static_cast<int>(time),1,kCCMAccumWaveformStartPulseCutID,kCCMPMTWaveformCountID,key);
      fAccumWaveform->FillIndex(static_cast<int>(time),pulseIntegral,kCCMAccumWaveformStartPulseCutID,kCCMPMTWaveformID,key);
      fAccumWaveform->FillIndex(static_cast<int>(time),1.0,kCCMAccumWaveformStartPulseCutID,kCCMPulsesTimeID);
      fAccumWaveform->FillIndex(static_cast<int>(time),pulseIntegral,kCCMAccumWaveformStartPulseCutID,kCCMIntegralTimeID);
    }

    /////////////////////////////////////////////
    // Fill the accumulated waveform based on
    // the triangle representation
    // of the pulse
    /////////////////////////////////////////////
    firstBin = time;
    lastBin  = std::min(time+length,static_cast<double>(Utility::fgkNumBins));
    middle = (lastBin + firstBin)/2.0;
   

    // add the pulse to the integral and pmt histograms
    // reconstruct the pulse for the integral histogram
    // as a triangle where the area of the triangle is
    // equal to to the pulseIntegral

    if (kBuildTriangle) {
      fAccumWaveform->FillIndex(firstBin,1,kCCMAccumWaveformTriangleID,kCCMPMTWaveformCountID,key);
    }
    if (passPulseTimeCut && kBuildTrianglePulseCut) {
      fAccumWaveform->FillIndex(firstBin,1,kCCMAccumWaveformTrianglePulseCutID,kCCMPMTWaveformCountID,key);
    }

    integral = 0.0;
    hitFraction = 0.0;
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
      
 
      if (kBuildTriangle) {
        fAccumWaveform->FillIndex(bin,integral,kCCMAccumWaveformTriangleID,kCCMPMTWaveformID,key);
        fAccumWaveform->FillIndex(bin,hitFraction,kCCMAccumWaveformTriangleID,kCCMPulsesTimeID);
        fAccumWaveform->FillIndex(bin,integral,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
      }

      if (passPulseTimeCut && kBuildTrianglePulseCut) {
        fAccumWaveform->FillIndex(bin,integral,kCCMAccumWaveformTrianglePulseCutID,kCCMPMTWaveformID,key);
        fAccumWaveform->FillIndex(bin,hitFraction,kCCMAccumWaveformTrianglePulseCutID,kCCMPulsesTimeID);
        fAccumWaveform->FillIndex(bin,integral,kCCMAccumWaveformTrianglePulseCutID,kCCMIntegralTimeID);
      }

    } // end for bin = firstBin
  } // for over each pulse

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMBuildAccumWaveform::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("TriggerType").Get(fTriggerType);

  std::string tempString = "";
  c("BuildMethods").Get(tempString);

  if (tempString.find("PulseCut") != std::string::npos) {
    c("PulseTimeLowValue").Get(fPulseTimeLowValue);
    c("PulseTimeHighValue").Get(fPulseTimeHighValue);

    fPulseTimeLowValue = std::max(fPulseTimeLowValue,Utility::fgkWindowStartTime);
    fPulseTimeHighValue = std::min(fPulseTimeHighValue,Utility::fgkWindowEndTime);

  } else {
    fPulseTimeLowValue = Utility::fgkWindowStartTime;
    fPulseTimeHighValue= Utility::fgkWindowEndTime;
  }

  std::stringstream ss(tempString);
  while (ss >> tempString) {
    fBuildMethods.push_back(Utility::ConvertStringToCCMAccumWaveformMethod(tempString));
  }

  MsgInfo("Input parameter values");
  MsgInfo(MsgLog::Form("-TriggerType: %s",fTriggerType.c_str()));
  MsgInfo(MsgLog::Form("-NumBins: %d",Utility::fgkNumBins));
  MsgInfo(MsgLog::Form("-BinWidth: %f",Utility::fgkBinWidth));
  MsgInfo(MsgLog::Form("-PulseTimeLowValue: %f",fPulseTimeLowValue));
  MsgInfo(MsgLog::Form("-PulseTimeHighValue: %f",fPulseTimeHighValue));
  MsgInfo("-Will build accumulated waveforms with the following methods");
  for (auto & name : fBuildMethods) {
    MsgInfo(MsgLog::Form("\t-%s",Utility::ConvertCCMAccumWaveformMethodToString(name).c_str()));
  }

  fIsInit = true;

}

//_______________________________________________________________________________________
CCMResult_t CCMBuildAccumWaveform::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Num Triggers %ld passed trigger type cut",fNumTriggers));

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMBuildAccumWaveform::Dump()
{
  for (auto & method : fBuildMethods) {
    MsgInfo(MsgLog::Form("Method %s",Utility::ConvertCCMAccumWaveformMethodToString(method).c_str()));
  }
  fAccumWaveform->DumpInfo();
}

