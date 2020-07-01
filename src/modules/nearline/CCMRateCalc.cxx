/*!**********************************************
 * \file CCMRateCalc.cc
 * \brief Main function for the NearlineDiag executable
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date February 24, 2020
 * 
 * Main code for the NearlineDiag executable.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMRateCalc.h"
#include "CCMModuleTable.h"

#include "RawData.h"
#include "Pulses.h"
#include "Events.h"
#include "PMTInfoMap.h"
#include "SinglePulse.h"
#include "PMTInformation.h"
#include "SimplifiedEvent.h"
#include "MsgLog.h"
#include "Utility.h"
#include "NearlineSPEDiag.h"

#include "TSystemFile.h"
#include "TSQLResult.h"
#include "TSQLRow.h"
#include "TSQLServer.h"

#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <map>
#include <locale>

//See CCMModuleTable for info
MODULE_DECL(CCMRateCalc);

//_______________________________________________________________________________________
CCMRateCalc::CCMRateCalc(const char* version) 
  : CCMModule("CCMRateCalc"),
    fRawData(nullptr),
    fPulses(nullptr),
    fEvents(nullptr),
    fTriggerType("BEAM"),
    fWriteDBEntry(false),
    fDBHost(""),
    fDBUser(""),
    fDBPwd(""),
    fFirstTriggerTime(0),
    fLastTriggerTime(0),
    fTotalTriggers(0),
    fPreBeamTriggers(0),
    fInBeamIntegral(0.f)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMRateCalc::CCMRateCalc(const CCMRateCalc& clufdr) 
: CCMModule(clufdr),
  fRawData(clufdr.fRawData),
  fPulses(clufdr.fPulses),
  fEvents(clufdr.fEvents),
  fTriggerType(clufdr.fTriggerType),
  fWriteDBEntry(clufdr.fWriteDBEntry),
  fDBHost(clufdr.fDBHost),
  fDBUser(clufdr.fDBUser),
  fDBPwd(clufdr.fDBPwd),
  fSPECount(clufdr.fSPECount),
  fPMTType(clufdr.fPMTType),
  fFirstTriggerTime(clufdr.fFirstTriggerTime),
  fLastTriggerTime(clufdr.fLastTriggerTime),
  fTotalTriggers(clufdr.fTotalTriggers),
  fPreBeamTriggers(clufdr.fPreBeamTriggers),
  fInBeamIntegral(clufdr.fInBeamIntegral)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMRateCalc::~CCMRateCalc()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMRateCalc::ProcessTrigger()
{

  // Check which trigger occured in DAQ window
  // make sure it is strobe
  if (!fRawData->IsTriggerPresent(fTriggerType)) {
    return kCCMDoNotWrite;
  }

  bool shiftTime = true;
  if (fRawData->IsTriggerPresent("BEAM")) {
    shiftTime = false;
  }

  // keep track of the number of triggers/events are looked at
  // for rate calculation
  if (fPulses->GetTriggerTime()*2e-3 - 9.92 > -0.5) {
    ++fTotalTriggers;

    AddPulses();

    if (fEvents->GetNumEvents() > 0 ) {
      AddEvents(shiftTime);
    }
  }

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMRateCalc::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("TriggerType").Get(fTriggerType);

  std::locale loc;
  for (auto & c : fTriggerType) {
    std::toupper(c,loc);
  }

  fIsInit = true;

}

//_______________________________________________________________________________________
void CCMRateCalc::AddPulses()
{

  // grab the time of the trigger
  std::chrono::time_point <std::chrono::system_clock,std::chrono::duration<unsigned int>> tp_seconds (
      std::chrono::duration<unsigned int>(fPulses->GetComputerSecIntoEpoch()));
  std::chrono::system_clock::time_point tp (tp_seconds);
  if (fFirstTriggerTime == 0) {
    fFirstTriggerTime = std::chrono::system_clock::to_time_t(tp);
  }

  fLastTriggerTime = std::chrono::system_clock::to_time_t(tp);


  // loop over each pulse in the DAQ window
  const size_t kNumPulses = fPulses->GetNumPulses();
  for (size_t loc = 0; loc < kNumPulses; ++loc) {
    int key = fPulses->GetKey(loc);
    if (!PMTInfoMap::IsActive(key)) {
      continue;
    }

    auto pmtInfo = PMTInfoMap::GetPMTInfo(key);
    if (!pmtInfo) {
      continue;
    }

    double time = fPulses->GetPulseTime(loc);
    if (std::isnan(time) || std::isinf(time)) {
      continue;
    }

    double threshold = pmtInfo->GetADCThreshold();
    double spe = pmtInfo->GetADCToPE();

    double pulseIntegral = fPulses->GetPulseIntegral(loc);
    if (std::isnan(pulseIntegral) || std::isinf(pulseIntegral)) {
      continue;
    }

    if (pmtInfo->IsVeto() && pulseIntegral > 5) {
      ++fSPECount[key];
      if (pmtInfo->Is1in()) {
        fPMTType[key] = 1;
      } else {
        fPMTType[key] = 0;
      }
    } else {
      if (pulseIntegral < threshold || pulseIntegral < 0 || spe <= 0) {
        continue;
      } 

      if (time*2e-3 - 9.92 <= -1) {
        ++fSPECount[key];
        if (pmtInfo->IsUncoated()){
          fPMTType[key]=2;
        }else{
          fPMTType[key]=3;
        }
        // if PulseIntegral is above threshold and 0 and the spe if not less than 0, add 1 to the fSPECount for the [key]
        // then use the pmtInfo functions Is1in and IsUncoated to determine the type of pmt and define that
        // redefined for every pulse, but this should not be an issue. Will not define a type if there are no spes (deemed acceptable
      }
    }
  } // end for size_t loc = 0
} // end void AddPulses

//_______________________________________________________________________________________
void CCMRateCalc::AddEvents(bool shiftTime)
{
  const size_t kNumEvents = fEvents->GetNumEvents();
  for (size_t e = 0; e < kNumEvents; ++e) {
    auto simplifiedEvent = fEvents->GetSimplifiedEvent(e);

    double st = simplifiedEvent.GetStartTime();// - 9.92;
    double length = simplifiedEvent.GetLength();

    double promptLength = 0.090;
    bool longerThanPrompt = length > promptLength;

    double energy = simplifiedEvent.GetIntegralTank(longerThanPrompt);// - randomRate*promptLength;
    double hits = simplifiedEvent.GetNumCoated(longerThanPrompt);

    if (shiftTime) {
      st -= 9.92;
    }

    if (hits < 3) {
      continue;
    }

    if (st < -1) {
      ++fPreBeamTriggers;
    } else if (st > -0.35) {
      fInBeamIntegral += energy;
    }
  }
}


//_______________________________________________________________________________________
CCMResult_t CCMRateCalc::EndOfJob()
{
  ///////////////////////////////////////////////
  // Save various information to text files
  // and to the database if connection exists
  ///////////////////////////////////////////////

  std::string speCountOutput = "Board\tChannel\tSPEcount\tSPErate\tType\n";
  int digit = 0;
  int channel = 0;
  double speRate = 0.0;

  std::array<int,3> pmtsOfType;
  std::array<int,3> spesPerType;
  std::array<std::string,3> pmtsTypeN;

  for (int pmtc = 0; pmtc<160; ++pmtc) {
    PMTInfoMap::DecodeKey(pmtc,digit,channel);

    if (channel == 15){
      continue;
    }
    if (fSPECount[pmtc] == 0) {
      continue;
    }

    speRate = fSPECount[pmtc]/(8.92*static_cast<double>(fTotalTriggers))*1000.0;
    spesPerType[fPMTType[pmtc]-1] += speRate;

    //MsgInfo(MsgLog::Form("Board %d Channel %d SPE count %d SPE rate %g",digit,channel,fSPECount[pmtc],speRate));                           
    speCountOutput = speCountOutput+Form("%d\t%d\t%d\t%g\t%d\n",digit,channel,fSPECount[pmtc],speRate,fPMTType[pmtc]);

    if (fPMTType[pmtc]==1){
      ++pmtsOfType[0];
      pmtsTypeN[0] = "Veto";
    } else if (fPMTType[pmtc]==2) {
      ++pmtsOfType[1];
      pmtsTypeN[1] = "Uncoated";
    } else if (fPMTType[pmtc]==3) {
      ++pmtsOfType[2];
      pmtsTypeN[2] = "Coated";
    }
  }

  std::array<float,3> speRates;

  for (int type=0;type<3;++type){
    speRates[type] = spesPerType[type]/pmtsOfType[type];
    speCountOutput = speCountOutput+pmtsTypeN[type]+Form("\t%d\t%g\t%d\n",spesPerType[type],speRates[type],pmtsOfType[type]);
    MsgInfo(MsgLog::Form("%s\t%d\t%g kHz\n",pmtsTypeN[type].c_str(),pmtsOfType[type],speRates[type]));
  }

  MsgInfo(MsgLog::Form("Number of pre beam events %g fInBeamIntegral %g total triggers = %zu overlap triggers = %zu",
        fPreBeamTriggers,fInBeamIntegral,fTotalTriggers));

    //  return EXIT_SUCCESS;

    struct tm * timeInfo = localtime(&fFirstTriggerTime);
    char bufferFirstTriggerTime[80];
    strftime(bufferFirstTriggerTime,80,"%F %T",timeInfo);
    timeInfo = localtime(&fLastTriggerTime);
    char bufferLastTriggerTime[80];
    strftime(bufferLastTriggerTime,80,"%F %T",timeInfo);

    MsgInfo(MsgLog::Form("First Time %s Last Time %s",bufferFirstTriggerTime,bufferLastTriggerTime));

    float preBeamRate = fPreBeamTriggers/static_cast<double>(fTotalTriggers)/8.92*1000.;

    if (fWriteDBEntry) {
      std::string insertCommand = "insert into ccmdb.nearline_diag (start_time, end_time, prebeam_event_rate, spe_rate_1in, spe_rate_uncoated, spe_rate_coated) value (";
      insertCommand = Form("%s'%s','%s','%f','%f','%f','%f');",
          insertCommand.c_str(),
          bufferFirstTriggerTime,bufferLastTriggerTime,
          preBeamRate,speRates[0],speRates[1],speRates[2]);

      TSQLServer *db = TSQLServer::Connect(fDBHost.c_str(),fDBUser.c_str(),fDBPwd.c_str());
      MsgInfo(MsgLog::Form("Insert Command = %s",insertCommand.c_str()));
      TSQLResult *res = db->Query(insertCommand.c_str());
      delete res;
      delete db;
    } // end fWriteDBEntry

    return kCCMSuccess;
}

