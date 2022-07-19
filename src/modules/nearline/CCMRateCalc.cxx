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
    fWriteDBEntry(0),
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
  if (Utility::ShiftTime(fPulses->GetTriggerTime(),0,false)> -500) {
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

  c("WriteDBEntry").Get(fWriteDBEntry);
  if (fWriteDBEntry) {
    c("DBHost").Get(fDBHost);
    c("DBUser").Get(fDBUser);
    c("DBPwd").Get(fDBPwd);
  }

  fIsInit = true;
  std::fill(std::begin(fPMTType),std::end(fPMTType),-1);
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

    // TODO Change the veto threshold / make sure these are defined
    if (pmtInfo->IsVeto() && pulseIntegral > threshold) {
      ++fSPECount[key];
      fPMTType[key] = 1;
    } else {
      if (pulseIntegral < threshold || pulseIntegral < 0 || spe <= 0) {
        continue;
      } 

      // TODO Change this time shift
      if (Utility::ShiftTime(time,0,false) <= -1000) {
        ++fSPECount[key];
        if (pmtInfo->IsUncoated()) {
          fPMTType[key]=2;
        } else{
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

    double st = simplifiedEvent.GetStartTime();

    double energy = simplifiedEvent.GetIntegralTank(true);
    double hits = simplifiedEvent.GetNumCoated(true);

    if (hits < 3) {
      continue;
    }

    if (st < -1000) {
      ++fPreBeamTriggers;
    } else if (st > -350) {
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

  //fill array with all zeros
  std::fill(std::begin(pmtsOfType),std::end(pmtsOfType),0);
  std::fill(std::begin(spesPerType),std::end(spesPerType),0);

  const size_t kNumChannels = fSPECount.size();
  for (size_t pmtc = 0; pmtc<kNumChannels; ++pmtc) {
    // only look at channels that correspond to veto, or tank PMTs
    auto pmtInfo = PMTInfoMap::GetPMTInfo(pmtc);
    if (!pmtInfo) {
      continue;
    }

    if (fSPECount.at(pmtc) == 0) {
      continue;
    }

    PMTInfoMap::DecodeKey(pmtc,digit,channel);

    // R. T. Thornton - 9/21/2020
    // the following works since 8.92 is in us and 1000 is converting it to 
    // whichever Hz order it is suppose to be in

    /* This worked for CCM120
    speRate = fSPECount.at(pmtc)/(8.92*static_cast<double>(fTotalTriggers))*1000.0;
    spesPerType.at(fPMTType.at(pmtc)) += speRate;
    speCountOutput = speCountOutput+Form("%d\t%d\t%d\t%g\t%d\n",digit,channel,(int)fSPECount.at(pmtc),speRate,(int)fPMTType.at(pmtc));
    */
    // TODO Fix this for CCM200
    speRate = fSPECount[pmtc]*1e6/(-(Utility::fgkWindowStartTime + 1000.0)*static_cast<double>(fTotalTriggers));
    spesPerType[fPMTType[pmtc]-1] += speRate;
    if (speRate > 1000){
        MsgInfo(MsgLog::Form("Board %d Channel %d SPE count %d SPE rate %g",digit,channel,fSPECount[pmtc],speRate));
        speCountOutput = speCountOutput+Form("%d\t%d\t%d\t%g\t%d\n",digit,channel,fSPECount[pmtc],speRate,fPMTType[pmtc]);
    }

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

  for (size_t type=0;type<speRates.size();++type){
    speRates.at(type) = spesPerType.at(type)/pmtsOfType.at(type);
    speCountOutput = speCountOutput+pmtsTypeN.at(type)+Form("\t%d\t%g\t%d\n",spesPerType.at(type),speRates.at(type),pmtsOfType.at(type));
    //MsgInfo(MsgLog::Form("%s\t%d\t%g kHz\n",pmtsTypeN.at(type).c_str(),pmtsOfType.at(type),speRates.at(type)));
  }

  MsgInfo(MsgLog::Form("Number of pre beam events %g fInBeamIntegral %g total triggers = %zu overlap triggers = %zu",
        fPreBeamTriggers,fInBeamIntegral,fTotalTriggers,fTotalTriggers));

  MsgInfo(MsgLog::Form("\n%s",speCountOutput.c_str()));

    struct tm * timeInfo = localtime(&fFirstTriggerTime);
    char bufferFirstTriggerTime[80];
    strftime(bufferFirstTriggerTime,80,"%F %T",timeInfo);
    timeInfo = localtime(&fLastTriggerTime);
    char bufferLastTriggerTime[80];
    strftime(bufferLastTriggerTime,80,"%F %T",timeInfo);

    MsgInfo(MsgLog::Form("First Time %s Last Time %s",bufferFirstTriggerTime,bufferLastTriggerTime));
    // TODO Check if this rate calculation is correct
    // float preBeamRate = fPreBeamTriggers/static_cast<double>(fTotalTriggers)/8.92*1000.;
    float preBeamRate = fPreBeamTriggers*1e6/static_cast<double>(fTotalTriggers)/-(Utility::fgkWindowStartTime + 1000);

    if (fWriteDBEntry) {
      MsgWarning("Going to write to data base");
      std::string insertCommand = "insert into ccmdb.nearline_diag (start_time, end_time, prebeam_event_rate, spe_rate_1in, spe_rate_uncoated, spe_rate_coated) value (";
      insertCommand = Form("%s'%s','%s','%f','%f','%f','%f');",
          insertCommand.c_str(),
          bufferFirstTriggerTime,bufferLastTriggerTime,
          preBeamRate,speRates.at(0),speRates.at(1),speRates.at(2));

      TSQLServer *db = TSQLServer::Connect(fDBHost.c_str(),fDBUser.c_str(),fDBPwd.c_str());
      MsgInfo(MsgLog::Form("Insert Command = %s",insertCommand.c_str()));
      TSQLResult *res = db->Query(insertCommand.c_str());
      delete res;
      delete db;
    } // end fWriteDBEntry

    return kCCMSuccess;
}

