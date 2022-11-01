/*!**********************************************
 * \file CCMStrobeOverlay.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to model the response of the PMTs
 * in the simulation.
 ***********************************************/

#include <map>
#include <array>
#include <cmath>
#include <locale>
#include <memory>
#include <vector>
#include <iostream>

#include "CCMAnalysis/CCMSimulationUtils/CCMStrobeOverlay.h"

#include "CCMAnalysis/CCMDataStructures/AccumWaveform.h"
#include "CCMAnalysis/CCMIO/CCMRootIO.h"
#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMModuleTable.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMUtils/PMTInfoMap.h"
#include "CCMAnalysis/CCMUtils/PMTInformation.h"

#include "TROOT.h"
#include "TFile.h"

//See CCMModuleTable for info
MODULE_DECL(CCMStrobeOverlay);

//_______________________________________________________________________________________
CCMStrobeOverlay::CCMStrobeOverlay(const char* version) 
  : CCMModule("CCMStrobeOverlay"),
    fAccumWaveform(nullptr),
    fStrobeData(nullptr),
    fTotalTriggers(0),
    fNumOverLap(0)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMStrobeOverlay::CCMStrobeOverlay(const CCMStrobeOverlay& clufdr) 
: CCMModule(clufdr),
  fAccumWaveform(clufdr.fAccumWaveform),
  fStrobeData(clufdr.fStrobeData),
  fTotalTriggers(clufdr.fTotalTriggers),
  fNumOverLap(clufdr.fNumOverLap)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMStrobeOverlay::~CCMStrobeOverlay()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMStrobeOverlay::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting StrobeOverlay for MC \"Event\"");
  }

  double before = fAccumWaveform->Integrate(0,Utility::fgkNumBins,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
  if (before == 0) {
    return kCCMFailure;
  }

  if (fStrobeData == nullptr) {
    return kCCMSuccess;
  }
  fStrobeData->GoToRandom();
  auto rhs = fStrobeData->GetAccumWaveform();
  double strobe = rhs.Integrate(0,Utility::fgkNumBins,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
  if (strobe == 0) {
    MsgWarning("Strobe waveform is empty");
    return kCCMSuccess;
  }

  /*
  bool overLap = false;
  int triggerTime = std::max(0,static_cast<int>(fAccumWaveform->GetTriggerTime()));
  for (size_t bin = triggerTime; bin < Utility::fgkNumBins; ++bin) {
    auto simContent = fAccumWaveform->GetIndex(bin,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
    if (simContent == 0) {
      continue;
    }

    auto strobeContent = rhs.GetIndex(bin,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
    if (strobeContent != 0) {
      overLap = true;
      break;
    }
  }

  ++fTotalTriggers;
  if (overLap) {
    ++fNumOverLap;
    return kCCMFailure;
  }
  */

  fAccumWaveform->operator+=(rhs);
  //double after = fAccumWaveform->Integrate(0,Utility::fgkNumBins,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
  //MsgInfo(MsgLog::Form("before %.2f strobe %.2f after %.2f",before,strobe,after));

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMStrobeOverlay::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  std::string tempString = "";
  c("FileList").Get(tempString);

  if (!tempString.empty()) {
    std::vector<std::string> fileList = Utility::IndirectFileList(tempString.c_str());

    fStrobeData = std::make_shared<CCMRootIO>();
    fStrobeData->SetParameter("ReadBranches","accumWaveform");
    fStrobeData->SetInFileList(fileList);
    fStrobeData->SetupInputFile();
    std::string env = std::getenv("CCMINSTALL");
    std::string fileName = env + "/build/txtFiles/strobeEventNumbers.txt"; 
    double totalEvents = fStrobeData->GetNumOfEvents(fileName);
    MsgInfo(MsgLog::Form("Number of total Strobe Events = %.0f",totalEvents));
  }

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMStrobeOverlay::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Total Trigger %ld Num Over Lap %ld",fTotalTriggers,fNumOverLap));
  return kCCMSuccess;
}


