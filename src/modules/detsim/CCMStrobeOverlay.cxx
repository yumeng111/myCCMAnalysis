/*!**********************************************
 * \file CCMStrobeOverlay.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to modle the response of the PMTs
 * in the simulation.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMStrobeOverlay.h"
#include "CCMModuleTable.h"

#include "CCMRootIO.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"
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
MODULE_DECL(CCMStrobeOverlay);

//_______________________________________________________________________________________
CCMStrobeOverlay::CCMStrobeOverlay(const char* version) 
  : CCMModule("CCMStrobeOverlay"),
    fAccumWaveform(nullptr),
    fStrobeData(nullptr)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMStrobeOverlay::CCMStrobeOverlay(const CCMStrobeOverlay& clufdr) 
: CCMModule(clufdr),
  fAccumWaveform(clufdr.fAccumWaveform),
  fStrobeData(clufdr.fStrobeData)
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

  if (fStrobeData == nullptr) {
    return kCCMSuccess;
  }
  
  if (!fStrobeData->ReadOK()) {
    if (!fStrobeData->AdvanceFile()) {
      MsgFatal("Did not pass enough strobe events");
    }
  }

  double before = fAccumWaveform->Integrate(0,Utility::fgkNumBins,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);

  auto rhs = fStrobeData->GetAccumWaveform();
  double strobe = rhs.Integrate(0,Utility::fgkNumBins,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);
  fAccumWaveform->operator+=(rhs);

  double after = fAccumWaveform->Integrate(0,Utility::fgkNumBins,kCCMAccumWaveformTriangleID,kCCMIntegralTimeID);

  MsgInfo(MsgLog::Form("Before %.2f After %.2f Strobe %.2f",before,after,strobe));

  fStrobeData->Advance();

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
    std::vector<std::string> fileList;
    Utility::IndirectFileList(tempString.c_str(),fileList);

    fStrobeData = std::make_shared<CCMRootIO>();
    fStrobeData->SetParameter("ReadBranches","accumWaveform");
    fStrobeData->SetInFileList(fileList);
    fStrobeData->SetupInputFile();
  }

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMStrobeOverlay::EndOfJob() 
{ 
  return kCCMSuccess;
}


