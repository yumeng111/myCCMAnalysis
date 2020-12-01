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
    fAccumWaveform(nullptr)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMStrobeOverlay::CCMStrobeOverlay(const CCMStrobeOverlay& clufdr) 
: CCMModule(clufdr),
  fAccumWaveform(clufdr.fAccumWaveform)
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


  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMStrobeOverlay::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMStrobeOverlay::EndOfJob() 
{ 
  return kCCMSuccess;
}


