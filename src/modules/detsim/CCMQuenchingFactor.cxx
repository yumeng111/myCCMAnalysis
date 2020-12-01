/*!**********************************************
 * \file CCMQuenchingFactor.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to model the LiqAr quenching factor.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMQuenchingFactor.h"
#include "CCMModuleTable.h"

#include "MCTruth.h"
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
MODULE_DECL(CCMQuenchingFactor);

//_______________________________________________________________________________________
CCMQuenchingFactor::CCMQuenchingFactor(const char* version) 
  : CCMModule("CCMQuenchingFactor"),
    fMCTruth(nullptr)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMQuenchingFactor::CCMQuenchingFactor(const CCMQuenchingFactor& clufdr) 
: CCMModule(clufdr),
  fMCTruth(clufdr.fMCTruth)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMQuenchingFactor::~CCMQuenchingFactor()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMQuenchingFactor::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting QuenchingFactor for MC \"Event\"");
  }


  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMQuenchingFactor::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMQuenchingFactor::EndOfJob() 
{ 
  return kCCMSuccess;
}


