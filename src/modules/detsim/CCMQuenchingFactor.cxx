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
    fMCTruth(nullptr),
    fRD(),
    fMT(),
    fUniform(0,1)    
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMQuenchingFactor::CCMQuenchingFactor(const CCMQuenchingFactor& clufdr) 
: CCMModule(clufdr),
  fMCTruth(clufdr.fMCTruth),
  fRD(),
  fMT(),
  fUniform(0,1)
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

  // first get the number of hits observed by the MC event
  const size_t kNumHits = fMCTruth->GetHitNumber();

  // if the number of hits is equal to 0 then return failure as nothing will be changed
  if (kNumHits == 0) {
    return kCCMFailure;
  }

  double quench = fMCTruth->GetQuenchingFactor();

  for (size_t hit = 0; hit < kNumHits; ++hit) {
    bool quenched = TestQuench(quench);
    
    fMCTruth->SetHitQuench(hit, quenched);
  }

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMQuenchingFactor::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.
  int seed = 0;
  c("RANSEED").Get(seed);
  if (seed == 0) {
    fMT.seed(fRD());
  } else {
    fMT.seed(seed);
  }

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMQuenchingFactor::EndOfJob() 
{ 
  return kCCMSuccess;
}


//_______________________________________________________________________________________
bool CCMQuenchingFactor::TestQuench(double quench)
{
  //the following will sample [0,1)
  //if you want [0,1] please let me know
  double testval = fUniform(fMT);

  if (testval < quench) {
    return true;
  }

  return false;
}

void CCMQuenchingFactor::ResetQuenchingFactor(std::shared_ptr<MCTruth> mcTruth, double quench)
{
  fMCTruth = mcTruth;
  fMCTruth->SetQuenchingFactor(quench);
  ProcessTrigger();
}
