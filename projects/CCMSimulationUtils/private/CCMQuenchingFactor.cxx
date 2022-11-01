/*!**********************************************
 * \file CCMQuenchingFactor.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to model the LiqAr quenching factor.
 ***********************************************/

#include <map>
#include <array>
#include <cmath>
#include <locale>
#include <memory>
#include <vector>
#include <iostream>

#include "CCMAnalysis/CCMSimulationUtils/CCMQuenchingFactor.h"

#include "CCMAnalysis/CCMDataStructures/MCTruth.h"
#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMModuleTable.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

#include "TROOT.h"
#include "TFile.h"

//See CCMModuleTable for info
MODULE_DECL(CCMQuenchingFactor);

//_______________________________________________________________________________________
CCMQuenchingFactor::CCMQuenchingFactor(const char* version) 
  : CCMModule("CCMQuenchingFactor"),
    fMCTruth(nullptr),
    fRD(),
    fMT(),
    fUniform(0,1),
    fQF(0.25),
    fTotalHits(0),
    fAfterHits(0)
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
  fUniform(0,1),
  fQF(clufdr.fQF),
  fTotalHits(clufdr.fTotalHits),
  fAfterHits(clufdr.fAfterHits)
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

  fMCTruth->SetQuenchingFactor(fQF);

  fTotalHits += kNumHits;
  for (size_t hit = 0; hit < kNumHits; ++hit) {
    fMCTruth->SetHitQuench(hit, TestQuench(fQF));
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

  c("QF").Get(fQF);

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMQuenchingFactor::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Total Hits %ld After QF %ld",fTotalHits,fAfterHits));
  return kCCMSuccess;
}


//_______________________________________________________________________________________
bool CCMQuenchingFactor::TestQuench(double quench)
{
  //the following will sample [0,1)
  //if you want [0,1] please let me know
  double testval = fUniform(fMT);

  if (testval < quench) {
    ++fAfterHits;
    return true;
  }

  return false;
}

//_______________________________________________________________________________________
void CCMQuenchingFactor::ResetQuenchingFactor(double quench)
{
  fQF = quench;
}

