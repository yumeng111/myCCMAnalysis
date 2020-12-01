/*!**********************************************
 * \file CCMPMTResponse.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to modle the response of the PMTs
 * in the simulation.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMPMTResponse.h"
#include "CCMModuleTable.h"

#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "MCTruth.h"
#include "SinglePulse.h"
#include "Pulses.h"
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
MODULE_DECL(CCMPMTResponse);

//_______________________________________________________________________________________
CCMPMTResponse::CCMPMTResponse(const char* version) 
  : CCMModule("CCMPMTResponse"),
    fMCTruth(nullptr),
    fPulses(nullptr),
    fRD(),
    fMT(),
    fUniform(0,1),
    fSPEWeights()
{
  //Default constructor
  this->SetCfgVersion(version);

  fMT.seed(fRD());
}

//_______________________________________________________________________________________
CCMPMTResponse::CCMPMTResponse(const CCMPMTResponse& clufdr) 
: CCMModule(clufdr),
  fMCTruth(clufdr.fMCTruth),
  fPulses(clufdr.fPulses),
  fRD(),
  fMT(),
  fUniform(0,1),
  fSPEWeights(clufdr.fSPEWeights)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMPMTResponse::~CCMPMTResponse()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMPMTResponse::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting PMTResponse for MC \"Event\"");
  }


  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMPMTResponse::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMPMTResponse::EndOfJob() 
{ 
  return kCCMSuccess;
}

//_______________________________________________________________________________________
bool CCMPMTResponse::PMTQE(double energy, double angle)
{
  if (energy > 4.0 || energy < 2.0) {
    return false;
  }

  //the following will sample [0,1)
  //if you want [0,1] please let me know
  double testval = fUniform(fMT); 

  if (angle < 0) {
    angle = 0.0;
  }

  //0.2 is the QE
  double prob = 0.2*std::sqrt(angle);

  if (testval < prob) {
    return true;
  }

  return false;
}

//_______________________________________________________________________________________
void CCMPMTResponse::FillSPEWeights()
{
  const size_t kMinKey = PMTInfoMap::GetMinKey();
  const size_t kMaxKey = PMTInfoMap::GetMaxKey();

  for (size_t key = kMinKey; key < kMaxKey; ++key) {
    if (!PMTInfoMap::IsActive(key)) {
      continue;
    }

    auto pmt = PMTInfoMap::GetPMTInfo(key);
    if (!pmt) {
      continue;
    }

    if (pmt->IsVeto()) {
      continue;
    }

    double adcToPE = pmt->GetADCToPE();
    int row = pmt->GetRow();
    int col = pmt->GetColumn();

    fSPEWeights.emplace(std::make_pair(row,col),adcToPE);
  } // end for size_t key

  return;
}

//_______________________________________________________________________________________
double CCMPMTResponse::GetADCValue(int row, int col, double initialCharge)
{
  auto itSPEWeight = fSPEWeights.find(std::make_pair(row,col));
  if (itSPEWeight == fSPEWeights.end()) {
    return 0;
  }

  double adcToPE = itSPEWeight->second;
  double error = std::sqrt(adcToPE);
  std::normal_distribution<double> gaus(adcToPE,error);

  double value = gaus(fMT);

  return initialCharge*value;
}


