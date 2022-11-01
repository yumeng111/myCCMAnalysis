/*!**********************************************
 * \file CCMSPECalc.cc
 * \brief Main function for the SPECalc executable
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date February 24, 2020
 * 
 * Main code for the SPECalc executable.
 ***********************************************/

#include <map>
#include <array>
#include <cmath>
#include <locale>
#include <memory>
#include <vector>
#include <iostream>

#include "CCMAnalysis/CCMNearline/CCMSPECalc.h"

#include "CCMAnalysis/CCMDataStructures/Pulses.h"
#include "CCMAnalysis/CCMDataStructures/RawData.h"
#include "CCMAnalysis/CCMDataStructures/SinglePulse.h"
#include "CCMAnalysis/CCMDataStructures/SimplifiedEvent.h"
#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMModuleTable.h"
#include "CCMAnalysis/CCMSPE/NearlineSPEDiag.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMUtils/PMTInfoMap.h"
#include "CCMAnalysis/CCMUtils/PMTInformation.h"
#include "CCMAnalysis/CCMUtils/MakeUniquePatch.h"

//See CCMModuleTable for info
MODULE_DECL(CCMSPECalc);

//_______________________________________________________________________________________
CCMSPECalc::CCMSPECalc(const char* version) 
  : CCMModule("CCMSPECalc"),
    fTriggerType("BEAM"),
    fCalibrationFile(""),
    fRedoLEDCalib(false),
    fDAQWindowStart(-9.0), // mus
    fDAQWindowEnd(-1.0), // mus
    fSaveParameters("")
{
  //Default constructor
  this->SetCfgVersion(version);

  fSPEFinder = std::make_unique<NearlineSPEDiag>();
  // TODO Make this more generic (should not assue just one more board)
  fSPEFinder->SetClassVec(Utility::fgkNumDigitizers*Utility::fgkNumChannels);
  fSPEFinder->CreatePEHists();
}

//_______________________________________________________________________________________
CCMSPECalc::CCMSPECalc(const CCMSPECalc& clufdr) 
: CCMModule(clufdr),
  fTriggerType(clufdr.fTriggerType),
  fCalibrationFile(clufdr.fCalibrationFile),
  fRedoLEDCalib(clufdr.fRedoLEDCalib),
  fDAQWindowStart(clufdr.fDAQWindowStart),
  fDAQWindowEnd(clufdr.fDAQWindowEnd),
  fSaveParameters(clufdr.fSaveParameters)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMSPECalc::~CCMSPECalc()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMSPECalc::ProcessTrigger()
{

  if (!fRedoLEDCalib) {
    // Check which trigger occured in DAQ window
    // make sure it is strobe
    if (!fRawData->IsTriggerPresent(fTriggerType)) {
      return kCCMDoNotWrite;
    }

    fSPEFinder->FillPulses(*fPulses,fDAQWindowStart,fDAQWindowEnd);
  }


  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMSPECalc::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("TriggerType").Get(fTriggerType);

  c("CalibrationFile").Get(fCalibrationFile);

  c("RedoLEDCalib").Get(fRedoLEDCalib);
  c("DAQWindowStart").Get(fDAQWindowStart);
  c("DAQWindowEnd").Get(fDAQWindowEnd);
  c("SaveParameters").Get(fSaveParameters);

  fIsInit = true;

}

//_______________________________________________________________________________________
CCMResult_t CCMSPECalc::EndOfJob()
{
  bool isLED = fTriggerType.compare("LED") == 0;
  if (fRedoLEDCalib) {
    fSPEFinder->GetHistsToAdjust(fCalibrationFile,isLED,fDAQWindowStart,fDAQWindowEnd);
    return kCCMSuccess;
  }

  //can activate the above line to have the program analyze the histograms in that file 
  //and produce fits and such in the similarly name pdf.
  fSPEFinder->CalculateRates();
  fSPEFinder->FitPEHists();

/*
  list of possible parameters in fSaveParameters
  rate2DHist = false;
  sPE2DHist = false;
  pEHists = true;
  aDC2DHist = false;
  flangeHists = false;
  sPE1DHist = true;
  rate1DHist = true;
  noisePeakHist = false;
  noiseEndHist = false;
  fitEndHist = false;
  fitRMSHist = false;
  noiseIntegralHist = false;
  tailIntegralHist = false;
  rootTree = true;
  //output file options below
  rootOnly = false;
  pdfOnly = false;
  rootAndPDF = true;
  */

  std::string localCopy = fCalibrationFile;
  size_t loc = localCopy.find(".");
  if (loc != std::string::npos) {
    localCopy.erase(localCopy.begin()+loc,localCopy.end());
  }
  fSPEFinder->FillHist(localCopy, fSaveParameters);

  PMTInfoMap::LoadCalibrationFile(fCalibrationFile);

  return kCCMSuccess;

}


