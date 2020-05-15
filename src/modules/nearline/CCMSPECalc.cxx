/*!**********************************************
 * \file CCMSPECalc.cc
 * \brief Main function for the SPECalc executable
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date February 24, 2020
 * 
 * Main code for the SPECalc executable.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMSPECalc.h"
#include "CCMModuleTable.h"

#include "RawData.h"
#include "Pulses.h"
#include "PMTInfoMap.h"
#include "SinglePulse.h"
#include "PMTInformation.h"
#include "SimplifiedEvent.h"
#include "MsgLog.h"
#include "Utility.h"
#include "NearlineSPEDiag.h"

#include <memory>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <map>
#include <locale>

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
  fSPEFinder->SetClassVec(176);
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
CCMResult_t CCMSPECalc::ProcessEvent()
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


