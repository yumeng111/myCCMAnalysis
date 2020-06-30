/*!**********************************************
 * \file CCMSumWaveforms.cxx
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 * \brief Read in events and apply some cuts
 *
 * Read in previously found events and apply
 * some cuts.
 *
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMSumWaveforms.h"
#include "CCMModuleTable.h"

#include "Pulses.h"
#include "RawData.h"
#include "SinglePulse.h"
#include "MsgLog.h"
#include "PMTInfoMap.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <memory>

#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"

//See CCMModuleTable for info
MODULE_DECL(CCMSumWaveforms);

//_______________________________________________________________________________________
CCMSumWaveforms::CCMSumWaveforms(const char* version) 
  : CCMModule("CCMSumWaveforms"),
    fRawData(nullptr),
    fPulses(nullptr),
    fOutFileName(""),
    fOutfile(nullptr),
    fTimeHists(),
    fTimeHistsInt()
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMSumWaveforms::CCMSumWaveforms(const CCMSumWaveforms& clufdr) 
: CCMModule(clufdr),
  fRawData(clufdr.fRawData),
  fPulses(clufdr.fPulses),
  fOutFileName(clufdr.fOutFileName),
  fOutfile(clufdr.fOutfile),
  fTimeHists(clufdr.fTimeHists),
  fTimeHistsInt(clufdr.fTimeHistsInt)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMSumWaveforms::~CCMSumWaveforms()
{ 
  // destructor
  fTimeHists.clear();
  fTimeHistsInt.clear();
}

//_______________________________________________________________________________________
CCMResult_t CCMSumWaveforms::ProcessEvent()
{
  if (fTimeHists.empty()) {
    for (size_t i=0; i < fRawData->GetNumBoards(); ++i) {
      for (size_t j=0; j < fRawData->GetNumChannels(); ++j) {
        int key = PMTInfoMap::CreateKey(i,j);
        fTimeHists.emplace_back(std::make_shared<TH1D>(Form("timeHist_%d",key),
              Form("Channel %d;Time #mus;Count",key),8000,-9.92,6.08));
        fTimeHistsInt.emplace_back(std::make_shared<TH1D>(Form("timeHistInt_%d",key),
              Form("Channel %d;Time #mus;Integral",key),8000,-9.92,6.08));
      }
    }
  } // end if fTimeHists.empty()

  auto beamTime = fRawData->GetBCMTime();

  int key = 0;
  int pulseTime = 0;
  float pulseTimeF = 0.f;
  float pulseIntegral = 0.f;

  const size_t kNumPulses = fPulses->GetNumPulses();
  for (size_t pulse = 0; pulse < kNumPulses; ++pulse) {
    key = fPulses->GetKey(pulse);

    auto singlePulse = fPulses->GetSinglePulse(pulse);
    if (singlePulse->GetLength()*2.0 < 20) {
      continue;
    }

    pulseTime = singlePulse->GetTime() - beamTime;
    pulseTimeF = static_cast<double>(pulseTime)*2e-3;
    pulseIntegral = singlePulse->GetIntegral();

    fTimeHists.at(key)->Fill(pulseTimeF+1e-4);
    fTimeHistsInt.at(key)->Fill(pulseTimeF+1e-4,pulseIntegral);
  }
  
  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMSumWaveforms::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Coonfiguration file");

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMSumWaveforms::SetupOutFile()
 * \brief Setup the output file that gets saved
 ***********************************************/
void CCMSumWaveforms::SetupOutFile()
{
  if (fOutfile) {
    if (fOutfile->GetName() != fOutFileName) {
      MsgWarning(MsgLog::Form("New outfile name %s is different than old %s, this should not happen, sticking with ld",
            fOutFileName.c_str(),fOutfile->GetName()));
    }
    return;
  }

  if (!fOutFileName.empty()) {
    fOutfile = gROOT->GetFile(fOutFileName.c_str());
    if (!fOutfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fOutFileName.c_str()));
      return;
    }
    fOutfile->cd();
  }
}

//------------------------------------------------------
CCMResult_t CCMSumWaveforms::EndOfJob()
{
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  if (fOutfile != nullptr) {
    fOutfile->cd();
    for (auto & hist : fTimeHists) {
      if (hist != nullptr) {
        hist->Write();
      }
    }
    for (auto & hist : fTimeHistsInt) {
      if (hist != nullptr) {
        hist->Write();
      }
    }
  }

  return kCCMSuccess;
}


