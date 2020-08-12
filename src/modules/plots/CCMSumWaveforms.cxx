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
#include "PMTInformation.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <memory>
#include <iterator>

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
CCMResult_t CCMSumWaveforms::ProcessTrigger()
{
  if (fTimeHists.empty()) {
    fTimeHists.emplace_back(std::make_shared<TH1D>("timeHist_CCM",
          "CCM;Time (ns);Count",Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime));
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_CCM")));
    fTimeHistsInt.back()->SetTitle("CCM");
    fTimeHistsInt.back()->SetYTitle("Integral");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_CCMVeto")));
    fTimeHists.back()->SetTitle("CCMVeto");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_CCMVeto")));
    fTimeHistsInt.back()->SetTitle("CCMVeto");
    fTimeHistsInt.back()->SetYTitle("Integral");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015A")));
    fTimeHists.back()->SetTitle("EJ3015A");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015A")));
    fTimeHistsInt.back()->SetTitle("EJ3015A");
    fTimeHistsInt.back()->SetYTitle("Integral");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ301B")));
    fTimeHists.back()->SetTitle("EJ301B");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ301B")));
    fTimeHistsInt.back()->SetTitle("EJ301B");
    fTimeHistsInt.back()->SetYTitle("Integral");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015C_FP3")));
    fTimeHists.back()->SetTitle("EJ3015C_FP3");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015C_FP3")));
    fTimeHistsInt.back()->SetTitle("EJ3015C_FP3");
    fTimeHistsInt.back()->SetYTitle("Integral");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015B")));
    fTimeHists.back()->SetTitle("EJ3015B");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015B")));
    fTimeHistsInt.back()->SetTitle("EJ3015B");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ301A")));
    fTimeHists.back()->SetTitle("EJ301A");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ301A")));
    fTimeHistsInt.back()->SetTitle("EJ301A");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsInt.back()->SetYTitle("Integral");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015D")));
    fTimeHists.back()->SetTitle("EJ3015D");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015D")));
    fTimeHistsInt.back()->SetTitle("EJ3015D");
    fTimeHistsInt.back()->SetYTitle("Integral");
  } // end if fTimeHists.empty()

  auto beamTime = fRawData->GetBCMTime();

  int key = 0;
  int pulseTime = 0;
  float pulseTimeF = 0.f;
  float pulseIntegral = 0.f;
  float pe = 0;
  float threshold = 5;

  const size_t kNumPulses = fPulses->GetNumPulses();
  for (size_t pulse = 0; pulse < kNumPulses; ++pulse) {
    key = fPulses->GetKey(pulse);
    if (!PMTInfoMap::IsActive(key) && key < 169) {
      continue;
    }

    auto singlePulse = fPulses->GetSinglePulse(pulse);
    pulseTime = singlePulse->GetTime();
    if (pulseTime > 7960) {
      continue;
    }

    if (singlePulse->GetLength()*2.0 < 20) {
      continue;
    }

    bool isVeto = false;
    if (key < 160) {
      auto pmtInfo = PMTInfoMap::GetPMTInfo(key);
      if (!pmtInfo) {
        continue;
      }

      isVeto = pmtInfo->IsVeto();

      pulseTimeF = Utility::ShiftTime(pulseTime,beamTime);
      if (isVeto) {
        threshold = 5.0;
        pe = 1.0;
      } else {
        pe = pmtInfo->GetADCToPE();
        if (pe < 0) {
          continue;
        }
        threshold = pmtInfo->GetADCThreshold();
      }

    } else {
      pulseTime = pulseTime - beamTime;
      pulseTimeF = static_cast<double>(pulseTime)*Utility::fgkBinWidth;
      threshold = 5;
      pe = 1.0;
    }

    pulseIntegral = singlePulse->GetIntegral();
    if (pulseIntegral < threshold) {
      continue;
    }

    pulseIntegral /= pe;

    auto it = fTimeHists.begin();
    auto itInt = fTimeHistsInt.begin();

    if (isVeto) {
      std::advance(it,1);
      std::advance(itInt,1);
    } else if (key >= 169 && key <= 174) {
      std::advance(it,key-167);
      std::advance(itInt,key-167);
    } else if (key > 160) {
      continue;
    }

    //int bin = (-Utility::fgkWindowStartTime+pulseTimeF)/Utility::fgkBinWidth+1;

    //MsgInfo(MsgLog::Form("pulseTimeF %.0f bin %d FindBin %d",pulseTimeF,bin,(*it)->FindBin(pulseTimeF+1e-4)));

    //if (bin > Utility::fgkNumBins) {
    //  bin = Utility::fgkNumBins+1;
    //}

    (*it)->Fill(pulseTimeF+1e-4);
    (*itInt)->Fill(pulseTimeF+1e-4,pulseIntegral);
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


