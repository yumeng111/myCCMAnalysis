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
  : CCMModule("CCMSumWaveforms"), fRawData(nullptr), fPulses(nullptr), fOutFileName(""),
    fOutfile(nullptr), fTriggerType("ALL"), fDoBCMCut(false), fApplyFP3Offset(false),
    fBCMTimeLowCut(Utility::fgkWindowStartTime), fBCMTimeHighCut(Utility::fgkWindowEndTime),
    fBCMLengthLowCut(0), fBCMLengthHighCut(Utility::fgkNumBins*Utility::fgkBinWidth),
    fBCMIntegralLowCut(0), fBCMIntegralHighCut(std::numeric_limits<double>::max()),
    fTimeHistsHit(),
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
  fTriggerType(clufdr.fTriggerType),
  fDoBCMCut(clufdr.fDoBCMCut),
  fApplyFP3Offset(clufdr.fApplyFP3Offset),
  fBCMTimeLowCut(clufdr.fBCMTimeLowCut), fBCMTimeHighCut(clufdr.fBCMTimeHighCut), 
  fBCMLengthLowCut(clufdr.fBCMLengthLowCut), fBCMLengthHighCut(clufdr.fBCMLengthHighCut), 
  fBCMIntegralLowCut(clufdr.fBCMIntegralLowCut), fBCMIntegralHighCut(clufdr.fBCMIntegralHighCut),
  fTimeHistsHit(clufdr.fTimeHistsHit),
  fTimeHists(clufdr.fTimeHists),
  fTimeHistsInt(clufdr.fTimeHistsInt)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMSumWaveforms::~CCMSumWaveforms()
{ 
  // destructor
  fTimeHistsHit.clear();
  fTimeHists.clear();
  fTimeHistsInt.clear();
}

//_______________________________________________________________________________________
CCMResult_t CCMSumWaveforms::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting SumWaveforms for Trigger");
  }

  if (fTimeHists.empty()) {
    fTimeHists.emplace_back(std::make_shared<TH1D>("timeHist_CCM",
          "CCM;Time (ns);Count",Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime));
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_CCM")));
    fTimeHistsInt.back()->SetTitle("CCM");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_CCM")));
    fTimeHistsHit.back()->SetTitle("CCM");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_CCMVeto")));
    fTimeHists.back()->SetTitle("CCMVeto");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_CCMVeto")));
    fTimeHistsInt.back()->SetTitle("CCMVeto");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_CCMVeto")));
    fTimeHistsHit.back()->SetTitle("CCMVeto");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015A")));
    fTimeHists.back()->SetTitle("EJ3015A");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015A")));
    fTimeHistsInt.back()->SetTitle("EJ3015A");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_EJ3015A")));
    fTimeHistsHit.back()->SetTitle("EJ3015A");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ301B")));
    fTimeHists.back()->SetTitle("EJ301B");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ301B")));
    fTimeHistsInt.back()->SetTitle("EJ301B");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_EJ301B")));
    fTimeHistsHit.back()->SetTitle("EJ301B");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015C_FP3")));
    fTimeHists.back()->SetTitle("EJ3015C_FP3");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015C_FP3")));
    fTimeHistsInt.back()->SetTitle("EJ3015C_FP3");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_EJ3015C_FP3")));
    fTimeHistsHit.back()->SetTitle("EJ3015C_FP3");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015B")));
    fTimeHists.back()->SetTitle("EJ3015B");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015B")));
    fTimeHistsInt.back()->SetTitle("EJ3015B");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_EJ3015B")));
    fTimeHistsHit.back()->SetTitle("EJ3015B");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ301A")));
    fTimeHists.back()->SetTitle("EJ301A");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ301A")));
    fTimeHistsInt.back()->SetTitle("EJ301A");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_EJ301A")));
    fTimeHistsHit.back()->SetTitle("EJ301A");
    fTimeHistsHit.back()->SetYTitle("Hit");

    fTimeHists.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHist_EJ3015D")));
    fTimeHists.back()->SetTitle("EJ3015D");
    fTimeHists.back()->SetYTitle("Count");
    fTimeHistsInt.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistInt_EJ3015D")));
    fTimeHistsInt.back()->SetTitle("EJ3015D");
    fTimeHistsInt.back()->SetYTitle("Integral");
    fTimeHistsHit.emplace_back(dynamic_cast<TH1D*>(fTimeHists.back()->Clone("timeHistHit_EJ3015D")));
    fTimeHistsHit.back()->SetTitle("EJ3015D");
    fTimeHistsHit.back()->SetYTitle("Hit");
  } // end if fTimeHists.empty()

  if (!fRawData->IsTriggerPresent(fTriggerType)) {
    return kCCMSuccess;
  }

  double beamIntegral = 0;
  double beamLength = 0;
  auto beamTime = fRawData->GetBCMTime(&beamIntegral,&beamLength);

  double beamTimeShifted = Utility::ShiftTime(beamTime,0,false);
  beamLength *= Utility::fgkBinWidth;

  if (fDoBCMCut) {
    //if (beamTimeShifted < -345.77-10.54*3.0 || beamTimeShifted > -345.77+10.54*3.0 ||
    //    beamLength < 576.0-23.5*3.0 || beamLength > 576.0+23.5*3.0 ||
    //    beamIntegral < 17811.4 || beamIntegral > 30978.4)
  if (beamTimeShifted < fBCMTimeLowCut || beamTimeShifted > fBCMTimeHighCut ||
      beamLength < fBCMLengthLowCut || beamLength > fBCMLengthHighCut ||
      beamIntegral < fBCMIntegralLowCut || beamIntegral > fBCMIntegralHighCut) {
      return kCCMSuccess;
    }
  } else {
    beamTime = 0;
  }

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

    double length = singlePulse->GetLength();
    if (length*2.0 < 20) {
      continue;
    }

    bool isVeto = false;
    if (key < 160) {
      auto pmtInfo = PMTInfoMap::GetPMTInfo(key);
      if (!pmtInfo) {
        continue;
      }

      isVeto = pmtInfo->IsVeto();

      pulseTimeF = Utility::ShiftTime(pulseTime,beamTime,fApplyFP3Offset);
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
      pulseTimeF = Utility::ShiftTime(pulseTime,beamTime,false);
      threshold = 5;
      pe = 1.0;
    }

    pulseIntegral = singlePulse->GetIntegral();
    if (pulseIntegral < threshold) {
      continue;
    }

    pulseIntegral /= pe;

    auto it = fTimeHists.begin();
    auto itHit = fTimeHistsHit.begin();
    auto itInt = fTimeHistsInt.begin();

    if (isVeto) {
      std::advance(it,1);
      std::advance(itInt,1);
      std::advance(itHit,1);
    } else if (key >= 169 && key <= 174) {
      std::advance(it,key-167);
      std::advance(itInt,key-167);
      std::advance(itHit,key-167);
    } else if (key > 160) {
      continue;
    }

    //int bin = (-Utility::fgkWindowStartTime+pulseTimeF)/Utility::fgkBinWidth+1;

    //MsgInfo(MsgLog::Form("pulseTimeF %.0f bin %d FindBin %d",pulseTimeF,bin,(*it)->FindBin(pulseTimeF+1e-4)));

    //if (bin > Utility::fgkNumBins) {
    //  bin = Utility::fgkNumBins+1;
    //}

    (*it)->Fill(pulseTimeF+1e-4,pulseIntegral);
    (*itHit)->Fill(pulseTimeF+1e-4,1);

    double middle = length/2.0;
    double integral = 0.0;
    for (int bin = 0; bin < length; ++bin) {
      integral = 0.0;
      if (bin < middle) {
        integral = 2.0*pulseIntegral/length*(static_cast<double>(bin)/middle);
      } else if (bin != middle) {
        integral = 2.0*pulseIntegral/length*(static_cast<double>(length - bin)/(length - middle));
      } else {
        integral = 2.0*pulseIntegral/length;
      }

      (*itInt)->Fill(pulseTimeF+1e-4+static_cast<double>(bin)*Utility::fgkBinWidth,integral);
    } // end for int bin < length
  }
  
  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMSumWaveforms::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Coonfiguration file");

  c("TriggerType").Get(fTriggerType);
  MsgInfo(MsgLog::Form("\t-Trigger Type %s",fTriggerType.c_str()));

  c("DoBCMCut").Get(fDoBCMCut);
  MsgInfo(MsgLog::Form("\t-Apply BCM Cut %s",fDoBCMCut? "true" : "false"));
  if (fDoBCMCut) {
    c("BCMTimeLowCut").Get(fBCMTimeLowCut);
    c("BCMTimeHighCut").Get(fBCMTimeHighCut);
    c("BCMLengthLowCut").Get(fBCMLengthLowCut);
    c("BCMLengthHighCut").Get(fBCMLengthHighCut);
    c("BCMIntegralLowCut").Get(fBCMIntegralLowCut);
    c("BCMIntegralHighCut").Get(fBCMIntegralHighCut);
    MsgInfo(MsgLog::Form("\t\t-BCMTime Cut Low %g high %g",fBCMTimeLowCut,fBCMTimeHighCut));
    MsgInfo(MsgLog::Form("\t\t-BCMLength Cut Low %g high %g",fBCMLengthLowCut,fBCMLengthHighCut));
    MsgInfo(MsgLog::Form("\t\t-BCMIntegral Cut Low %g high %g",fBCMIntegralLowCut,fBCMIntegralHighCut));
  }

  c("ApplyFP3Offset").Get(fApplyFP3Offset);
  MsgInfo(MsgLog::Form("\t-Apply FP3Offset %s",fApplyFP3Offset? "true" : "false"));

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
    for (auto & hist : fTimeHistsHit) {
      if (hist != nullptr) {
        hist->Write();
      }
    }
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


