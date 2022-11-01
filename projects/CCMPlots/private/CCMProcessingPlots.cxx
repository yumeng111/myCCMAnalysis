/*!**********************************************
 * \file CCMProcessingPlots.cxx
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 * \brief Read in events and apply some cuts
 *
 * Read in previously found events and apply
 * some cuts.
 *
 ***********************************************/

#include <cmath>
#include <limits>
#include <memory>
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>

#include "CCMAnalysis/CCMPlots/CCMProcessingPlots.h"

#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMModuleTable.h"
#include "CCMAnalysis/CCMDataStructures/Events.h"
#include "CCMAnalysis/CCMDataStructures/SimplifiedEvent.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/PMTInfoMap.h"
#include "CCMAnalysis/CCMUtils/CCMBeamInfo.h"
#include "CCMAnalysis/CCMUtils/PMTInformation.h"

#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TTree.h"

//See CCMModuleTable for info
MODULE_DECL(CCMProcessingPlots);

//_______________________________________________________________________________________
CCMProcessingPlots::CCMProcessingPlots(const char* version) 
  : CCMModule("CCMProcessingPlots"), fEvents(nullptr), fOutFileName(""), fInFileName(""), fPrevInFileName(""),
    fOutfile(nullptr), fInfile(nullptr), fTree(nullptr), fCCMPulsesTimeHist(nullptr), fCCMPulsesTimeIntHist(nullptr),
    fCCMVetoTimeHist(nullptr), fCCMVetoTimeIntHist(nullptr), fFP3TimeHist(nullptr), fFP3TimeIntHist(nullptr),
    fCCMEventsTimeHist(nullptr), fCCMEventsTimeIntHist(nullptr), fBCM3DHist(nullptr),
    fPromptTopVeto(nullptr), fPromptBottomVeto(nullptr), fPromptCFrontVeto(nullptr),
    fPromptCBackVeto(nullptr), fPromptCLeftVeto(nullptr), fPromptCRightVeto(nullptr), fPromptTotalVeto(nullptr),
    fPrevCurrentInfo(), fNextCurrentInfo(), fNumberOfTriggersInSum(0.0),fBCMIntegralTime(),fBCMTimeTime(),
    fBCMWidthTime(),fPromptTopVetoTime(),fPromptBottomVetoTime(),fPromptCLeftVetoTime(),fPromptCRightVetoTime(),
    fPromptCFrontVetoTime(),fPromptCBackVetoTime(),fPromptTotalVetoTime()
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMProcessingPlots::CCMProcessingPlots(const CCMProcessingPlots& clufdr) 
: CCMModule(clufdr), fEvents(clufdr.fEvents), fOutFileName(clufdr.fOutFileName),
  fInFileName(clufdr.fInFileName), fPrevInFileName(clufdr.fPrevInFileName),
  fOutfile(clufdr.fOutfile), fInfile(clufdr.fInfile),
  fCCMPulsesTimeHist(clufdr.fCCMPulsesTimeHist), fCCMPulsesTimeIntHist(clufdr.fCCMPulsesTimeIntHist),
  fCCMVetoTimeHist(clufdr.fCCMVetoTimeHist), fCCMVetoTimeIntHist(clufdr.fCCMVetoTimeIntHist),
  fFP3TimeHist(clufdr.fFP3TimeHist), fFP3TimeIntHist(clufdr.fFP3TimeIntHist),
  fCCMEventsTimeHist(clufdr.fCCMEventsTimeHist), fCCMEventsTimeIntHist(clufdr.fCCMEventsTimeIntHist),
  fBCM3DHist(clufdr.fBCM3DHist), fPromptTopVeto(clufdr.fPromptTopVeto),
  fPromptBottomVeto(clufdr.fPromptBottomVeto), fPromptCFrontVeto(clufdr.fPromptCFrontVeto),
  fPromptCBackVeto(clufdr.fPromptCBackVeto), fPromptCLeftVeto(clufdr.fPromptCLeftVeto),
  fPromptCRightVeto(clufdr.fPromptCRightVeto), fPromptTotalVeto(clufdr.fPromptTotalVeto),
  fPrevCurrentInfo(clufdr.fPrevCurrentInfo), fNextCurrentInfo(clufdr.fNextCurrentInfo),
  fNumberOfTriggersInSum(clufdr.fNumberOfTriggersInSum),
  fBCMIntegralTime(clufdr.fBCMIntegralTime), fBCMTimeTime(clufdr.fBCMTimeTime),
  fBCMWidthTime(clufdr.fBCMWidthTime), fPromptTopVetoTime(clufdr.fPromptTopVetoTime),
  fPromptBottomVetoTime(clufdr.fPromptBottomVetoTime), fPromptCLeftVetoTime(clufdr.fPromptCLeftVetoTime),
  fPromptCRightVetoTime(clufdr.fPromptCRightVetoTime), fPromptCFrontVetoTime(clufdr.fPromptCFrontVetoTime), 
  fPromptCBackVetoTime(clufdr.fPromptCBackVetoTime), fPromptTotalVetoTime(clufdr.fPromptTotalVetoTime)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMProcessingPlots::~CCMProcessingPlots()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMProcessingPlots::ProcessTrigger()
{
  auto beamTime = fEvents->GetBeamTime()*Utility::fgkBinWidth + Utility::fgkWindowStartTime;
  auto beamLength = fEvents->GetBeamLength();
  auto beamIntegral = fEvents->GetBeamIntegral();

  auto wallTime = fEvents->GetComputerSecIntoEpoch();

  auto wallTimeC = CCMBeamInfo::ConvertTimeSinceEPOCHtoTime(wallTime);

  if (wallTimeC > fNextCurrentInfo.first) {
    if (fNumberOfTriggersInSum >= 1.0) {
      fOutfile->cd();
      fTree->Fill();
    }

    CCMBeamInfo::Find(wallTimeC);
    fPrevCurrentInfo = CCMBeamInfo::GetBeamInfo();
    CCMBeamInfo::Next();
    fNextCurrentInfo = CCMBeamInfo::GetBeamInfo();

    fCurrentTime = static_cast<unsigned int>(fNextCurrentInfo.first);
    fBCMIntegralTime.clear();
    fBCMTimeTime.clear();
    fBCMWidthTime.clear();
    fPromptTopVetoTime.clear();
    fPromptBottomVetoTime.clear();
    fPromptCLeftVetoTime.clear();
    fPromptCRightVetoTime.clear();
    fPromptCFrontVetoTime.clear();
    fPromptCBackVetoTime.clear();
    fPromptTopVetoTime.clear();
    fNumberOfTriggersInSum = 0.0;

  } // end if in new time region

  fBCMIntegralTime.emplace_back(beamIntegral);
  fBCMWidthTime.emplace_back(beamLength);
  fBCMTimeTime.emplace_back(beamTime);

  ++fNumberOfTriggersInSum;

  double bcmVariables[3] = {beamTime+1e-4,beamIntegral,beamLength+1e-4};

  fBCM3DHist->Fill(bcmVariables);

  size_t numEvents = fEvents->GetNumEvents();
  for (size_t event = 0; event < numEvents; ++event) {
    auto simplifiedEvent = fEvents->GetSimplifiedEvent(event);
    if (simplifiedEvent.GetEventFinderMethod() != kCCMDynamicLengthEventID ||
        simplifiedEvent.GetAccumWaveformMethod() != kCCMAccumWaveformTriangleID) {
      continue;
    }
    auto st = simplifiedEvent.GetStartTime();

    fCCMEventsTimeHist->Fill(st+1e-4);
    fCCMEventsTimeIntHist->Fill(st+1e-4,simplifiedEvent.GetIntegralTank());

    if (st > -600) {
      break;
    }

    auto topVeto = simplifiedEvent.GetNumVetoTop(true);
    auto bottomVeto = simplifiedEvent.GetNumVetoBottom(true);
    auto leftVeto = simplifiedEvent.GetNumVetoLeft(true);
    auto rightVeto = simplifiedEvent.GetNumVetoRight(true);
    auto frontVeto = simplifiedEvent.GetNumVetoFront(true);
    auto backVeto = simplifiedEvent.GetNumVetoBack(true);
    auto totalVeto = simplifiedEvent.GetNumVeto(true);

    fPromptTopVeto->Fill(topVeto);
    fPromptBottomVeto->Fill(bottomVeto);
    fPromptCLeftVeto->Fill(leftVeto);
    fPromptCRightVeto->Fill(rightVeto);
    fPromptCFrontVeto->Fill(frontVeto);
    fPromptCBackVeto->Fill(backVeto);
    fPromptTotalVeto->Fill(totalVeto);

    fPromptTopVetoTime.emplace_back(topVeto);
    fPromptBottomVetoTime.emplace_back(bottomVeto);
    fPromptCLeftVetoTime.emplace_back(leftVeto);
    fPromptCRightVetoTime.emplace_back(rightVeto);
    fPromptCFrontVetoTime.emplace_back(frontVeto);
    fPromptCBackVetoTime.emplace_back(backVeto);
    fPromptTotalVetoTime.emplace_back(totalVeto);


    //auto wfInt = simplifiedEvent.GetWaveformInt();
    //for (size_t loc = 0; loc < wfInt.size(); ++loc) {
    //  fCCMEventsTimeIntHist->Fill(st+loc*Utility::fgkBinWidth+1e-4,wfInt.at(loc));
    //}
  } // end for event

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMProcessingPlots::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Coonfiguration file");

  std::string temp = "";
  c("BeamInfoFile").Get(temp);
  CCMBeamInfo::LoadTable(temp);

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMProcessingPlots::SetupOutFile()
 * \brief Setup the output file that gets saved
 ***********************************************/
void CCMProcessingPlots::SetupOutFile()
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
    fTree = new TTree("timeVariables","timeVariables");
    fTree->Branch("time",&fCurrentTime);
    fTree->Branch("numTriggers",&fNumberOfTriggersInSum);
    fTree->Branch("bcmTime",&fBCMTimeTime);
    fTree->Branch("bcmLength",&fBCMWidthTime);
    fTree->Branch("bcmIntegral",&fBCMIntegralTime);
    fTree->Branch("vetoTop",&fPromptTopVetoTime);
    fTree->Branch("vetoBottom",&fPromptBottomVetoTime);
    fTree->Branch("vetoFront",&fPromptCFrontVetoTime);
    fTree->Branch("vetoBack",&fPromptCBackVetoTime);
    fTree->Branch("vetoLeft",&fPromptCLeftVetoTime);
    fTree->Branch("vetoRight",&fPromptCRightVetoTime);
    fTree->Branch("vetoTotal",&fPromptTotalVetoTime);
  }
}

/*!**********************************************
 * \fn void CCMProcessingPlots::SetupInputFile()
 * \brief Setup the input file that was read in to get the hitograms
 ***********************************************/
void CCMProcessingPlots::SetupInputFile()
{
  if (fInFileName== fPrevInFileName) {
    return;
  }

  if (!fInFileName.empty()) {
    fInfile= gROOT->GetFile(fInFileName.c_str());
    if (!fInfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fInFileName.c_str()));
      return;
    }
    fInfile->cd();

    TH1D * tempHist = 0;

    fInfile->GetObject("timeHist_CCM",tempHist);
    if (tempHist) {
      if (fCCMPulsesTimeHist == nullptr) {
        fCCMPulsesTimeHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(tempHist->Clone("timeHist_CCM_day")));
      } else {
        fCCMPulsesTimeHist->Add(tempHist);
      }
      delete tempHist;
    }

    fInfile->GetObject("timeHistInt_CCM",tempHist);
    if (tempHist) {
      if (fCCMPulsesTimeIntHist == nullptr) {
        fCCMPulsesTimeIntHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(tempHist->Clone("timeHistInt_CCM_day")));
      } else {
        fCCMPulsesTimeIntHist->Add(tempHist);
      }
      delete tempHist;
    }

    fInfile->GetObject("timeHist_CCMVeto",tempHist);
    if (tempHist) {
      if (fCCMVetoTimeHist == nullptr) {
        fCCMVetoTimeHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(tempHist->Clone("timeHist_CCMVeto_day")));
      } else {
        fCCMVetoTimeHist->Add(tempHist);
      }
      delete tempHist;
    }

    fInfile->GetObject("timeHistInt_CCMVeto",tempHist);
    if (tempHist) {
      if (fCCMVetoTimeIntHist == nullptr) {
        fCCMVetoTimeIntHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(tempHist->Clone("timeHistInt_CCMVeto_day")));
      } else {
        fCCMVetoTimeIntHist->Add(tempHist);
      }
      delete tempHist;
    }

    fInfile->GetObject("timeHist_EJ3015C_FP3",tempHist);
    if (tempHist) {
      if (fFP3TimeHist == nullptr) {
        fFP3TimeHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(tempHist->Clone("timeHist_FP3_day")));
      } else {
        fFP3TimeHist->Add(tempHist);
      }
      delete tempHist;
    }

    fInfile->GetObject("timeHistInt_EJ3015C_FP3",tempHist);
    if (tempHist) {
      if (fFP3TimeIntHist == nullptr) {
        fFP3TimeIntHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(tempHist->Clone("timeHistInt_FP3_day")));
      } else {
        fFP3TimeIntHist->Add(tempHist);
      }
      delete tempHist;
    }
  }

  if (fCCMEventsTimeHist == nullptr) {
    fCCMEventsTimeHist = std::make_shared<TH1D>("timeHist_CCMEvents_day",
          "CCM;Time (ns);Count",Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime);
    fCCMEventsTimeIntHist = std::make_shared<TH1D>("timeHistInt_CCMEvents_day",
          "CCM;Time (ns);Count",Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime);
    fPromptTopVeto = std::make_shared<TH1D>("promptTopVeto","Prompt Top Veto;Number in prompt;Count",100,0,100);
    fPromptBottomVeto = std::make_shared<TH1D>("promptBottomVeto","Prompt Bottom Veto;Number in prompt;Count",100,0,100);
    fPromptCLeftVeto = std::make_shared<TH1D>("promptCLeftVeto","Prompt CLeft Veto;Number in prompt;Count",100,0,100);
    fPromptCRightVeto = std::make_shared<TH1D>("promptCRightVeto","Prompt CRight Veto;Number in prompt;Count",100,0,100);
    fPromptCFrontVeto = std::make_shared<TH1D>("promptCFrontVeto","Prompt CFront Veto;Number in prompt;Count",100,0,100);
    fPromptCBackVeto = std::make_shared<TH1D>("promptCBackVeto","Prompt CBack Veto;Number in prompt;Count",100,0,100);
    fPromptTotalVeto = std::make_shared<TH1D>("promptTotalVeto","Prompt Total Veto;Number in prompt;Count",1000,0,1000);

    int nBins[3] = {Utility::fgkNumBins,static_cast<int>(std::pow(2,16)),Utility::fgkNumBins};
    double xmin[3] = {Utility::fgkWindowStartTime,0,0};
    double xmax[3] = {Utility::fgkWindowEndTime,std::pow(2,16),Utility::fgkNumBins*Utility::fgkBinWidth};

    fBCM3DHist = std::make_shared<THnSparseF>("bcm3DHist",";Time (ns);Integral (adc);Time (ns);Count",3,nBins,xmin,xmax);
  }
}

//------------------------------------------------------
CCMResult_t CCMProcessingPlots::EndOfJob()
{
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  if (fOutfile != nullptr) {
    fOutfile->cd();
    if (fTree) {
      fTree->Fill();
      fTree->Write();
      MsgInfo(MsgLog::Form("Number of Trigger times = %ld",fTree->GetEntries()));
    }

    if (fCCMPulsesTimeHist) {
      fCCMPulsesTimeHist->Write();
    }
    if (fCCMPulsesTimeIntHist) {
      fCCMPulsesTimeIntHist->Write();
    }
    if (fCCMVetoTimeHist) {
      fCCMVetoTimeHist->Write();
    }
    if (fCCMVetoTimeIntHist) {
      fCCMVetoTimeIntHist->Write();
    }
    if (fFP3TimeHist) {
      fFP3TimeHist->Write();
    }
    if (fFP3TimeIntHist) {
      fFP3TimeIntHist->Write();
    }
    if (fCCMEventsTimeHist) {
      fCCMEventsTimeHist->Write();
    }
    if (fCCMEventsTimeIntHist) {
      fCCMEventsTimeIntHist->Write();
    }
    if (fBCM3DHist) {
      fBCM3DHist->Write();
    }
    if (fPromptTopVeto) {
      fPromptTopVeto->Write();
    }
    if (fPromptBottomVeto) {
      fPromptBottomVeto->Write();
    }
    if (fPromptCFrontVeto) {
      fPromptCFrontVeto->Write();
    }
    if (fPromptCBackVeto) {
      fPromptCBackVeto->Write();
    }
    if (fPromptCLeftVeto) {
      fPromptCLeftVeto->Write();
    }
    if (fPromptCRightVeto) {
      fPromptCRightVeto->Write();
    }
    if (fPromptTotalVeto) {
      fPromptTotalVeto->Write();
    }
  }

  return kCCMSuccess;
}


