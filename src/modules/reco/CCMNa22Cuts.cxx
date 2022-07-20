/*!**********************************************
 * \file CCMNa22Cuts.cxx
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

#include "CCMAnalysis/modules/reco/CCMNa22Cuts.h"

#include "CCMAnalysis/ds/Events.h"
#include "CCMAnalysis/ds/SimplifiedEvent.h"
#include "CCMAnalysis/modules/framework/CCMConfig.h"
#include "CCMAnalysis/modules/framework/CCMConfigParam.h"
#include "CCMAnalysis/modules/framework/CCMModuleTable.h"
#include "CCMAnalysis/utils/MsgLog.h"
#include "CCMAnalysis/utils/PMTInfoMap.h"
#include "CCMAnalysis/utils/PMTInformation.h"

#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile.h"
#include "TVector3.h"

//See CCMModuleTable for info
MODULE_DECL(CCMNa22Cuts);

//_______________________________________________________________________________________
CCMNa22Cuts::CCMNa22Cuts(const char* version) 
  : CCMModule("CCMNa22Cuts")
{
  //Default constructor
  this->SetCfgVersion(version);

  fEnergyLengthFile = nullptr;
  fEnergyLengthProf = nullptr;

  fEvents = nullptr;
  fEventFinderID = kCCMDynamicLengthEventID;

  fAccumWaveformMethodID = kCCMAccumWaveformTriangleID;
  fOutFileName = "";
  fTreeName = "tree";

  fRemovePrimaryEvents = true;
  fRemoveOtherEvents = false;
  fDoBCMCut = false;
  fDoVetoCut = false;
  fDoPositionCut = false;

  fDoEnergyCut = false;
  fDoPrevCut = false;
  fDoLengthCut = false;
  fDoTimeCut = false;
  fDoNicenessCut = false;
  fReverseVeto = false;
  fWaveformMaxCut = false;

  fPassedVetoCut = true;
  fPassedPositionCut = true;
  fPassedEnergyCut = true;
  fPassedPrevCut = true;
  fPassedLengthCut = true;
  fPassedTimeCut = true;
  fPassedNicenessCut = true;
  fPassedWaveformMaxCut = true;
  fPassedNumHitsCut = true;
  fPassedBCMCut = false;

  fNumFailedVetoCut = 0;
  fNumFailedPositionCut = 0;
  fNumFailedEnergyCut = 0;
  fNumFailedPrevCut = 0;
  fNumFailedLengthCut = 0;
  fNumFailedTimeCut = 0;
  fNumFailedNicenessCut = 0;
  fNumFailedWaveformMaxCut = 0;
  fNumFailedNumHitsCut = 0;
  fNumFailedBCMCut = 0;

  fNumVetoCut = 3;
  fBCMTimeLowCut = Utility::fgkWindowStartTime;
  fBCMTimeHighCut = Utility::fgkWindowEndTime;

  fBCMLengthLowCut = 0;
  fBCMLengthHighCut = Utility::fgkNumBins*Utility::fgkBinWidth;

  fBCMIntegralLowCut = 0;
  fBCMIntegralHighCut = std::numeric_limits<double>::max();

  fRadiusCutValueLow = 0;
  fRadiusCutValueHigh = 0.95;
  fZCutValueLow = -0.4;
  fZCutValueHigh = 0.4;

  fEnergyCutValueLow = 0;
  fEnergyCutValueHigh = std::numeric_limits<double>::max();
  fPrevCutTime = 1600;
  fPrevCutHoldOff = 90.0;

  fPrevCutEnergyFrac = 0;
  fLengthCutValueLow = -std::numeric_limits<double>::max();

  fLengthCutValueHigh = std::numeric_limits<double>::max();
  fTimeCutValueLow = -7000;
  fTimeCutValueHigh = 4000;

  fNicenessCutValueLow = 0;
  fNicenessCutValueHigh = 2.5;
  fOutfile = nullptr;
  fTree = nullptr;
  fEnergy = 0;
  fTotalEnergy = 0;
  fLength = 0;

  fHits = 0;
  fNumPMTs = 0;
  fTime = 0;
  fVetoTop = 0;
  fVetoBottom = 0;
  fVetoCLeft = 0;
  fVetoCRight = 0;
  fVetoCFront = 0;

  fVetoCBack = 0;
  fVetoPromptTop = 0;
  fVetoPromptBottom = 0;
  fVetoPromptCLeft = 0;
  fVetoPromptCRight = 0;

  fVetoPromptCFront = 0;
  fVetoPromptCBack = 0;
  fWeight = 0;
  fX = 0;
  fY = 0;
  fZ = 0;
  fLargestPMTFraction = 0;

  fEpochSec = 0;
  fEpochNSSec = 0;
  fBCMTime = 0;
  fBCMLength = 0;
  fBCMIntegral = 0;
  fNumInitEvents = 0;
  fNumFinalEvents = 0;
}

//_______________________________________________________________________________________
CCMNa22Cuts::CCMNa22Cuts(const CCMNa22Cuts& clufdr) 
  : CCMModule(clufdr)
{
  // copy constructor
  fEvents = clufdr.fEvents;
  fEventFinderID = clufdr.fEventFinderID;
  fAccumWaveformMethodID = clufdr.fAccumWaveformMethodID;
  fOutFileName = clufdr.fOutFileName;
  fTreeName = clufdr.fTreeName;
  fRemovePrimaryEvents = clufdr.fRemovePrimaryEvents;
  fRemoveOtherEvents = clufdr.fRemoveOtherEvents;
  fDoBCMCut = clufdr.fDoBCMCut;
  fDoVetoCut = clufdr.fDoVetoCut;
  fDoPositionCut = clufdr.fDoPositionCut;
  fDoEnergyCut = clufdr.fDoEnergyCut;
  fDoPrevCut = clufdr.fDoPrevCut;
  fDoLengthCut = clufdr.fDoLengthCut;
  fDoTimeCut = clufdr.fDoTimeCut;
  fDoNicenessCut = clufdr.fDoNicenessCut;
  fReverseVeto = clufdr.fReverseVeto;
  fWaveformMaxCut = clufdr.fWaveformMaxCut;
  fPassedVetoCut = clufdr.fPassedVetoCut;
  fPassedPositionCut = clufdr.fPassedPositionCut;
  fPassedEnergyCut = clufdr.fPassedEnergyCut;
  fPassedPrevCut = clufdr.fPassedPrevCut;
  fPassedLengthCut = clufdr.fPassedLengthCut;
  fPassedTimeCut = clufdr.fPassedTimeCut;
  fPassedNicenessCut = clufdr.fPassedNicenessCut;
  fPassedWaveformMaxCut = clufdr.fPassedWaveformMaxCut;
  fPassedNumHitsCut = clufdr.fPassedNumHitsCut;
  fPassedBCMCut = clufdr.fPassedBCMCut;
  fNumFailedVetoCut = clufdr.fPassedVetoCut;
  fNumFailedPositionCut = clufdr.fPassedPositionCut;
  fNumFailedEnergyCut = clufdr.fPassedEnergyCut;
  fNumFailedPrevCut = clufdr.fPassedPrevCut;
  fNumFailedLengthCut = clufdr.fPassedLengthCut;
  fNumFailedTimeCut = clufdr.fPassedTimeCut;
  fNumFailedNicenessCut = clufdr.fPassedNicenessCut;
  fNumFailedWaveformMaxCut = clufdr.fPassedWaveformMaxCut;
  fNumFailedNumHitsCut = clufdr.fPassedNumHitsCut;
  fNumFailedBCMCut = clufdr.fPassedBCMCut;
  fNumVetoCut = clufdr.fNumVetoCut;
  fBCMTimeLowCut = clufdr.fBCMTimeLowCut;
  fBCMTimeHighCut = clufdr.fBCMTimeHighCut;
  fBCMLengthLowCut = clufdr.fBCMLengthLowCut;
  fBCMLengthHighCut = clufdr.fBCMLengthHighCut;
  fBCMIntegralLowCut = clufdr.fBCMIntegralLowCut;
  fBCMIntegralHighCut = clufdr.fBCMIntegralHighCut;
  fRadiusCutValueLow = clufdr.fRadiusCutValueLow;
  fRadiusCutValueHigh = clufdr.fRadiusCutValueHigh;
  fZCutValueLow = clufdr.fZCutValueLow;
  fZCutValueHigh = clufdr.fZCutValueHigh;
  fEnergyCutValueLow = clufdr.fEnergyCutValueLow;
  fEnergyCutValueHigh = clufdr.fEnergyCutValueHigh;
  fPrevCutTime = clufdr.fPrevCutTime;
  fPrevCutHoldOff = clufdr.fPrevCutHoldOff;
  fPrevCutEnergyFrac = clufdr.fPrevCutEnergyFrac;
  fLengthCutValueLow = clufdr.fLengthCutValueLow;
  fLengthCutValueHigh = clufdr.fLengthCutValueHigh;
  fTimeCutValueLow = clufdr.fTimeCutValueLow;
  fTimeCutValueHigh = clufdr.fTimeCutValueHigh;
  fNicenessCutValueLow = clufdr.fNicenessCutValueLow;
  fNicenessCutValueHigh = clufdr.fNicenessCutValueHigh;
  fOutfile = clufdr.fOutfile;
  fTree = clufdr.fTree;
  fEnergy = clufdr.fEnergy;
  fTotalEnergy = clufdr.fTotalEnergy;
  fLength = clufdr.fLength;
  fHits = clufdr.fHits;
  fNumPMTs = clufdr.fNumPMTs;
  fTime = clufdr.fTime;
  fVetoTop = clufdr.fVetoTop;
  fVetoBottom = clufdr.fVetoBottom;
  fVetoCLeft = clufdr.fVetoCLeft;
  fVetoCRight = clufdr.fVetoCRight;
  fVetoCFront = clufdr.fVetoCFront;
  fVetoCBack = clufdr.fVetoCBack;
  fVetoPromptTop = clufdr.fVetoPromptTop;
  fVetoPromptBottom = clufdr.fVetoPromptBottom;
  fVetoPromptCLeft = clufdr.fVetoPromptCLeft;
  fVetoPromptCRight = clufdr.fVetoPromptCRight;
  fVetoPromptCFront = clufdr.fVetoPromptCFront;
  fVetoPromptCBack = clufdr.fVetoPromptCBack;
  fWeight = clufdr.fWeight;
  fX = clufdr.fX;
  fY = clufdr.fY;
  fZ = clufdr.fZ;
  fLargestPMTFraction = clufdr.fLargestPMTFraction;
  fEpochSec = clufdr.fEpochSec;
  fEpochNSSec = clufdr.fEpochNSSec;
  fBCMTime = clufdr.fBCMTime;
  fBCMLength = clufdr.fBCMLength;
  fBCMIntegral = clufdr.fBCMIntegral;
  fNumInitEvents = clufdr.fNumInitEvents;
  fNumFinalEvents = clufdr.fNumFinalEvents;
}

//_______________________________________________________________________________________
CCMNa22Cuts::~CCMNa22Cuts()
{ 
  // destructor
  if (fEnergyLengthProf) {
    delete fEnergyLengthProf;
  }
}

//_______________________________________________________________________________________
CCMResult_t CCMNa22Cuts::ProcessTrigger()
{
  // since the class does not own the TFile, always check to see if it is present
  fOutfile = gROOT->GetFile(fOutFileName.c_str());

  std::vector<size_t> locationsToRemove;

  unsigned long int tempCounter = 0;

  const size_t kNumEvents = fEvents->GetNumEvents();
  if (kNumEvents == 0) {
    return kCCMDoNotWrite;
  }

  fEpochSec = fEvents->GetComputerSecIntoEpoch();
  fEpochNSSec = fEvents->GetComputerNSIntoSec();

  fBCMTime = fEvents->GetBeamTime()*Utility::fgkBinWidth + Utility::fgkWindowStartTime;
  fBCMLength = fEvents->GetBeamLength();
  fBCMIntegral = fEvents->GetBeamIntegral();

  fPassedBCMCut = PassedBCMCut(fBCMTime,fBCMLength,fBCMIntegral);
  if (fDoBCMCut && !fPassedBCMCut) {
    return kCCMDoNotWrite;
  }

  for (size_t e = 0; e < kNumEvents; ++e) {
    fPassedVetoCut = true;
    fPassedPositionCut = true;
    fPassedEnergyCut = true;
    fPassedPrevCut = true;
    fPassedLengthCut = true;
    fPassedTimeCut = true;
    fPassedNicenessCut = true;
    fPassedWaveformMaxCut = true;
    fPassedNumHitsCut = true;

    auto simplifiedEvent = fEvents->GetSimplifiedEvent(e);

    if (simplifiedEvent.GetEventFinderMethod() != fEventFinderID ||
        simplifiedEvent.GetAccumWaveformMethod() != fAccumWaveformMethodID) {
      if (fRemoveOtherEvents) {
        locationsToRemove.emplace_back(e);
      }
      continue;
    }

    ++fNumInitEvents;
    ++tempCounter;

    fTime = simplifiedEvent.GetStartTime();
    fLength = simplifiedEvent.GetLength();

    fEnergy = simplifiedEvent.GetIntegralTank(true);
    fTotalEnergy = simplifiedEvent.GetIntegralTank(false);
    fHits = simplifiedEvent.GetNumTank(true);
    fNumPMTs = simplifiedEvent.GetPMTHits();

    fX = simplifiedEvent.GetXPosition();
    fY = simplifiedEvent.GetYPosition();
    fZ = simplifiedEvent.GetZPosition();

    fVetoTop = simplifiedEvent.GetNumVetoTop(false);
    fVetoBottom = simplifiedEvent.GetNumVetoBottom(false);
    fVetoCLeft = simplifiedEvent.GetNumVetoLeft(false);
    fVetoCBack = simplifiedEvent.GetNumVetoBack(false);
    fVetoCFront = simplifiedEvent.GetNumVetoFront(false);
    fVetoCRight = simplifiedEvent.GetNumVetoRight(false);

    fVetoPromptTop = simplifiedEvent.GetNumVetoTop(true);
    fVetoPromptBottom = simplifiedEvent.GetNumVetoBottom(true);
    fVetoPromptCLeft = simplifiedEvent.GetNumVetoLeft(true);
    fVetoPromptCBack = simplifiedEvent.GetNumVetoBack(true);
    fVetoPromptCFront = simplifiedEvent.GetNumVetoFront(true);
    fVetoPromptCRight = simplifiedEvent.GetNumVetoRight(true);

    fLargestPMTFraction = simplifiedEvent.GetLargestPMTFraction();


    if (fWaveformMaxCut) {
      fPassedWaveformMaxCut = PassedWaveformMaxCut(fLength,simplifiedEvent.GetMaxAccumWaveformTime());
      if (!fPassedWaveformMaxCut && std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    }

    if (fDoTimeCut) {
      fPassedTimeCut = PassedTimeCut(fTime);
      if (!fPassedTimeCut && std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    } // end fDoTimeCut

    if (fTime + fLength < fTime) {
      if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    }

    if (fDoLengthCut) {
      fPassedLengthCut = PassedLengthCut(fLength);
      if (!fPassedLengthCut && std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    } // end if fDoLengthCut

    if (fDoPositionCut) {
      fPassedPositionCut = PassedPositionCut(fX,fY,fZ);
      if (!fPassedPositionCut && 
          std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    } // end if fDoPositionCut

    if (fDoEnergyCut) {
      fPassedEnergyCut = PassedEnergyCut(fEnergy);
      if (!fPassedEnergyCut && std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    }// end if fDoEnergyCut

    fPassedNumHitsCut = PassedNumHitsCut(fHits);
    if (!fPassedNumHitsCut) {
      if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    }

    if (fDoNicenessCut) {
      fPassedNicenessCut = PassedNicenessCut(fEnergy,fNumPMTs,fLargestPMTFraction);
      if (!fPassedNicenessCut && std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    } // end if fDoNicenessCut

    if (fDoVetoCut) {
      fPassedVetoCut = PassedVetoCut(simplifiedEvent.GetNumVeto(true));
      if (!fPassedVetoCut && 
          std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    }// end if fDoVetoCut


    if (fDoPrevCut) {
      fPassedPrevCut = PassedPrevCut(fTime,e,fEnergy);
      if (!fPassedPrevCut && std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && 
          fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
    }

    fWeight = 1.0;// /weightFunction.at(energyRegion)->Eval(fTime);
    //MsgInfo(MsgLog::Form("fTime %f weighFunction %g weight %g",fTime,weightFunction->Eval(fTime),weight));

    if (fOutfile != nullptr) {
      //MsgInfo("Filling the tree");
      fOutfile->cd();
      fTree->Fill();
    }

    if (!fPassedBCMCut) {
      ++fNumFailedBCMCut;
    } else if (!fPassedTimeCut) {
      ++fNumFailedTimeCut;
    } else if (!fPassedNumHitsCut) {
      ++fNumFailedNumHitsCut;
    } else if (!fPassedPrevCut) {
      ++fNumFailedPrevCut;
    } else if (!fPassedVetoCut) {
      ++fNumFailedVetoCut;
    } else if (!fPassedPositionCut) {
      ++fNumFailedPositionCut;
    } else if (!fPassedNicenessCut) {
      ++fNumFailedNicenessCut;
    } else if (!fPassedWaveformMaxCut) {
      ++fNumFailedWaveformMaxCut;
    } else if (!fPassedLengthCut) {
      ++fNumFailedLengthCut;
    } else if (!fPassedEnergyCut) {
      ++fNumFailedEnergyCut;
    }
  } // end e

  fNumFinalEvents += tempCounter - locationsToRemove.size();

  if (locationsToRemove.empty()) {
    return kCCMSuccess;
  }

  for (auto it = locationsToRemove.rbegin();
      it != locationsToRemove.rend();
      ++it) {
    fEvents->RemoveSimplifiedEvent(*it);
  }

  if (fEvents->GetNumEvents() == 0) {
    return kCCMDoNotWrite;
  }

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMNa22Cuts::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Coonfiguration file");

  std::string tempString = "";
  c("EventFinderID").Get(tempString);
  MsgInfo(MsgLog::Form("\t-EventFinderID %s",tempString.c_str()));
  fEventFinderID = Utility::ConvertStringToCCMEventFinderID(tempString);

  c("AccumWaveformMethodID").Get(tempString);
  MsgInfo(MsgLog::Form("\t-AccumWaveformMethodID %s",tempString.c_str()));
  fAccumWaveformMethodID = Utility::ConvertStringToCCMAccumWaveformMethod(tempString);

  c("RemovePrimaryEvents").Get(fRemovePrimaryEvents);
  MsgInfo(MsgLog::Form("\t-Remove Primary Events %s",fRemovePrimaryEvents ? "true" : "false"));
  c("RemoveOtherEvents").Get(fRemoveOtherEvents);
  MsgInfo(MsgLog::Form("\t-Remove Other Events %s",fRemoveOtherEvents ? "true" : "false"));
  c("TreeName").Get(fTreeName);
  MsgInfo(MsgLog::Form("\t-Output Tree Name %s",fTreeName.c_str()));

  c("DoLengthCut").Get(fDoLengthCut);
  MsgInfo(MsgLog::Form("\t-Apply Length Cut %s",fDoLengthCut? "true" : "false"));
  if (fDoLengthCut) {
    c("LengthCutValueLow").Get(fLengthCutValueLow);
    c("LengthCutValueHigh").Get(fLengthCutValueHigh);
    MsgInfo(MsgLog::Form("\t\t-Length Cut Low %g high %g",fLengthCutValueLow,fLengthCutValueHigh));
  }

  c("DoTimeCut").Get(fDoTimeCut);
  MsgInfo(MsgLog::Form("\t-Apply Time Cut %s",fDoTimeCut? "true" : "false"));
  if (fDoTimeCut) {
    c("TimeCutValueLow").Get(fTimeCutValueLow);
    c("TimeCutValueHigh").Get(fTimeCutValueHigh);
    MsgInfo(MsgLog::Form("\t\t-Time Cut Low %g high %g",fTimeCutValueLow,fTimeCutValueHigh));
  }

  c("DoNicenessCut").Get(fDoNicenessCut);
  MsgInfo(MsgLog::Form("\t-Apply Niceness Cut %s",fDoNicenessCut? "true" : "false"));
  if (fDoNicenessCut) {
    c("NicenessCutValueLow").Get(fNicenessCutValueLow);
    c("NicenessCutValueHigh").Get(fNicenessCutValueHigh);
    MsgInfo(MsgLog::Form("\t\t-Niceness Cut Low %g high %g",fNicenessCutValueLow,fNicenessCutValueHigh));
  }

  c("DoPrevCut").Get(fDoPrevCut);
  MsgInfo(MsgLog::Form("\t-Apply Prev Cut %s",fDoPrevCut? "true" : "false"));
  if (fDoPrevCut) {
    c("PrevCutTime").Get(fPrevCutTime);
    c("PrevCutHoldOff").Get(fPrevCutHoldOff);
    //c("PrevCutEnergyFrac").Get(fPrevCutEnergyFrac);
    MsgInfo(MsgLog::Form("\t\t-Prev Cut Time %g Required Energy Fraction %g Hold off (ns) %g",
          fPrevCutTime,fPrevCutEnergyFrac,fPrevCutHoldOff));
  }

  c("WaveformMaxCut").Get(fWaveformMaxCut);
  MsgInfo(MsgLog::Form("\t-Apply Waveform Max Cut %s",fWaveformMaxCut? "true" : "false"));

  c("DoPositionCut").Get(fDoPositionCut);
  MsgInfo(MsgLog::Form("\t-Apply Position Cut %s",fDoPositionCut? "true" : "false"));
  if (fDoPositionCut) {
    c("RadiusCutValueLow").Get(fRadiusCutValueLow);
    c("RadiusCutValueHigh").Get(fRadiusCutValueHigh);
    c("ZCutValueLow").Get(fZCutValueLow);
    c("ZCutValueHigh").Get(fZCutValueHigh);
    MsgInfo(MsgLog::Form("\t\t-Radius Cut Low %g high %g",fRadiusCutValueLow,fRadiusCutValueHigh));
    MsgInfo(MsgLog::Form("\t\t-Z Cut Low %g high %g",fZCutValueLow,fZCutValueHigh));
  }

  c("DoEnergyCut").Get(fDoEnergyCut);
  MsgInfo(MsgLog::Form("\t-Apply Energy Cut %s",fDoEnergyCut? "true" : "false"));
  if (fDoEnergyCut) {
    c("EnergyCutValueLow").Get(fEnergyCutValueLow);
    c("EnergyCutValueHigh").Get(fEnergyCutValueHigh);
    MsgInfo(MsgLog::Form("\t\t-Energy Cut Low %g high %g",fEnergyCutValueLow,fEnergyCutValueHigh));
  }

  c("DoVetoCut").Get(fDoVetoCut);
  MsgInfo(MsgLog::Form("\t-Apply Veto Cut %s",fDoVetoCut? "true" : "false"));
  if (fDoVetoCut) {
    c("NumVetoCut").Get(fNumVetoCut);
    c("ReverseVeto").Get(fReverseVeto);
    MsgInfo(MsgLog::Form("\t\t-Veto Cut %d Reverse %s",fNumVetoCut,fReverseVeto ? "true" : "false"));
  }

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

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMNa22Cuts::SetupOutFile()
 * \brief Setup the output file that gets saved
 ***********************************************/
void CCMNa22Cuts::SetupOutFile()
{
  if (fOutfile) {
    if (fOutfile->GetName() != fOutFileName) {
      MsgWarning(MsgLog::Form("New outfile name %s is different than old %s, this should not happen, sticking with ld",
            fOutFileName.c_str(),fOutfile->GetName()));
    }
    return;
  }

  fEnergy = 0;
  fTotalEnergy = 0;
  fLength = 0;
  fHits = 0;
  fNumPMTs = 0;
  fTime = 0;
  fVetoTop = 0;
  fVetoBottom = 0;
  fVetoCLeft = 0;
  fVetoCRight = 0;
  fVetoCFront = 0;
  fVetoCBack = 0;
  fVetoPromptTop = 0;
  fVetoPromptBottom = 0;
  fVetoPromptCLeft = 0;
  fVetoPromptCRight = 0;
  fVetoPromptCFront = 0;
  fVetoPromptCBack = 0;
  fWeight = 0;
  fX = 0;
  fY = 0;
  fZ = 0;
  fLargestPMTFraction = 0;
  fEpochSec = 0;
  fEpochNSSec = 0;
  fBCMTime = 0;
  fBCMLength = 0;
  fBCMIntegral = 0;
  fMaxPrevEnergy1600 = 0;
  fMaxPrevEnergy3200 = 0;
  fMaxPrevEnergy4800 = 0;
  fMaxPrevEnergy = 0;
  fMaxPrevEnergyTime1600 = 0;
  fMaxPrevEnergyTime3200 = 0;
  fMaxPrevEnergyTime4800 = 0;
  fMaxPrevEnergyTime = 0;
  fMaxPrevEnergyLength1600 = 0;
  fMaxPrevEnergyLength3200 = 0;
  fMaxPrevEnergyLength4800 = 0;
  fMaxPrevEnergyLength = 0;

  fPassedVetoCut = false;
  fPassedPositionCut = false;
  fPassedEnergyCut = false;
  fPassedPrevCut = false;
  fPassedLengthCut = false;
  fPassedTimeCut = false;
  fPassedNicenessCut = false;
  fPassedWaveformMaxCut = false;
  fPassedNumHitsCut = false;
  fPassedBCMCut = false;

  if (!fOutFileName.empty()) {
    fOutfile = gROOT->GetFile(fOutFileName.c_str());
    if (!fOutfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fOutFileName.c_str()));
      return;
    }
    fOutfile->cd();
    fTree = new TTree(fTreeName.c_str(),fTreeName.c_str());
    fTree->Branch("epochSec",&fEpochSec);
    fTree->Branch("epochNS",&fEpochNSSec);
    fTree->Branch("bcmTime",&fBCMTime);
    fTree->Branch("bcmLength",&fBCMLength);
    fTree->Branch("bcmIntegral",&fBCMIntegral);
    fTree->Branch("energy",&fEnergy);
    fTree->Branch("totalenergy",&fTotalEnergy);
    fTree->Branch("length",&fLength);
    fTree->Branch("hits",&fHits);
    fTree->Branch("numPMTs",&fNumPMTs);
    fTree->Branch("time",&fTime);
    fTree->Branch("vetoTop",&fVetoTop);
    fTree->Branch("vetoBottom",&fVetoBottom);
    fTree->Branch("vetoCLeft",&fVetoCLeft);
    fTree->Branch("vetoCRight",&fVetoCRight);
    fTree->Branch("vetoCFront",&fVetoCFront);
    fTree->Branch("vetoCBack",&fVetoCBack);
    fTree->Branch("vetoPromptTop",&fVetoPromptTop);
    fTree->Branch("vetoPromptBottom",&fVetoPromptBottom);
    fTree->Branch("vetoPromptCLeft",&fVetoPromptCLeft);
    fTree->Branch("vetoPromptCRight",&fVetoPromptCRight);
    fTree->Branch("vetoPromptCFront",&fVetoPromptCFront);
    fTree->Branch("vetoPromptCBack",&fVetoPromptCBack);
    fTree->Branch("weight",&fWeight);
    fTree->Branch("x",&fX);
    fTree->Branch("y",&fY);
    fTree->Branch("z",&fZ);
    fTree->Branch("largestPMTFraction",&fLargestPMTFraction);
    fTree->Branch("maxPrevEnergy1600",&fMaxPrevEnergy1600);
    fTree->Branch("maxPrevEnergy3200",&fMaxPrevEnergy3200);
    fTree->Branch("maxPrevEnergy4800",&fMaxPrevEnergy4800);
    fTree->Branch("maxPrevEnergy",&fMaxPrevEnergy);
    fTree->Branch("maxPrevEnergyTime1600",&fMaxPrevEnergyTime1600);
    fTree->Branch("maxPrevEnergyTime3200",&fMaxPrevEnergyTime3200);
    fTree->Branch("maxPrevEnergyTime4800",&fMaxPrevEnergyTime4800);
    fTree->Branch("maxPrevEnergyTime",&fMaxPrevEnergyTime);
    fTree->Branch("maxPrevEnergyLength1600",&fMaxPrevEnergyLength1600);
    fTree->Branch("maxPrevEnergyLength3200",&fMaxPrevEnergyLength3200);
    fTree->Branch("maxPrevEnergyLength4800",&fMaxPrevEnergyLength4800);
    fTree->Branch("maxPrevEnergyLength",&fMaxPrevEnergyLength);

    fTree->Branch("passedVetoCut",&fPassedVetoCut);
    fTree->Branch("passedPositionCut",&fPassedPositionCut);
    fTree->Branch("passedEnergyCut",&fPassedEnergyCut);
    fTree->Branch("passedPrevCut",&fPassedPrevCut);
    fTree->Branch("passedLengthCut",&fPassedLengthCut);
    fTree->Branch("passedTimeCut",&fPassedTimeCut);
    fTree->Branch("passedNicenessCut",&fPassedNicenessCut);
    fTree->Branch("passedWaveformMaxCut",&fPassedWaveformMaxCut);
    fTree->Branch("passedNumHitsCut",&fPassedNumHitsCut);
    fTree->Branch("passedBCMCut",&fPassedBCMCut);
  }
}

//------------------------------------------------------
CCMResult_t CCMNa22Cuts::EndOfJob()
{
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  if (fOutfile != nullptr) {
    fOutfile->cd();
    fTree->Write();
  }

  MsgInfo(MsgLog::Form("Number of events %ld passed from %ld total (config %s)",fNumFinalEvents,fNumInitEvents,fCfgVersion.c_str()));
  MsgInfo("Failed Breakdown:");
  MsgInfo(MsgLog::Form("- Failed BCM %ld",fNumFailedBCMCut));
  MsgInfo(MsgLog::Form("- Failed Time %ld",fNumFailedTimeCut));
  MsgInfo(MsgLog::Form("- Failed NumHits %ld",fNumFailedNumHitsCut));
  MsgInfo(MsgLog::Form("- Failed Prev %ld",fNumFailedPrevCut));
  MsgInfo(MsgLog::Form("- Failed Veto %ld",fNumFailedVetoCut));
  MsgInfo(MsgLog::Form("- Failed Position %ld",fNumFailedPositionCut));
  MsgInfo(MsgLog::Form("- Failed Niceness %ld",fNumFailedNicenessCut));
  MsgInfo(MsgLog::Form("- Failed WaveformMax %ld",fNumFailedWaveformMaxCut));
  MsgInfo(MsgLog::Form("- Failed Length %ld",fNumFailedLengthCut));
  MsgInfo(MsgLog::Form("- Failed Energy %ld",fNumFailedEnergyCut));

  return kCCMSuccess;
}

/*!**********************************************
 * \fn void CCMNa22Cuts::RecalculatePosition(const SimplifiedEvent & simplifiedEvent, const double fitLength, double & x, double & y, double & z)
 * \brief Recalcualte the position of the event based on the \p fitLength value passed to the function
 * \param[in] simplifiedEvent The #SimplifiedEvent to recalculate the position for
 * \param[in] fitLength The length of time to use to calculate the position
 * \param[out] x The resulting x position
 * \param[out] y The resulting y position
 * \param[out] z The resulting z position
 *
 * Recalculate the position of the event based on the \p fitLength value passed to the function
 * The current algorithm is a charge fitter
 ***********************************************/
void CCMNa22Cuts::RecalculatePosition(const SimplifiedEvent & simplifiedEvent, 
    const double fitLength, double & x, double & y, double & z)
{
  int key = 0;
  std::vector<float> pmtInt;
  std::vector<int> pmtCount;

  const size_t kNumBins = fitLength/2e-3;

  std::vector<TVector3> pmtContribution;

  simplifiedEvent.ResetPMTWaveformItr();
  TVector3 pos;
  double totalWeight = 0.0;
  do
  {
    simplifiedEvent.GetPMTWaveform(key,pmtInt,pmtCount);
    if (!PMTInfoMap::IsActive(key)) {
      MsgWarning("PMT not Active: should not be here because this should have already been checked");
      continue;
    }

    //skipping every 16th channel.                                                                                                                           
    if((key+1)%16 == 0){
      continue;
    }

    auto pmtInfo = PMTInfoMap::GetPMTInfo(key);
    if (!pmtInfo) {
      MsgWarning("PMTInformation does not exists: should not be here because this should have already been checked");
      continue;
    }

    if (pmtInfo->IsVeto()) {
      continue;
    }
    if(key >= 240){
      continue;
    }

    auto start = pmtInt.begin();
    auto end = pmtInt.end();
    if (kNumBins < pmtInt.size()) {
      end = pmtInt.begin()+kNumBins;
    }
    float pmtEnergy = std::accumulate(start,end,0.f);

    pos += *(pmtInfo->GetPosition())*pmtEnergy*pmtEnergy;
    pmtContribution.push_back(*(pmtInfo->GetPosition())*pmtEnergy*pmtEnergy);
    totalWeight += pmtEnergy*pmtEnergy;

  } while (simplifiedEvent.NextPMTWaveform());

  pos *= 1.0/totalWeight;

  x = pos.X();
  y = pos.Y();
  z = pos.Z();

  if (x == 0 && y == 0 && z == 0) {
    for (const auto & vec : pmtContribution) {
      vec.Print();
    }
    MsgInfo(MsgLog::Form("Total weight = %.3f",totalWeight));
  }
  return;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedVetoCut(int vetoTotal)
{
  if ((vetoTotal >= fNumVetoCut && !fReverseVeto) || 
      (fReverseVeto && vetoTotal < fNumVetoCut)) {
    return false;
  } // end veto cut check

  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedPositionCut(double x, double y, double z)
{
  double radius = std::sqrt(x*x+y*y);
  if (radius < fRadiusCutValueLow || radius > fRadiusCutValueHigh ||
      z < fZCutValueLow || z > fZCutValueHigh) {
    return false;
  } // end position cut check

  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedEnergyCut(double energy)
{
  if (energy < fEnergyCutValueLow || energy > fEnergyCutValueHigh) {
    return false;
  }// end energy cut check

  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedPrevCut(const double kStartTime, const long kStartingIndex, const double kEnergy)
{
  if (fEnergyLengthFile == nullptr) {
    fEnergyLengthFile = std::make_shared<TFile>("$CCMPROJECT/calibrationFiles/2019/energyLength_0p2Thresh.root","READ");
    fEnergyLengthFile->GetObject("eneLength_pfx",fEnergyLengthProf);
  }
  std::vector<float> currentWF;
  //std::vector<float> totalWF;
  bool passed = true;
  fMaxPrevEnergy1600 = 0;
  fMaxPrevEnergy3200 = 0;
  fMaxPrevEnergy4800 = 0;
  fMaxPrevEnergy = 0;
  fMaxPrevEnergyTime1600 = 0;
  fMaxPrevEnergyTime3200 = 0;
  fMaxPrevEnergyTime4800 = 0;
  fMaxPrevEnergyTime = 0;
  fMaxPrevEnergyLength1600 = 0;
  fMaxPrevEnergyLength3200 = 0;
  fMaxPrevEnergyLength4800 = 0;
  fMaxPrevEnergyLength = 0;
  for (long e2 = kStartingIndex-1; e2 >= 0; --e2) {
    auto simplifiedEvent2 = fEvents->GetSimplifiedEvent(e2);

    if (simplifiedEvent2.GetEventFinderMethod() != fEventFinderID ||
        simplifiedEvent2.GetAccumWaveformMethod() != fAccumWaveformMethodID) {
      continue;
    }

    double promptHits2 = simplifiedEvent2.GetNumCoated(true);
    double st2 = simplifiedEvent2.GetStartTime();
    double length = simplifiedEvent2.GetLength();

    //currentWF = simplifiedEvent2.GetWaveformInt();

    //if (toatlWF.size() == 0) {
    //  std::assign(currentWF.begin(),currentWF.end());
    //} else {
    //  double numZeros = std::fabs(kStartTime - st2+length);
    //  if (numZeros > 0) {
    //    totalWF.insert(totalWF.begin(),numZeros,0.0);
    //  }
    //  totalWF.insert(totalWF.begin(),currentWF.begin(),currentWF.end());
    //}

    //if (e2 == kStartingIndex) {
    //  continue;
    //}

    if (promptHits2 < 6) {
      //continue;
    }

    //double x2 = simplifiedEvent2.GetXPosition();
    //double y2 = simplifiedEvent2.GetYPosition();
    //double z2 = simplifiedEvent2.GetZPosition();
    //if (fDoPositionCut) {
    //  bool passedPosition = PassedPositionCut(simplifiedEvent2);
    //}

    double stDiff = std::fabs(st2 - kStartTime);

    double energy2 = simplifiedEvent2.GetIntegralTank(true);

    if (stDiff < 4800) {
      if (fMaxPrevEnergy4800 < energy2) {
        fMaxPrevEnergy4800 = energy2;
        fMaxPrevEnergyTime4800 = st2;
        fMaxPrevEnergyLength4800 = length;
      }
      if (stDiff < 3200) {
        if (fMaxPrevEnergy3200 < energy2) {
          fMaxPrevEnergy3200 = energy2;
          fMaxPrevEnergyTime3200 = st2;
          fMaxPrevEnergyLength3200 = length;
        }
        if (stDiff < 1600) {
          if (fMaxPrevEnergy1600 < energy2) {
            fMaxPrevEnergy1600 = energy2;
            fMaxPrevEnergyTime1600 = st2;
            fMaxPrevEnergyLength1600 = length;
          }
        }
      }
    }

    if (fMaxPrevEnergyLength != 0) {
      continue;
    }

    if (energy2 >= kEnergy/2.0) {
      if (std::fabs(st2+length - kStartTime) <= fPrevCutHoldOff) {
        fMaxPrevEnergy = energy2;
        fMaxPrevEnergyTime = st2;
        fMaxPrevEnergyLength = length;
        passed = false;
      } else {
        int bin2 = fEnergyLengthProf->FindBin(energy2);
        double content = fEnergyLengthProf->GetBinContent(bin2);
        double error = fEnergyLengthProf->GetBinError(bin2);
        if (stDiff < content+error*3.0+fPrevCutHoldOff) {
          fMaxPrevEnergy = energy2;
          fMaxPrevEnergyTime = st2;
          fMaxPrevEnergyLength = length;
          passed = false;
        }
      }
    } else if (length > 90) {
      currentWF = simplifiedEvent2.GetWaveformInt();
      double firstInt = currentWF.front();
      double newST = st2;
      size_t loc = 45;
      for (; loc < currentWF.size(); ++loc) {
        energy2 -= firstInt;
        energy2 += currentWF.at(loc);
        firstInt = currentWF.at(loc-44);
        newST += Utility::fgkBinWidth;
        if (energy2 < kEnergy/2.0) {
          continue;
        }
        int bin2 = fEnergyLengthProf->FindBin(energy2);
        double content = fEnergyLengthProf->GetBinContent(bin2);
        double error = fEnergyLengthProf->GetBinError(bin2);
        if (std::fabs(newST-kStartTime) < content+error*3.0+fPrevCutHoldOff) {
          fMaxPrevEnergy = energy2;
          fMaxPrevEnergyTime = newST;
          fMaxPrevEnergyLength = length;
          passed = false;
          break;
        }
      }
    }

    continue;

    if (stDiff < fPrevCutTime) {
      passed = false;
    } else {
      break;
    }
  } // end for e2

  return passed;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedLengthCut(double length)
{
  if (length < fLengthCutValueLow || length > fLengthCutValueHigh) {
    return false;
  } // end length cut condition

  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedTimeCut(double time)
{
  if (time < fTimeCutValueLow || time > fTimeCutValueHigh) {
    return false;
  } // end time cut condition

  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedNicenessCut(double energy, int hits, double largestFrac)
{
  double uniform = energy/static_cast<double>(hits);
  double nice = (largestFrac - uniform)/uniform;
  if (nice < fNicenessCutValueLow || nice > fNicenessCutValueHigh) {
    return false;
  }

  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedWaveformMaxCut(double length, double maxWFTime)
{
  double promptLength = std::min(90.0,length);
  if (maxWFTime > promptLength) {
    return false;
  }
  return true;
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedNumHitsCut(double hits)
{
  return !(hits < 6);
}

/*!**********************************************
 ***********************************************/
bool CCMNa22Cuts::PassedBCMCut(double bcmTime, double bcmLength, double bcmIntegral)
{
  //bcmTime < -345.77-10.54*3.0 || bcmTime > -345.77+10.54*3.0 ||
  //bcmLength < 576.0-23.5*3.0 || bcmLength > 576.0+23.5*3.0
  
  if (bcmTime < fBCMTimeLowCut || bcmTime > fBCMTimeHighCut ||
      bcmLength < fBCMLengthLowCut || bcmLength > fBCMLengthHighCut ||
      bcmIntegral < fBCMIntegralLowCut || bcmIntegral > fBCMIntegralHighCut) {
    return false;
  }

  return true;
}

/*!**********************************************
 * \fn void CCMNa22Cuts::RecalculateStartTime(const SimplifiedEvent & events, double & st, double & charge, double & hits, double & length) 
 * \brief Recalcualte the start time of an event given a different method
 * \param[in] events The event to look at
 * \param[out] st The new start time
 * \param[out] charge The new integral of the event
 * \param[out] hits The new number of hits in the event
 * \param[out] length The new length of the event
 ***********************************************/
void CCMNa22Cuts::RecalculateStartTime(const SimplifiedEvent & events, double & st, double & charge, double & hits, double & length) 
{
  const std::vector<float> & waveform = events.GetWaveformInt();
  const std::vector<float> & waveformCount = events.GetWaveformCount();
  bool crossed = false;
  size_t stCount = 0;
  charge = 0;
  hits = 0;
  for (size_t count = 0; count < waveform.size(); ++count) {
    if (!crossed && waveform[count] >= 0.1) {
      crossed = true;
      st = events.GetStartTime()+static_cast<double>(count)*Utility::fgkBinWidth;
      stCount = count;
      charge = 0;
      hits = 0;
      length = (waveform.size() - count)*Utility::fgkBinWidth;
    }

    if (crossed && count < stCount+45) {
      charge += waveform[count];
      hits += waveformCount[count];
    } else if (crossed && count >= stCount+45) {
      break;
    }
    ++count;
  }

  return;

}

/*!**********************************************
 * \fn int CCMNa22Cuts::WaveformMaxPosition(const SimplifiedEvent & simplifiedEvent)
 * \brief Find where the waveform of the event is the largest
 * \param[in] simplifiedEvent The #SimplifiedEvent to recalculate the position for
 * \return The position of when the waveform is the largest
 ***********************************************/
/*
   int CCMNa22Cuts::WaveformMaxPosition(const SimplifiedEvent & simplifiedEvent)
   {
   auto waveform = simplifiedEvent.GetWaveformInt();
   auto max_iter = std::max_element(waveform.begin(),waveform.end());
   return std::distance(waveform.begin(),max_iter);
   }
   */

/*
   std::vector<std::shared_ptr<TF1>> weightFunction;
   for (int i=0; i < 4; ++i) {
   weightFunction.push_back(std::make_shared<TF1>(Form("weightFunction%d",i),"[0]*TMath::Exp(-(x-[1])*[2])",-9.92,6.08));
   switch (i) {
   case 0:
   weightFunction.back()->FixParameter(0,8.63796e-7);
   weightFunction.back()->FixParameter(1,2.45729e1);
   weightFunction.back()->FixParameter(2,9.84429e-2);
   break;
   case 1:
   weightFunction.back()->FixParameter(0,7.90765e-5);
   weightFunction.back()->FixParameter(1,3.73707e1);
   weightFunction.back()->FixParameter(2,8.74439e-2);
   break;
   case 2:
   weightFunction.back()->FixParameter(0,2.14018e-4);
   weightFunction.back()->FixParameter(1,3.70502e1);
   weightFunction.back()->FixParameter(2,8.85969e-2);
   break;
   case 3:
   weightFunction.back()->FixParameter(0,6.36795e-6);
   weightFunction.back()->FixParameter(1,3.16871e1);
   weightFunction.back()->FixParameter(2,9.26635e-2);
   break;
   default: break;
   }// end switch
   } // end for i < 4;
   */

