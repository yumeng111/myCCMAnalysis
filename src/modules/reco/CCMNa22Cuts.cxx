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

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMNa22Cuts.h"
#include "CCMModuleTable.h"

#include "Events.h"
#include "SimplifiedEvent.h"
#include "MsgLog.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <memory>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TVector3.h"

//See CCMModuleTable for info
MODULE_DECL(CCMNa22Cuts);

//_______________________________________________________________________________________
CCMNa22Cuts::CCMNa22Cuts(const char* version) 
  : CCMModule("CCMNa22Cuts"),
    fEvents(nullptr),
    fEventFinderID(kCCMDynamicLengthEventID),
    fAccumWaveformMethodID(kCCMAccumWaveformTriangleID),
    fOutFileName(""),
    fTreeName("tree"),
    fRemovePrimaryEvents(true),
    fRemoveOtherEvents(false),
    fDoVetoCut(false),
    fDoPositionCut(false),
    fDoEnergyCut(false),
    fDoPrevCut(false),
    fDoLengthCut(false),
    fDoTimeCut(false),
    fDoNicenessCut(false),
    fReverseVeto(false),
    fWaveformMaxCut(false),
    fPassedVetoCut(true),
    fPassedPositionCut(true),
    fPassedEnergyCut(true),
    fPassedPrevCut(true),
    fPassedLengthCut(true),
    fPassedTimeCut(true),
    fPassedNicenessCut(true),
    fPassedWaveformMaxCut(true),
    fPassedNumHitsCut(true),
    fNumVetoCut(3),
    fRadiusCutValueLow(0),
    fRadiusCutValueHigh(0.95),
    fZCutValueLow(-0.4),
    fZCutValueHigh(0.4),
    fEnergyCutValueLow(0),
    fEnergyCutValueHigh(std::numeric_limits<double>::max()),
    fPrevCutTime(1600),
    fLengthCutValueLow(-std::numeric_limits<double>::max()),
    fLengthCutValueHigh(std::numeric_limits<double>::max()),
    fTimeCutValueLow(-7000),
    fTimeCutValueHigh(4000),
    fNicenessCutValueLow(0),
    fNicenessCutValueHigh(2.5),
    fOutfile(nullptr),
    fTree(nullptr),
    fEnergy(0),
    fLength(0),
    fHits(0),
    fNumPMTs(0),
    fTime(0),
    fVetoTop(0),
    fVetoBottom(0),
    fVetoCLeft(0),
    fVetoCRight(0),
    fVetoCFront(0),
    fVetoCBack(0),
    fVetoPromptTop(0),
    fVetoPromptBottom(0),
    fVetoPromptCLeft(0),
    fVetoPromptCRight(0),
    fVetoPromptCFront(0),
    fVetoPromptCBack(0),
    fWeight(0),
    fX(0),
    fY(0),
    fZ(0),
    fLargestPMTFraction(0),
    fEpochSec(0)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMNa22Cuts::CCMNa22Cuts(const CCMNa22Cuts& clufdr) 
: CCMModule(clufdr),
  fEvents(clufdr.fEvents),
  fEventFinderID(clufdr.fEventFinderID),
  fAccumWaveformMethodID(clufdr.fAccumWaveformMethodID),
  fOutFileName(clufdr.fOutFileName),
  fTreeName(clufdr.fTreeName),
  fRemovePrimaryEvents(clufdr.fRemovePrimaryEvents),
  fRemoveOtherEvents(clufdr.fRemoveOtherEvents),
  fDoVetoCut(clufdr.fDoVetoCut),
  fDoPositionCut(clufdr.fDoPositionCut),
  fDoEnergyCut(clufdr.fDoEnergyCut),
  fDoPrevCut(clufdr.fDoPrevCut),
  fDoLengthCut(clufdr.fDoLengthCut),
  fDoTimeCut(clufdr.fDoTimeCut),
  fDoNicenessCut(clufdr.fDoNicenessCut),
  fReverseVeto(clufdr.fReverseVeto),
  fWaveformMaxCut(clufdr.fWaveformMaxCut),
  fPassedVetoCut(clufdr.fPassedVetoCut),
  fPassedPositionCut(clufdr.fPassedPositionCut),
  fPassedEnergyCut(clufdr.fPassedEnergyCut),
  fPassedPrevCut(clufdr.fPassedPrevCut),
  fPassedLengthCut(clufdr.fPassedLengthCut),
  fPassedTimeCut(clufdr.fPassedTimeCut),
  fPassedNicenessCut(clufdr.fPassedNicenessCut),
  fPassedWaveformMaxCut(clufdr.fPassedWaveformMaxCut),
  fPassedNumHitsCut(clufdr.fPassedNumHitsCut),
  fNumVetoCut(clufdr.fNumVetoCut),
  fRadiusCutValueLow(clufdr.fRadiusCutValueLow),
  fRadiusCutValueHigh(clufdr.fRadiusCutValueHigh),
  fZCutValueLow(clufdr.fZCutValueLow),
  fZCutValueHigh(clufdr.fZCutValueHigh),
  fEnergyCutValueLow(clufdr.fEnergyCutValueLow),
  fEnergyCutValueHigh(clufdr.fEnergyCutValueHigh),
  fPrevCutTime(clufdr.fPrevCutTime),
  fLengthCutValueLow(clufdr.fLengthCutValueLow),
  fLengthCutValueHigh(clufdr.fLengthCutValueHigh),
  fTimeCutValueLow(clufdr.fTimeCutValueLow),
  fTimeCutValueHigh(clufdr.fTimeCutValueHigh),
  fNicenessCutValueLow(clufdr.fNicenessCutValueLow),
  fNicenessCutValueHigh(clufdr.fNicenessCutValueHigh),
  fOutfile(clufdr.fOutfile),
  fTree(clufdr.fTree),
  fEnergy(clufdr.fEnergy),
  fLength(clufdr.fLength),
  fHits(clufdr.fHits),
  fNumPMTs(clufdr.fNumPMTs),
  fTime(clufdr.fTime),
  fVetoTop(clufdr.fVetoTop),
  fVetoBottom(clufdr.fVetoBottom),
  fVetoCLeft(clufdr.fVetoCLeft),
  fVetoCRight(clufdr.fVetoCRight),
  fVetoCFront(clufdr.fVetoCFront),
  fVetoCBack(clufdr.fVetoCBack),
  fVetoPromptTop(clufdr.fVetoPromptTop),
  fVetoPromptBottom(clufdr.fVetoPromptBottom),
  fVetoPromptCLeft(clufdr.fVetoPromptCLeft),
  fVetoPromptCRight(clufdr.fVetoPromptCRight),
  fVetoPromptCFront(clufdr.fVetoPromptCFront),
  fVetoPromptCBack(clufdr.fVetoPromptCBack),
  fWeight(clufdr.fWeight),
  fX(clufdr.fX),
  fY(clufdr.fY),
  fZ(clufdr.fZ),
  fLargestPMTFraction(clufdr.fLargestPMTFraction),
  fEpochSec(clufdr.fEpochSec)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMNa22Cuts::~CCMNa22Cuts()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMNa22Cuts::ProcessTrigger()
{
  // since the class does not own the TFile, always check to see if it is present
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  
  std::vector<size_t> locationsToRemove;

  const size_t kNumEvents = fEvents->GetNumEvents();
  fEpochSec = fEvents->GetComputerSecIntoEpoch();

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

    double st = simplifiedEvent.GetStartTime();
    fLength = simplifiedEvent.GetLength();

    // time of length is in ns
    double promptLength = std::min(90.0,fLength);

    fEnergy = simplifiedEvent.GetIntegralTank(true);
    fHits = simplifiedEvent.GetNumTank(true);
    fNumPMTs = simplifiedEvent.GetPMTHits();
    //MsgInfo(MsgLog::Form("e %zu energy %.2f hits %f numPTs %d length %.3f",e,fEnergy,fHits,fNumPMTs,fLength));

    fX = simplifiedEvent.GetXPosition();
    fY = simplifiedEvent.GetYPosition();
    fZ = simplifiedEvent.GetZPosition();

    if (fWaveformMaxCut) {
      if (simplifiedEvent.GetMaxAccumWaveformTime() > promptLength) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedWaveformMaxCut = false;
        //continue;
      }
    }


    fTime = st;

    if (fDoTimeCut) {
      if (st < fTimeCutValueLow || st > fTimeCutValueHigh) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedTimeCut = false;
        //continue;
      } // end time cut condition
    } // end fDoTimeCut

    double et = st + fLength;
    if (et <= st) {
      if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
      //continue;
    }

    if (fDoLengthCut) {
      if (fLength < fLengthCutValueLow || fLength > fLengthCutValueHigh) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedLengthCut = false;
        //continue;
      } // end length cut condition
    } // end if fDoLengthCut

    double radius = std::sqrt(fX*fX+fY*fY);
    if (fDoPositionCut) {
      if (radius < fRadiusCutValueLow || radius > fRadiusCutValueHigh ||
          fZ < fZCutValueLow || fZ > fZCutValueHigh) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedPositionCut = false;
        //continue;
      } // end position cut check
    } // end if fDoPositionCut

    if (fDoEnergyCut) {
      if (fEnergy < fEnergyCutValueLow || fEnergy > fEnergyCutValueHigh) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedEnergyCut = false;
        //continue;
      }// end energy cut check
    }// end if fDoEnergyCut

    if (fHits < 3) {
      if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
      fPassedNumHitsCut = false;
      //continue;
    }

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
    int vetoTotal = simplifiedEvent.GetNumVeto(true);

    fLargestPMTFraction = simplifiedEvent.GetLargestPMTFraction();
    if (fDoNicenessCut) {
      double uniform = fEnergy/static_cast<double>(fNumPMTs);
      double nice = (fLargestPMTFraction - uniform)/uniform;
      if (nice < fNicenessCutValueLow || nice > fNicenessCutValueHigh) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedNicenessCut = false;
        //continue;
      }
    } // end if fDoNicenessCut

    if (fDoVetoCut) {
      if ((vetoTotal >= fNumVetoCut && !fReverseVeto) || 
          (fReverseVeto && vetoTotal < fNumVetoCut)) {
        if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
          locationsToRemove.emplace_back(e);
        }
        fPassedVetoCut = false;
        //continue;
      } // end veto cut check
    }// end if fDoVetoCut


    bool reject = false;
    for (long e2 = e-1; e2 >= 0 && fDoPrevCut; --e2) {
      auto simplifiedEvent2 = fEvents->GetSimplifiedEvent(e2);

      if (simplifiedEvent2.GetEventFinderMethod() != fEventFinderID ||
          simplifiedEvent2.GetAccumWaveformMethod() != fAccumWaveformMethodID) {
        continue;
      }

      double promptHits2 = simplifiedEvent2.GetNumCoated(true);
      double st2 = simplifiedEvent2.GetStartTime();

      if (promptHits2 < 3) {
        //continue;
      }

      double x2 = simplifiedEvent2.GetXPosition();
      double y2 = simplifiedEvent2.GetYPosition();
      double z2 = simplifiedEvent2.GetZPosition();
      if (fDoPositionCut) {
        double radius2 = std::sqrt(x2*x2+y2*y2);
        if (radius2 < fRadiusCutValueLow || radius2 > fRadiusCutValueHigh ||
            z2 < fZCutValueLow || z2 > fZCutValueHigh) {
          //continue;
        } // end position cut check
      }

      double stDiff = std::fabs(st2 - st);

      //MsgInfo(MsgLog::Form("Event (%ld) ST = %.0f Previous (%ld) ST = %.0f Diff = %.0f",e,st,e2,st2,stDiff));

      //double pIDiff = promptIntegral - promptIntegral2;
      if (stDiff < fPrevCutTime) {
        reject = true;
        break;
      }
    } // end for e2

    if (reject) {
      fPassedPrevCut = false;
      if (std::find(locationsToRemove.begin(),locationsToRemove.end(),e) == locationsToRemove.end() && fRemovePrimaryEvents) {
        locationsToRemove.emplace_back(e);
      }
      //continue;
    }

    fWeight = 1.0;// /weightFunction.at(energyRegion)->Eval(st);
    //MsgInfo(MsgLog::Form("st %f weighFunction %g weight %g",st,weightFunction->Eval(st),weight));

    if (fOutfile != nullptr) {
      //MsgInfo("Filling the tree");
      fOutfile->cd();
      fTree->Fill();
    }
  } // end e

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
  fEventFinderID = Utility::ConvertStringToCCMEventFinderID(tempString);
  c("AccumWaveformMethodID").Get(tempString);
  fAccumWaveformMethodID = Utility::ConvertStringToCCMAccumWaveformMethod(tempString);

  c("RemovePrimaryEvents").Get(fRemovePrimaryEvents);
  c("RemoveOtherEvents").Get(fRemoveOtherEvents);
  c("TreeName").Get(fTreeName);

  c("DoLengthCut").Get(fDoLengthCut);
  if (fDoLengthCut) {
    c("LengthCutValueLow").Get(fLengthCutValueLow);
    c("LengthCutValueHigh").Get(fLengthCutValueHigh);
  }

  c("DoTimeCut").Get(fDoTimeCut);
  if (fDoTimeCut) {
    c("TimeCutValueLow").Get(fTimeCutValueLow);
    c("TimeCutValueHigh").Get(fTimeCutValueHigh);
  }

  c("DoNicenessCut").Get(fDoNicenessCut);
  if (fDoNicenessCut) {
    c("NicenessCutValueLow").Get(fNicenessCutValueLow);
    c("NicenessCutValueHigh").Get(fNicenessCutValueHigh);
  }

  c("DoPrevCut").Get(fDoPrevCut);
  if (fDoPrevCut) {
    c("PrevCutTime").Get(fPrevCutTime);
  }

  c("WaveformMaxCut").Get(fWaveformMaxCut);

  c("DoPositionCut").Get(fDoPositionCut);
  if (fDoPositionCut) {
    c("RadiusCutValueLow").Get(fRadiusCutValueLow);
    c("RadiusCutValueHigh").Get(fRadiusCutValueHigh);
    c("ZCutValueLow").Get(fZCutValueLow);
    c("ZCutValueHigh").Get(fZCutValueHigh);
  }

  c("DoEnergyCut").Get(fDoEnergyCut);
  if (fDoEnergyCut) {
    c("EnergyCutValueLow").Get(fEnergyCutValueLow);
    c("EnergyCutValueHigh").Get(fEnergyCutValueHigh);
  }

  c("DoVetoCut").Get(fDoVetoCut);
  if (fDoVetoCut) {
    c("NumVetoCut").Get(fNumVetoCut);
    c("ReverseVeto").Get(fReverseVeto);
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
  fPassedVetoCut = false;
  fPassedPositionCut = false;
  fPassedEnergyCut = false;
  fPassedPrevCut = false;
  fPassedLengthCut = false;
  fPassedTimeCut = false;
  fPassedNicenessCut = false;
  fPassedWaveformMaxCut = false;
  fPassedNumHitsCut = false;

  if (!fOutFileName.empty()) {
    fOutfile = gROOT->GetFile(fOutFileName.c_str());
    if (!fOutfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fOutFileName.c_str()));
      return;
    }
    fOutfile->cd();
    fTree = new TTree(fTreeName.c_str(),fTreeName.c_str());
    fTree->Branch("epochSec",&fEpochSec);
    fTree->Branch("energy",&fEnergy);
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
    fTree->Branch("passedVetoCut",&fPassedVetoCut);
    fTree->Branch("passedPositionCut",&fPassedPositionCut);
    fTree->Branch("passedEnergyCut",&fPassedEnergyCut);
    fTree->Branch("passedPrevCut",&fPassedPrevCut);
    fTree->Branch("passedLengthCut",&fPassedLengthCut);
    fTree->Branch("passedTimeCut",&fPassedTimeCut);
    fTree->Branch("passedNicenessCut",&fPassedNicenessCut);
    fTree->Branch("passedWaveformMaxCut",&fPassedWaveformMaxCut);
    fTree->Branch("passedNumHitsCut",&fPassedNumHitsCut);
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
    auto pmtInfo = PMTInfoMap::GetPMTInfo(key);
    if (!pmtInfo) {
      MsgWarning("PMTInformation does not exists: should not be here because this should have already been checked");
      continue;
    }
    if (pmtInfo->IsVeto()) {
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

