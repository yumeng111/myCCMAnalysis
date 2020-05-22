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

/*
 * Fit results to weight the event based on when it occured in the DAQ window
 * Function is [0]*TMath::Exp(-(x-[1])*[2])
 * FCN=37.4258 FROM HESSE     STATUS=NOT POSDEF     16 CALLS         341 TOTAL
 * EDM=1.13757e-07    STRATEGY= 1      ERR MATRIX NOT POS-DEF
 * EXT PARAMETER                APPROXIMATE        STEP         FIRST
 * NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
 * 1  p0           8.63796e-07   2.06143e-07   8.85969e-12  -3.24369e+04
 * 2  p1           2.45729e+01   2.42758e+00   1.04140e-04  -2.75826e-03
 * 3  p2           9.84429e-02   6.69477e-03   3.60225e-07  -7.89168e-01
 *
 * FCN=69.1019 FROM HESSE     STATUS=NOT POSDEF     16 CALLS         436 TOTAL
 * EDM=1.63512e-10    STRATEGY= 1      ERR MATRIX NOT POS-DEF
 * EXT PARAMETER                APPROXIMATE        STEP         FIRST
 * NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
 * 1  p0           7.90765e-05   1.34883e-06   7.65338e-11  -1.45708e+02
 * 2  p1           3.73707e+01   1.95086e-01   1.10682e-05  -1.00752e-03
 * 3  p2           8.74439e-02   3.66805e-04   2.34528e-08  -5.06801e-01
 *
 * FCN=53.9903 FROM HESSE     STATUS=NOT POSDEF     16 CALLS         418 TOTAL
 * EDM=3.36202e-08    STRATEGY= 1      ERR MATRIX NOT POS-DEF
 * EXT PARAMETER                APPROXIMATE        STEP         FIRST
 * NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
 * 1  p0           2.14018e-04   2.19993e-06   1.10589e-10  -1.69884e+03
 * 2  p1           3.70502e+01   1.16125e-01   8.83345e-06  -3.21872e-02
 * 3  p2           8.85283e-02   2.22341e-04   2.11068e-08  -1.48715e+01
 * 
 * FCN=40.4078 FROM HESSE     STATUS=NOT POSDEF     16 CALLS         380 TOTAL
 * EDM=1.50743e-10    STRATEGY= 1      ERR MATRIX NOT POS-DEF
 * EXT PARAMETER                APPROXIMATE        STEP         FIRST
 * NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
 * 1  p0           6.36795e-06   4.43500e-07   1.95030e-11  -5.56355e+02
 * 2  p1           3.16871e+01   7.51928e-01   3.30514e-05  -3.28286e-04
 * 3  p2           9.26635e-02   1.67588e-03   8.60329e-08  -1.26971e-01
 */


//See CCMModuleTable for info
MODULE_DECL(CCMNa22Cuts);

//_______________________________________________________________________________________
CCMNa22Cuts::CCMNa22Cuts(const char* version) 
  : CCMModule("CCMNa22Cuts"),
    fEvents(nullptr),
    fOutFileName(""),
    fDoVetoCut(false),
    fDoPositionCut(false),
    fDoEnergyCut(false),
    fDoPrevCut(false),
    fShiftTime(false),
    fDoLengthCut(false),
    fDoTimeCut(false),
    fDoNicenessCut(false),
    fReverseVeto(false),
    fNumVetoCut(std::numeric_limits<int>::max()),
    fRadiusCutValueLow(-std::numeric_limits<double>::max()),
    fZCutValueLow(-std::numeric_limits<double>::max()),
    fRadiusCutValueHigh(std::numeric_limits<double>::max()),
    fZCutValueHigh(std::numeric_limits<double>::max()),
    fEnergyCutValueLow(-std::numeric_limits<double>::max()),
    fEnergyCutValueHigh(std::numeric_limits<double>::max()),
    fPrevCutTime(-std::numeric_limits<double>::max()),
    fLengthCutValueLow(-std::numeric_limits<double>::max()),
    fLengthCutValueHigh(std::numeric_limits<double>::max()),
    fTimeCutValueLow(-std::numeric_limits<double>::max()),
    fTimeCutValueHigh(std::numeric_limits<double>::max()),
    fNicenessCutValueLow(-std::numeric_limits<double>::max()),
    fNicenessCutValueHigh(std::numeric_limits<double>::max()),
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
  fOutFileName(clufdr.fOutFileName),
  fDoVetoCut(clufdr.fDoVetoCut),
  fDoPositionCut(clufdr.fDoPositionCut),
  fDoEnergyCut(clufdr.fDoEnergyCut),
  fDoPrevCut(clufdr.fDoPrevCut),
  fShiftTime(clufdr.fShiftTime),
  fDoLengthCut(clufdr.fDoLengthCut),
  fDoTimeCut(clufdr.fDoTimeCut),
  fDoNicenessCut(clufdr.fDoNicenessCut),
  fReverseVeto(clufdr.fReverseVeto),
  fNumVetoCut(clufdr.fNumVetoCut),
  fRadiusCutValueLow(clufdr.fRadiusCutValueLow),
  fZCutValueLow(clufdr.fZCutValueLow),
  fRadiusCutValueHigh(clufdr.fRadiusCutValueHigh),
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
CCMResult_t CCMNa22Cuts::ProcessEvent()
{
  // since the class does not own the TFile, always check to see if it is present
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  
  std::vector<size_t> locationsToRemove;

  const double kPromptLength = 0.090;

  const size_t kNumEvents = fEvents->GetNumEvents();
  fEpochSec = fEvents->GetComputerSecIntoEpoch();
  for (size_t e = 0; e < kNumEvents; ++e) {
    auto simplifiedEvent = fEvents->GetSimplifiedEvent(e);

    double st = simplifiedEvent.GetStartTime();// - 9.92;
    fLength = simplifiedEvent.GetLength();

    bool longerThanPrompt = fLength > kPromptLength;

    fEnergy = simplifiedEvent.GetIntegralTank(longerThanPrompt);// - randomRate*kPromptLength;
    fHits = simplifiedEvent.GetNumTank(longerThanPrompt);
    fNumPMTs = simplifiedEvent.GetPMTHits();
    //MsgInfo(MsgLog::Form("e %zu energy %.2f hits %d numPTs %d length %.3f",e,fEnergy,fHits,fNumPMTs,fLength));

    if (fShiftTime) {
      st -= 9.92;
    }
    fTime = st;

    if (fDoTimeCut) {
      if (st < fTimeCutValueLow || st > fTimeCutValueHigh) {
        locationsToRemove.emplace_back(e);
        continue;
      } // end time cut condition
    } // end fDoTimeCut

    double et = st + fLength;
    if (et <= st) {
      locationsToRemove.emplace_back(e);
      continue;
    }

    if (fDoLengthCut) {
      if (fLength < fLengthCutValueLow || fLength > fLengthCutValueHigh) {
        locationsToRemove.emplace_back(e);
        continue;
      } // end length cut condition
    } // end if fDoLengthCut

    fX = simplifiedEvent.GetXPosition();
    fY = simplifiedEvent.GetYPosition();
    fZ = simplifiedEvent.GetZPosition();

    double radius = std::sqrt(fX*fX+fY*fY);
    if (fDoPositionCut) {
      if (radius < fRadiusCutValueLow || radius > fRadiusCutValueHigh ||
          fZ < fZCutValueLow || fZ > fZCutValueHigh) {
        locationsToRemove.emplace_back(e);
        continue;
      } // end position cut check
    } // end if fDoPositionCut

    if (fDoEnergyCut) {
      if (fEnergy < fEnergyCutValueLow || fEnergy > fEnergyCutValueHigh) {
        locationsToRemove.emplace_back(e);
        continue;
      }// end energy cut check
    }// end if fDoEnergyCut

    if (fHits < 3) {
        locationsToRemove.emplace_back(e);
      continue;
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
    int vetoTotal = simplifiedEvent.GetNumVeto(false);

    fLargestPMTFraction = simplifiedEvent.GetLargestPMTFraction();
    if (fDoNicenessCut) {
      double uniform = fEnergy/static_cast<double>(fNumPMTs);
      double nice = (fLargestPMTFraction - uniform)/uniform;
      if (nice < fNicenessCutValueLow || nice > fNicenessCutValueHigh) {
        locationsToRemove.emplace_back(e);
        continue;
      }
    } // end if fDoNicenessCut

    if (fDoVetoCut) {
      if ((vetoTotal >= fNumVetoCut && !fReverseVeto) || 
          (fReverseVeto && vetoTotal < fNumVetoCut)) {
        locationsToRemove.emplace_back(e);
        continue;
      } // end veto cut check
    }// end if fDoVetoCut


    bool reject = false;
    for (long e2 = e-1; e2 >= 0 && fDoPrevCut; --e2) {
      auto simplifiedEvent2 = fEvents->GetSimplifiedEvent(e2);

      bool longerThanPrompt2 = simplifiedEvent2.GetLength() > kPromptLength;
      double promptHits2 = simplifiedEvent2.GetNumCoated(longerThanPrompt2);
      double st2 = simplifiedEvent2.GetStartTime();// - 9.92;
      //double promptIntegral2 = simplifiedEvent2.GetIntegralTank(longerThanPrompt2);// - randomRate*kPromptLength;
      //double length2 = simplifiedEvent2.GetLength();

      if (fShiftTime) {
        st2 -= 9.92;
      }

      if (promptHits2 < 3) {
        continue;
      }

      double x2 = simplifiedEvent2.GetXPosition();
      double y2 = simplifiedEvent2.GetYPosition();
      double z2 = simplifiedEvent2.GetZPosition();
      if (std::sqrt(x2*x2+y2*y2) > 0.8 || std::fabs(z2) > 0.35) {
        continue;
      }

      double stDiff = st2 - st;

      //double pIDiff = promptIntegral - promptIntegral2;
      if (stDiff < fPrevCutTime) {
        reject = true;
        break;
      }
    } // end for e2

    if (reject) {
      locationsToRemove.emplace_back(e);
      continue;
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

  c("ShiftTime").Get(fShiftTime);

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

  if (!fOutFileName.empty()) {
    fOutfile = gROOT->GetFile(fOutFileName.c_str());
    if (!fOutfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fOutFileName.c_str()));
      return;
    }
    fOutfile->cd();
    fTree = new TTree("tree","tree");
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


