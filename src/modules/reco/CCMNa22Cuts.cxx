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
    fNicenessCutValueHigh(std::numeric_limits<double>::max())
{
  //Default constructor
  this->SetCfgVersion(version);
  if (!fInFileNames.empty()) {
    fInFileNames.clear();
  }
}

//_______________________________________________________________________________________
CCMNa22Cuts::CCMNa22Cuts(const CCMNa22Cuts& clufdr) 
: CCMModule(clufdr),
  fOutFileName(clufdr.fOutFileName),
  fInFileNames(clufdr.fInFileNames),
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
  fNicenessCutValueHigh(clufdr.fNicenessCutValueHigh)
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
  TH1::AddDirectory(0);
  TH1::SetDefaultSumw2(true);

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

  std::vector<double> energyBins;
  energyBins.push_back(0.0);
  for (double  value = 1; value < 20; value+=0.25) {
    energyBins.push_back(value);
  }
  for (double  value = 20; value < 400; ++value) {
    energyBins.push_back(value);
  }
  for (int exp = 0.0; exp < 5.0; ++exp) {
    for (int digit = 0; digit < 36; ++digit) {
      double value =(1.0+static_cast<double>(digit)/4.0)*std::pow(10.0,static_cast<double>(exp));
      if (value < 400) {
        continue;
      }
      energyBins.push_back(value);
      if (value >= 500.0) {
        ++digit;
      }
      //if (exp < 1) {
      //  digit += 0.9;
      //} else if (exp < 2) {
      //  digit += 0.15;
      //} else {
      //  digit += 0.9;
      //}
    }
  }

  double energy = 0;
  double length = 0;
  double hits = 0;
  int numPMTs = 0;
  double time = 0;
  int vetoTop = 0;
  int vetoBottom = 0;
  int vetoCT = 0;
  int vetoCB = 0;
  int vetoCFront = 0;
  int vetoCBack = 0;
  double weight = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double largestPMTFraction = 0;
  unsigned int epochSec = 0;

  auto outfile = TFile::Open(fOutFileName.c_str(),"RECREATE");
  auto tree = new TTree("tree","tree");
  tree->Branch("epochSec",&epochSec);
  tree->Branch("energy",&energy);
  tree->Branch("length",&length);
  tree->Branch("hits",&hits);
  tree->Branch("numPMTs",&numPMTs);
  tree->Branch("time",&time);
  tree->Branch("vetoTop",&vetoTop);
  tree->Branch("vetoBottom",&vetoBottom);
  tree->Branch("vetoCT",&vetoCT);
  tree->Branch("vetoCB",&vetoCB);
  tree->Branch("vetoCFront",&vetoCFront);
  tree->Branch("vetoCBack",&vetoCBack);
  tree->Branch("weight",&weight);
  tree->Branch("x",&x);
  tree->Branch("y",&y);
  tree->Branch("z",&z);
  tree->Branch("largestPMTFraction",&largestPMTFraction);

  unsigned int prevSec = 0;
  unsigned int prevNS = 0;
  //double prevEvent = 0;

  std::shared_ptr<TFile> file = nullptr;
  TChain * chain = new TChain("events","chain");
  for (const auto & fileName : fInFileNames) {
    chain->AddFile(fileName.c_str());
  }

  SimplifiedEvent* events = 0;

  chain->SetBranchAddress("events",&events);

  int nEntries = chain->GetEntries();
  MsgInfo(MsgLog::Form("Number of Events in Run is %d",nEntries));

  //double numUncoated = 22.0;
  //double numCoated = 77.0;

  long e = 0;
  for (e = 0; e < nEntries; ++e) {
    if (e % 100000 == 0) {
      MsgInfo(MsgLog::Form("%ld out of %ld (%.2f%%)",e,nEntries,static_cast<double>(e)/static_cast<double>(nEntries)*100.0));
    }
    chain->GetEntry(e);

    unsigned int seconds = events->GetComputerSecIntoEpoch();
    unsigned int nseconds = events->GetComputerNSIntoSec();

    if (seconds != prevSec && nseconds != prevNS) {
      prevSec = seconds;
      prevNS = nseconds;
    }

    epochSec = seconds;


    double st = events->GetStartTime();// - 9.92;
    length = events->GetLength();

    double promptLength = 0.090;
    bool longerThanPrompt = length > promptLength;
    double promptIntegral = events->GetIntegralTank(longerThanPrompt);// - randomRate*promptLength;
    double promptHits = events->GetNumCoated(longerThanPrompt);

    energy = promptIntegral;
    hits = promptHits;
    numPMTs = events->GetOtherEvents();

    //RecalculateStartTime(events,st,promptIntegral,promptHits,length);

    if (fShiftTime) {
      st -= 9.92;
    }
    time = st;

    if (fDoTimeCut) {
      if (st < fTimeCutValueLow || st > fTimeCutValueHigh) {
        continue;
      } // end time cut condition
    } // end fDoTimeCut

    //prevEvent = st;

    double et = st + length;
    if (et <= st) {
      continue;
    }

    if (fDoLengthCut) {
      if (length < fLengthCutValueLow || length > fLengthCutValueHigh) {
        continue;
      } // end length cut condition
    } // end if fDoLengthCut

    x = events->GetXPosition();
    y = events->GetYPosition();
    z = events->GetZPosition();

    double radius = std::sqrt(x*x+y*y);
    if (fDoPositionCut) {
      if (radius < fRadiusCutValueLow || radius > fRadiusCutValueHigh ||
          z < fZCutValueLow || z > fZCutValueHigh) {
        continue;
      } // end position cut check
    } // end if fDoPositionCut

    if (fDoEnergyCut) {
      if (promptIntegral < fEnergyCutValueLow || promptIntegral > fEnergyCutValueHigh) {
        continue;
      }// end energy cut check
    }// end if fDoEnergyCut

    if (promptHits < 3) {
      continue;
    }

    vetoTop = events->GetNumVetoTop();
    vetoBottom = events->GetNumVetoBottom();
    vetoCT = events->GetNumVetoFront();
    vetoCBack = events->GetNumVetoRight();
    vetoCFront = events->GetNumVetoLeft();
    vetoCB = events->GetNumVetoBack();
    int vetoTotal = events->GetNumVeto();

    largestPMTFraction = events->GetLargestPMTFraction();
    if (fDoNicenessCut) {
      double uniform = promptIntegral/static_cast<double>(numPMTs);
      double nice = (largestPMTFraction - uniform)/uniform;
      if (nice < fNicenessCutValueLow || nice > fNicenessCutValueHigh) {
        continue;
      }
    } // end if fDoNicenessCut

    if (fDoVetoCut) {
      if ((vetoTotal >= fNumVetoCut && !fReverseVeto) || 
          (fReverseVeto && vetoTotal < fNumVetoCut)) {
        continue;
      } // end veto cut check
    }// end if fDoVetoCut
    //MsgInfo(MsgLog::Form("Veto: T: %d B %d CT %d CB %d To: %d",
    //      vetoTop,vetoBottom,vetoCT,vetoCB,vetoTotal));


    bool reject = false;
    for (long e2 = e-1; e2 >= 0 && fDoPrevCut; --e2) {
      chain->GetEntry(e2);
      unsigned int sec2 = events->GetComputerSecIntoEpoch();
      unsigned int nsec2 = events->GetComputerNSIntoSec();

      if (sec2 != seconds || nsec2 != nseconds) {
        break;
      }

      bool longerThanPrompt2 = events->GetLength() > promptLength;
      double promptHits2 = events->GetNumCoated(longerThanPrompt2);
      double st2 = events->GetStartTime();// - 9.92;
      //double promptIntegral2 = events->GetIntegralTank(longerThanPrompt2);// - randomRate*promptLength;
      //double length2 = events->GetLength();

      //RecalculateStartTime(events,st2,promptIntegral2,promptHits2,length2);

      if (fShiftTime) {
        st2 -= 9.92;
      }

      if (promptHits2 < 3) {
        continue;
      }

      double x2 = events->GetXPosition();
      double y2 = events->GetYPosition();
      double z2 = events->GetZPosition();
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
      continue;
    }

    weight = 1.0;// /weightFunction.at(energyRegion)->Eval(st);
    //MsgInfo(MsgLog::Form("st %f weighFunction %g weight %g",st,weightFunction->Eval(st),weight));

    outfile->cd();
    tree->Fill();

  } // end e

  outfile->cd();
  tree->Write();
  delete outfile;

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMNa22Cuts::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Coonfiguration file");

  c("OutFileName").Get(fOutFileName);

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

  std::string temp = "";
  c("InFileName").Get(temp);

  MsgInfo(MsgLog::Form("File List Name ....%s....",temp.c_str()));

  std::ifstream infile(temp.c_str());
  MsgInfo(MsgLog::Form("infile.good(): %d",infile.good()));
  while (infile >> temp) {
    MsgInfo(MsgLog::Form("Adding file %s",temp.c_str()));
    fInFileNames.push_back(temp);
  }
  infile.close();
  MsgInfo(MsgLog::Form("Size of fInFileNames: %zu",fInFileNames.size()));

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMNa22Cuts::RecalculateStartTime(SimplifiedEvent * events, double & st, double & charge, double & hits, double & length) 
 * \brief Recalcualte the start time of an event given a different method
 * \param[in] events The event to look at
 * \param[out] st The new start time
 * \param[out] charge The new integral of the event
 * \param[out] hits The new number of hits in the event
 * \param[out] length The new length of the event
 ***********************************************/
void CCMNa22Cuts::RecalculateStartTime(SimplifiedEvent * events, double & st, double & charge, double & hits, double & length) 
{
  const std::vector<float> & waveform = events->GetWaveformInt();
  const std::vector<float> & waveformCount = events->GetWaveformCount();
  bool crossed = false;
  size_t stCount = 0;
  charge = 0;
  hits = 0;
  for (size_t count = 0; count < waveform.size(); ++count) {
    if (!crossed && waveform[count] >= 0.4) {
      crossed = true;
      st = events->GetStartTime()+static_cast<double>(count)*2e-3;
      stCount = count;
      charge = 0;
      hits = 0;
      length = (waveform.size() - count)*2e-3;
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

