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

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"


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
    fTimeCutValueHigh(std::numeric_limits<double>::max())
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
  fTimeCutValueHigh(clufdr.fTimeCutValueHigh)
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

  std::vector<double> energyBins2;
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

  std::vector<TH1D*> timeHist;
  std::vector<TH1D*> timeHistCount;
  TH1D* eventDist = new TH1D("eventDist",";Energy (PE);Count",energyBins.size()-1,&energyBins.front());
  TH1D* vetoCountTopHist = new TH1D("vetoCountTopHist",";Number Vetos;Count",30,0,30);
  TH1D* vetoCountBottomHist = dynamic_cast<TH1D*>(vetoCountTopHist->Clone("vetoCountBottomHist"));
  TH1D* vetoCountCTHist = dynamic_cast<TH1D*>(vetoCountTopHist->Clone("vetoCountCTHist"));
  TH1D* vetoCountCBHist = dynamic_cast<TH1D*>(vetoCountTopHist->Clone("vetoCountCBHist"));
  TH1D* vetoCountTotal = dynamic_cast<TH1D*>(vetoCountTopHist->Clone("vetoCountTotal"));

  TH2D* eventLengthDist = new TH2D("eventLengthDist",";Energy (PE);Length (#mus)",energyBins.size()-1,&energyBins.front(),8000,0,16.0);

  TH2D * positionXY = new TH2D("positionXY",";X (m);Y (m)",220,-1.2,1.2,220,-1.2,1.2);
  TH2D * positionXZ = new TH2D("positionXZ",";X (m);Z (m)",220,-1.2,1.2,80,-0.4,0.4);
  TH2D * positionYZ = new TH2D("positionYZ",";Y (m);Z (m)",220,-1.2,1.2,80,-0.4,0.4);

  auto vetoTopEnergy = new TH2D("vetoTopEnergy",";Energy (PE);Num Vetos;Count",energyBins.size()-1,&energyBins.front(),30,0,30);
  auto vetoBottomEnergy = dynamic_cast<TH2D*>(vetoTopEnergy->Clone("vetoBottomEnergy"));
  auto vetoCBEnergy = dynamic_cast<TH2D*>(vetoTopEnergy->Clone("vetoCBEnergy"));
  auto vetoCTEnergy = dynamic_cast<TH2D*>(vetoTopEnergy->Clone("vetoCTEnergy"));
  auto vetoTotalEnergy = dynamic_cast<TH2D*>(vetoTopEnergy->Clone("vetoTotalEnergy"));

  auto vetoTopLength = new TH2D("vetoTopLength",";Length (#mus);Num Vetos;Count",8000,0,16.0,30,0,30);
  auto vetoBottomLength = dynamic_cast<TH2D*>(vetoTopLength->Clone("vetoBottomLength"));
  auto vetoCBLength = dynamic_cast<TH2D*>(vetoTopLength->Clone("vetoCBLength"));
  auto vetoCTLength = dynamic_cast<TH2D*>(vetoTopLength->Clone("vetoCTLength"));
  auto vetoTotalLength = dynamic_cast<TH2D*>(vetoTopLength->Clone("vetoTotalLength"));

  for (int i=0; i < 4; ++i) {
    timeHist.push_back(new TH1D(Form("timeHist%d",i),";Time (#mus);Accumulated Integral",8000,-9.92,6.08));
    timeHistCount.push_back(new TH1D(Form("timeHistCount%d",i),";Time (#mus);Accumulated Number of Events",8000,-9.92,6.08));
  }

  unsigned int prevSec = 0;
  unsigned int prevNS = 0;
  //double prevEvent = 0;

  std::shared_ptr<TFile> file = nullptr;
  TChain * chain = new TChain("events","chain");
  for (const auto & fileName : fInFileNames) {
    chain->AddFile(fileName.c_str());
  }

  SimplifiedEvent* events = 0;

    //file = std::make_shared<TFile>(fileName.c_str(),"READ");

    //TTree * tree = 0;

    //file->GetObject("events",tree);

    //if (!tree) {
      //MsgError("rawData tree was not loaded");
    //}

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
        //prevEvent = 0.0;
      }


      double st = events->GetStartTime();// - 9.92;
      double length = events->GetLength()+1e-4;

      double promptLength = 0.090;
      bool longerThanPrompt = length-1e-4 > promptLength;
      double promptIntegral = events->GetIntegralTank(longerThanPrompt);// - randomRate*promptLength;
      double promptHits = events->GetNumCoated(longerThanPrompt);

      //RecalculateStartTime(events,st,promptIntegral,promptHits,length);

      if (fShiftTime) {
        st -= 9.92;
      }

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
        if (length-1e-4 < fLengthCutValueLow || length-1e-4 > fLengthCutValueHigh) {
          continue;
        } // end length cut condition
      } // end if fDoLengthCut

      double x = events->GetXPosition();
      double y = events->GetYPosition();
      double z = events->GetZPosition();

      if (fDoPositionCut) {
        double radius = std::sqrt(x*x+y*y);
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

      //if (promptIntegral < 20 || promptIntegral > 40) {
      //  continue;
      //}

      if (promptHits < 3) {
        continue;
      }

      int vetoTop = events->GetNumVetoTop();
      int vetoBottom = events->GetNumVetoBottom();
      int vetoCT = events->GetNumVetoFront();
      int vetoCB = events->GetNumVetoBack();
      int vetoTotal = events->GetNumVeto();

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
        //double length2 = events->GetLength()+1e-4;

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

      if (promptIntegral < 3) {
        timeHist[0]->Fill(st+1e-4,promptIntegral);
        timeHistCount[0]->Fill(st+1e-4);
      } else if (promptIntegral < 11) {
        timeHist[1]->Fill(st+1e-4,promptIntegral);
        timeHistCount[1]->Fill(st+1e-4);
      } else if (promptIntegral < 300) {
        timeHist[2]->Fill(st+1e-4,promptIntegral);
        timeHistCount[2]->Fill(st+1e-4);
      } else {
        timeHist[3]->Fill(st+1e-4,promptIntegral);
        timeHistCount[3]->Fill(st+1e-4);
      }

      eventDist->Fill(promptIntegral);
      eventLengthDist->Fill(promptIntegral,length);

      positionXY->Fill(x,y);
      positionXZ->Fill(x,z);
      positionYZ->Fill(y,z);

      vetoCountTopHist->Fill(vetoTop);
      vetoCountBottomHist->Fill(vetoBottom);
      vetoCountCTHist->Fill(vetoCT);
      vetoCountCBHist->Fill(vetoCB);
      vetoCountTotal->Fill(vetoTotal);

      vetoTopEnergy->Fill(promptIntegral,vetoTop);
      vetoBottomEnergy->Fill(promptIntegral,vetoBottom);
      vetoCBEnergy->Fill(promptIntegral,vetoCB);
      vetoCTEnergy->Fill(promptIntegral,vetoCT);
      vetoTotalEnergy->Fill(promptIntegral,vetoTotal);

      vetoTopLength->Fill(length,vetoTop);
      vetoBottomLength->Fill(length,vetoBottom);
      vetoCBLength->Fill(length,vetoCB);
      vetoCTLength->Fill(length,vetoCT);
      vetoTotalLength->Fill(length,vetoTotal);

      //MsgInfo(MsgLog::Form("Num Veto: Top %d, Bottom %d, CT %d CB %d",vetoTop,vetoBottom,vetoCT,vetoCB,vetoTotal));

    } // end e
  //} // end for fileName

  TFile* outfile = TFile::Open(fOutFileName.c_str(),"RECREATE");
  outfile->cd();
  for (auto & hist : timeHist) {
    hist->Write();
  }
  for (auto & hist : timeHistCount) {
    hist->Write();
  }
  eventDist->Write();
  eventLengthDist->Write();
  vetoCountTopHist->Write();
  vetoCountBottomHist->Write();
  vetoCountCTHist->Write();
  vetoCountCBHist->Write();
  vetoCountTotal->Write();
  vetoTopEnergy->Write();
  vetoBottomEnergy->Write();
  vetoCBEnergy->Write();
  vetoCTEnergy->Write();
  vetoTotalEnergy->Write();
  vetoTopLength->Write();
  vetoBottomLength->Write();
  vetoCBLength->Write();
  vetoCTLength->Write();
  vetoTotalLength->Write();
  positionXY->Write();
  positionXZ->Write();
  positionYZ->Write();
  delete outfile;

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMNa22Cuts::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Coonfiguration file");

  if (&c != 0)
  {
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
  } else {
    fIsInit = false;
  }

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

