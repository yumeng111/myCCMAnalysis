/*!**********************************************
 * \file CCMApplyPreCuts.cxx
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
#include "CCMApplyPreCuts.h"
#include "CCMModuleTable.h"

#include "SimplifiedEvent.h"
#include "MsgLog.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"


//See CCMModuleTable for info
MODULE_DECL(CCMApplyPreCuts);

//_______________________________________________________________________________________
CCMApplyPreCuts::CCMApplyPreCuts(const char* version) 
  : CCMModule("CCMApplyPreCuts"),
    fOutFileName(""),
    fInFileName("")
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMApplyPreCuts::CCMApplyPreCuts(const CCMApplyPreCuts& clufdr) 
: CCMModule(clufdr),
  fOutFileName(clufdr.fOutFileName),
  fInFileName(clufdr.fInFileName)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMApplyPreCuts::~CCMApplyPreCuts()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMApplyPreCuts::ProcessEvent()
{
  TH1::AddDirectory(0);

  TFile * file = TFile::Open(fInFileName.c_str(),"READ");

  TTree * tree = 0;

  file->GetObject("events",tree);

  if (!tree) {
    MsgError("rawData tree was not loaded");
  }

  SimplifiedEvent* events = 0;

  tree->SetBranchAddress("events",&events);

  int nEntries = tree->GetEntries();
  MsgInfo(MsgLog::Form("Number of Events in Run is %d",nEntries));

  //double numUncoated = 22.0;
  //double numCoated = 77.0;

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
  TH1D* preEventHist = new TH1D("preEventHist",";Energy (PE);Count",energyBins.size()-1,&energyBins.front());
  TH1D* beamWindow1Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow1Hist"));
  TH1D* beamWindow2Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow2Hist"));
  TH1D* beamWindow3Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow3Hist"));
  TH1D* beamWindow4Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow4Hist"));
  TH1D* beamWindow5Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow5Hist"));
  TH1D* beamWindow6Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow6Hist"));
  TH1D* beamWindow7Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow7Hist"));
  TH1D* beamWindow8Hist = dynamic_cast<TH1D*>(preEventHist->Clone("beamWindow8Hist"));

  for (int i=0; i < 3; ++i) {
    timeHist.push_back(new TH1D(Form("timeHist%d",i),";Time (#mus);Count",8000,-9.92,6.08));
  }

  unsigned int prevSec = 0;
  unsigned int prevNS = 0;
  //double prevEvent = 0;

  long e = 0;
  for (e = 0; e < nEntries; ++e) {
    if (e % 100000 == 0) {
      MsgInfo(MsgLog::Form("%ld out of %ld",e,nEntries));
    }
    tree->GetEntry(e);

    unsigned int seconds = events->GetComputerSecIntoEpoch();
    unsigned int nseconds = events->GetComputerNSIntoSec();

    if (seconds != prevSec && nseconds != prevNS) {

      //std::fill(integral.begin(),integral.end(),0.0);

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

    std::cout << st << std::endl;

    RecalculateStartTime(events,st,promptIntegral,promptHits,length);

    st -= 9.92;

    //prevEvent = st;

    double et = st + length;
    if (et <= st) {// || length > 0.1) {
      continue;
    }



    //double totalIntegral = events->GetIntegralTank(false);// - randomRate*(length-1e-4);

    //double promptIntegralC = events->GetIntegralCoated(true);// - randomPromptC;
    //double promptIntegralU = events->GetIntegralUncoated(true);// - randomPromptU;

    //double totalIntegralC = events->GetIntegralCoated(false);// - randomTotalC;
    //double totalIntegralU = events->GetIntegralUncoated(false);// - randomTotalU;

    /*
    const std::vector<float> & waveform = events->GetWaveformInt();
    double totalWeight = 0.0;
    double mean = 0;//std::accumulate(waveform.begin(),waveform.end(),0.f);
    double sum = 0;
    double count = 0;
    double peakLoc = 0;
    double peakValue = 0;
    for (const auto & point : waveform) {
      totalWeight += point*point;
      sum += point*point*count*count;
      mean += point*point*count;
      if (point > peakValue) {
        peakLoc = count;
        peakValue = point;
      }
      ++count;
    }
    mean  /= totalWeight;
    sum = std::sqrt(sum/totalWeight - mean*mean);
    */

    double x = events->GetXPosition();
    double y = events->GetYPosition();
    double z = events->GetZPosition();

    if (std::sqrt(x*x+y*y) > 0.85 || std::fabs(z) > 0.3) {
      continue;
    }

    
    if (promptHits < 3) {
      continue;
    }

    bool reject = false;
    for (long e2 = e-1; e2 >= 0 ; --e2) {
      tree->GetEntry(e2);
      unsigned int sec2 = events->GetComputerSecIntoEpoch();
      unsigned int nsec2 = events->GetComputerNSIntoSec();

      if (sec2 != seconds || nsec2 != nseconds) {
        break;
      }

      bool longerThanPrompt2 = events->GetLength() > promptLength;
      double promptHits2 = events->GetNumCoated(longerThanPrompt2);
      double st2 = events->GetStartTime();// - 9.92;
      double promptIntegral2 = events->GetIntegralTank(longerThanPrompt2);// - randomRate*promptLength;
      double length2 = events->GetLength()+1e-4;

      RecalculateStartTime(events,st2,promptIntegral2,promptHits2,length2);

      st2 -= 9.92;

      if (promptHits2 < 3) {
        continue;
      }

      //if (st2 > -1) {
      //  continue;
      //}

      double x2 = events->GetXPosition();
      double y2 = events->GetYPosition();
      double z2 = events->GetZPosition();
      if (std::sqrt(x2*x2+y2*y2) > 0.85 || std::fabs(z2) > 0.3) {
        continue;
      }

      //double stDiff = st2 - st;

      double pIDiff = promptIntegral - promptIntegral2;
      if (pIDiff < 1.6) {
        reject = true;
        break;
      }
    }

    if (reject) {
      continue;
    }

    if (promptIntegral < 3) {
      timeHist[0]->Fill(st+1e-4,promptIntegral);
    } else if (promptIntegral < 11) {
      timeHist[1]->Fill(st+1e-4,promptIntegral);
    } else if (promptIntegral < 300) {
      timeHist[2]->Fill(st+1e-4,promptIntegral);
    }

    if (st < -0.8 && st >= -1) {
      preEventHist->Fill(promptIntegral);
    } else if (st >= -0.7 && st < -0.65) {
      beamWindow1Hist->Fill(promptIntegral);
    } else if (st >= -0.65 && st < -0.6) {
      beamWindow2Hist->Fill(promptIntegral);
    } else if (st >= -0.6 && st < -0.55) {
      beamWindow3Hist->Fill(promptIntegral);
    } else if (st >= -0.55 && st < -0.5) {
      beamWindow4Hist->Fill(promptIntegral);
    } else if (st >= -0.5 && st < -0.45) {
      beamWindow5Hist->Fill(promptIntegral);
    } else if (st >= -0.45 && st < -0.4) {
      beamWindow6Hist->Fill(promptIntegral);
    } else if (st >= -0.4 && st < -0.35) {
      beamWindow7Hist->Fill(promptIntegral);
    } else if (st >= -0.35 && st < -0.3) {
      beamWindow8Hist->Fill(promptIntegral);
    }

  } // end e

  if (file) {
    delete file;
  }

  //preEventHist->Scale(0.05/5.3);
  preEventHist->Scale(50.0/200.0);

  TFile* outfile = TFile::Open(fOutFileName.c_str(),"RECREATE");
  outfile->cd();
  for (auto & hist : timeHist) {
    hist->Write();
  }
  preEventHist->Write();
  beamWindow1Hist->Write();
  beamWindow2Hist->Write();
  beamWindow3Hist->Write();
  beamWindow4Hist->Write();
  beamWindow5Hist->Write();
  beamWindow6Hist->Write();
  beamWindow7Hist->Write();
  beamWindow8Hist->Write();
  delete outfile;

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMApplyPreCuts::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("OutFileName").Get(fOutFileName);
  c("InFileName").Get(fInFileName);

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMApplyPreCuts::RecalculateStartTime(SimplifiedEvent * events, double & st, double & charge, double & hits, double & length) 
 * \brief Recalcualte the start time of an event given a different method
 * \param[in] events The event to look at
 * \param[out] st The new start time
 * \param[out] charge The new integral of the event
 * \param[out] hits The new number of hits in the event
 * \param[out] length The new length of the event
 ***********************************************/
void CCMApplyPreCuts::RecalculateStartTime(SimplifiedEvent * events, double & st, double & charge, double & hits, double & length) 
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

