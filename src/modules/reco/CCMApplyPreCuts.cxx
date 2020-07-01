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

#include "Events.h"
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
    fEvents()
    
{
  //Default constructor
  this->SetCfgVersion(version);

  SetupHists();
}

//_______________________________________________________________________________________
CCMApplyPreCuts::CCMApplyPreCuts(const CCMApplyPreCuts& clufdr) 
: CCMModule(clufdr),
  fOutFileName(clufdr.fOutFileName),
  fEvents(clufdr.fEvents)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMApplyPreCuts::~CCMApplyPreCuts()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMApplyPreCuts::ProcessTrigger()
{

  const size_t kNumEvents = fEvents->GetNumEvents();
  if (kNumEvents == 0) {
    MsgWarning("Number of Events in the fEvents object is equal to zero, doing nothing");
    return kCCMSuccess;
  }

  for (size_t e = 0; e < kNumEvents; ++e) {
    auto simplifiedEvent = fEvents->GetSimplifiedEvent(e);

    double st = simplifiedEvent.GetStartTime();// - 9.92;
    double length = simplifiedEvent.GetLength()+1e-4;

    double promptLength = 0.090;
    bool longerThanPrompt = length-1e-4 > promptLength;
    double promptIntegral = simplifiedEvent.GetIntegralTank(longerThanPrompt);// - randomRate*promptLength;
    double promptHits = simplifiedEvent.GetNumCoated(longerThanPrompt);

    std::cout << st << std::endl;

    RecalculateStartTime(simplifiedEvent,st,promptIntegral,promptHits,length);

    st -= 9.92;

    //prevEvent = st;

    double et = st + length;
    if (et <= st) {// || length > 0.1) {
      continue;
    }



    double x = simplifiedEvent.GetXPosition();
    double y = simplifiedEvent.GetYPosition();
    double z = simplifiedEvent.GetZPosition();

    if (std::sqrt(x*x+y*y) > 0.85 || std::fabs(z) > 0.3) {
      continue;
    }

    
    if (promptHits < 3) {
      continue;
    }

    bool reject = false;
    for (long e2 = e-1; e2 >= 0 ; --e2) {
      auto simplifiedEvent2 = fEvents->GetSimplifiedEvent(e2);
      
      bool longerThanPrompt2 = simplifiedEvent2.GetLength() > promptLength;
      double promptHits2 = simplifiedEvent2.GetNumCoated(longerThanPrompt2);
      double st2 = simplifiedEvent2.GetStartTime();// - 9.92;
      double promptIntegral2 = simplifiedEvent2.GetIntegralTank(longerThanPrompt2);// - randomRate*promptLength;
      double length2 = simplifiedEvent2.GetLength()+1e-4;

      RecalculateStartTime(simplifiedEvent2,st2,promptIntegral2,promptHits2,length2);

      st2 -= 9.92;

      if (promptHits2 < 3) {
        continue;
      }

      //if (st2 > -1) {
      //  continue;
      //}

      double x2 = simplifiedEvent2.GetXPosition();
      double y2 = simplifiedEvent2.GetYPosition();
      double z2 = simplifiedEvent2.GetZPosition();
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
      fTimeHist[0]->Fill(st+1e-4,promptIntegral);
    } else if (promptIntegral < 11) {
      fTimeHist[1]->Fill(st+1e-4,promptIntegral);
    } else if (promptIntegral < 300) {
      fTimeHist[2]->Fill(st+1e-4,promptIntegral);
    }

    if (st < -0.8 && st >= -1) {
      fPreEventHist->Fill(promptIntegral);
    } else if (st >= -0.7 && st < -0.65) {
      fBeamWindow1Hist->Fill(promptIntegral);
    } else if (st >= -0.65 && st < -0.6) {
      fBeamWindow2Hist->Fill(promptIntegral);
    } else if (st >= -0.6 && st < -0.55) {
      fBeamWindow3Hist->Fill(promptIntegral);
    } else if (st >= -0.55 && st < -0.5) {
      fBeamWindow4Hist->Fill(promptIntegral);
    } else if (st >= -0.5 && st < -0.45) {
      fBeamWindow5Hist->Fill(promptIntegral);
    } else if (st >= -0.45 && st < -0.4) {
      fBeamWindow6Hist->Fill(promptIntegral);
    } else if (st >= -0.4 && st < -0.35) {
      fBeamWindow7Hist->Fill(promptIntegral);
    } else if (st >= -0.35 && st < -0.3) {
      fBeamWindow8Hist->Fill(promptIntegral);
    }

  } // end e

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMApplyPreCuts::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  c("OutFileName").Get(fOutFileName);

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMApplyPreCuts::RecalculateStartTime(const SimplifiedEvent & events, double & st, double & charge, double & hits, double & length) 
 * \brief Recalcualte the start time of an event given a different method
 * \param[in] events The event to look at
 * \param[out] st The new start time
 * \param[out] charge The new integral of the event
 * \param[out] hits The new number of hits in the event
 * \param[out] length The new length of the event
 ***********************************************/
void CCMApplyPreCuts::RecalculateStartTime(const SimplifiedEvent & events, double & st, double & charge, double & hits, double & length) 
{
  const std::vector<float> & waveform = events.GetWaveformInt();
  const std::vector<float> & waveformCount = events.GetWaveformCount();
  bool crossed = false;
  size_t stCount = 0;
  charge = 0;
  hits = 0;
  for (size_t count = 0; count < waveform.size(); ++count) {
    if (!crossed && waveform[count] >= 0.4) {
      crossed = true;
      st = events.GetStartTime()+static_cast<double>(count)*2e-3;
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

/*!**********************************************
 * \fn void CCMApplyPreCuts::SetupHists()
 * \brief Setup the histograms to save the information
 ***********************************************/
void CCMApplyPreCuts::SetupHists()
{
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

  fPreEventHist = new TH1D("PreEventHist",";Energy (PE);Count",energyBins.size()-1,&energyBins.front());
  fBeamWindow1Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow1Hist"));
  fBeamWindow2Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow2Hist"));
  fBeamWindow3Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow3Hist"));
  fBeamWindow4Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow4Hist"));
  fBeamWindow5Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow5Hist"));
  fBeamWindow6Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow6Hist"));
  fBeamWindow7Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow7Hist"));
  fBeamWindow8Hist = dynamic_cast<TH1D*>(fPreEventHist->Clone("BeamWindow8Hist"));

  for (int i=0; i < 3; ++i) {
    fTimeHist.push_back(new TH1D(Form("TimeHist%d",i),";Time (#mus);Count",8000,-9.92,6.08));
  }
} // end CCMApplyPreCuts::SetupHists


CCMResult_t CCMApplyPreCuts::EndOfJob()
{ 
  //fPreEventHist->Scale(0.05/5.3);
  fPreEventHist->Scale(50.0/200.0);

  TFile* outfile = TFile::Open(fOutFileName.c_str(),"RECREATE");
  outfile->cd();
  for (auto & hist : fTimeHist) {
    hist->Write();
  }
  fPreEventHist->Write();
  fBeamWindow1Hist->Write();
  fBeamWindow2Hist->Write();
  fBeamWindow3Hist->Write();
  fBeamWindow4Hist->Write();
  fBeamWindow5Hist->Write();
  fBeamWindow6Hist->Write();
  fBeamWindow7Hist->Write();
  fBeamWindow8Hist->Write();
  delete outfile;

  return kCCMSuccess;
}

