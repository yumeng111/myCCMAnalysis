/*!**********************************************
 * \file CCMNuDMPlots.cxx
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
#include "CCMNuDMPlots.h"
#include "CCMModuleTable.h"

#include "Events.h"
#include "SimplifiedEvent.h"
#include "MsgLog.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "CCMBeamInfo.h"

#include <iterator>
#include <sstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"

//See CCMModuleTable for info
MODULE_DECL(CCMNuDMPlots);

//_______________________________________________________________________________________
CCMNuDMPlots::CCMNuDMPlots(const char* version) 
  : CCMModule("CCMNuDMPlots")
{
  //Default constructor
  this->SetCfgVersion(version);

  fOutFileName = "";
  fInFileName = "";
  fPrevInFileName = "";
  fOutfile = nullptr;
  fInfile = nullptr;
  fTreeName = "tree";
  fIncludeBCM = false;
  fTotalPOT = 0;
  fTotalTriggers = 0;

  fTimeWindows.emplace_back(-260);
  fTimeWindows.emplace_back(-340);
  while (fTimeWindows.back() > -4520) {
    fTimeWindows.emplace_back(fTimeWindows.back()-190);
  }
  std::reverse(fTimeWindows.begin(),fTimeWindows.end());

  fNumTimeRegions = fTimeWindows.size();

  fHists1D.clear();
  fHists2D.clear();

  /*
  fCutNames.emplace_back("None");
  if (fIncludeBCM) {
    fCutNames.emplace_back("BCM");
  }
  fCutNames.emplace_back("Time");
  fCutNames.emplace_back("NumHits");
  fCutNames.emplace_back("Prev"); // 1p0
  //fCutNames.emplace_back("Prev"); // 0p75
  //fCutNames.emplace_back("Prev"); // 0p50
  //fCutNames.emplace_back("Prev"); // 0p25
  //fCutNames.emplace_back("Prev"); // 0p10
  //fCutNames.emplace_back("Prev"); // 0p0
  fCutNames.emplace_back("Veto");
  //fCutNames.emplace_back("Position");
  fCutNames.emplace_back("Niceness");
  fCutNames.emplace_back("WaveformMax");
  fCutNames.emplace_back("Length"); // 40
  fCutNames.emplace_back("Length"); // 50
  fCutNames.emplace_back("Length"); // 52
  fCutNames.emplace_back("Length"); // 60
  fCutNames.emplace_back("Length"); // 70
  fCutNames.emplace_back("Length"); // 80
  fCutNames.emplace_back("Length"); // 90
  fCutNames.emplace_back("Length"); // 100
  //fCutNames.emplace_back("Energy");
  */

}

//_______________________________________________________________________________________
CCMNuDMPlots::CCMNuDMPlots(const CCMNuDMPlots& clufdr) 
: CCMModule(clufdr)
{
  // copy constructor
  fOutFileName = clufdr.fOutFileName;
  fInFileName = clufdr.fInFileName;
  fPrevInFileName = clufdr.fPrevInFileName;
  fOutfile = clufdr.fOutfile;
  fInfile = clufdr.fInfile;
  fTreeName = clufdr.fTreeName;
  fIncludeBCM = clufdr.fIncludeBCM;
  fTotalPOT = clufdr.fTotalPOT;
  fTotalTriggers = clufdr.fTotalTriggers;

  fTimeWindows = clufdr.fTimeWindows;


  fNumTimeRegions = clufdr.fNumTimeRegions;
  fCutNames = clufdr.fCutNames;
  fHists1D = clufdr.fHists1D;
  fHists2D = clufdr.fHists2D;
}

//_______________________________________________________________________________________
CCMNuDMPlots::~CCMNuDMPlots()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMNuDMPlots::ProcessTrigger()
{
  // everything in this module is done when the ConnectInFileName is run
  // because it gets its information from a tree other than the EventTree
  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMNuDMPlots::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.

  MsgInfo("Inside Configuration file");

  c("TreeName").Get(fTreeName);
  MsgInfo(MsgLog::Form("- TreeName %s",fTreeName.c_str()));

  std::string tempString = "";
  c("CutNames").Get(tempString);
  MsgInfo(MsgLog::Form("- CutNames string %s",tempString.c_str()));
  std::stringstream ss(tempString);
  fCutNames.clear();
  fCutNames.assign((std::istream_iterator<std::string>(ss)),
      std::istream_iterator<std::string>());
  if (!fCutNames.empty()) {
    for (auto & name : fCutNames) {
      MsgInfo(MsgLog::Form("- - Cut %s",name.c_str()));
    }
    DefineHists();

    if (std::find(fCutNames.begin(),fCutNames.end(),"BCM") != fCutNames.end()) {
      fIncludeBCM = true;
    } else {
      fIncludeBCM = false;
    }
  } else {
    MsgFatal("To use this module you must past at least 1 cut name");
  }

  fIsInit = true;

}

/*!**********************************************
 * \fn void CCMNuDMPlots::SetupOutFile()
 * \brief Setup the output file that gets saved
 ***********************************************/
void CCMNuDMPlots::SetupOutFile()
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

/*!**********************************************
 * \fn void CCMNuDMPlots::SetupInputFile()
 * \brief Setup the input file that was read in to get the hitograms
 ***********************************************/
void CCMNuDMPlots::SetupInputFile()
{
  if (fInFileName == fPrevInFileName) {
    return;
  }

  if (!fInFileName.empty()) {
    fPrevInFileName = fInFileName;
    fInfile = gROOT->GetFile(fInFileName.c_str());
    if (!fInfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fInFileName.c_str()));
      return;
    }
    fInfile->cd();

    TTree * tree = nullptr;
    fInfile->GetObject(fTreeName.c_str(),tree);
    if (tree == nullptr) {
      MsgWarning(MsgLog::Form("Could not find Tree with name %s inside of %s",fTreeName.c_str(),fInFileName.c_str()));
      return;
    }

    MsgInfo(MsgLog::Form("Getting tree for file %s",fInFileName.c_str()));
    Get(tree);
    delete tree;
  }

}

//------------------------------------------------------
CCMResult_t CCMNuDMPlots::EndOfJob()
{
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  if (fOutfile != nullptr) {
    fOutfile->cd();
    for (auto & pair : fHists1D) {
      if (pair.second != nullptr) {
        pair.second->Write();
      }
    }
    for (auto & pair : fHists2D) {
      if (pair.second != nullptr) {
        pair.second->Write();
      }
    }
  }

  return kCCMSuccess;
}


//-------------------------------------------------------------------------------------------------
void CCMNuDMPlots::Get(TTree * tree)
{
  if (!tree) {
    return;
  }

  double energy = 0;
  double length = 0;
  double time = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double bcmIntegral = 0;
  unsigned int sec = 0;
  unsigned int ns = 0;
  double maxEnergy1600 = 0;
  double maxEnergy3200 = 0;
  double maxEnergy4800 = 0;
  double maxEnergy = 0;
  double maxEnergyTime1600 = 0;
  double maxEnergyTime3200 = 0;
  double maxEnergyTime4800 = 0;
  double maxEnergyTime = 0;
  double maxEnergyLength1600 = 0;
  double maxEnergyLength3200 = 0;
  double maxEnergyLength4800 = 0;
  double maxEnergyLength = 0;

  //double hits = 0;
  //int numPMTs = 0;
  //std::vector<int> vetoValues(fgkNumVetoRegions,0);
  //std::vector<int> vetoPromptValues(fgkNumVetoRegions,0);
  //double weight = 0;
  //double largestPMTFraction = 0;

  std::vector<int> passedCut(fCutNames.size()-1,false);

  if (fIncludeBCM) {
    tree->SetBranchAddress("bcmIntegral",&bcmIntegral);
  }
  tree->SetBranchAddress("epochSec",&sec);
  tree->SetBranchAddress("epochNS",&ns);
  tree->SetBranchAddress("energy",&energy);
  tree->SetBranchAddress("length",&length);
  tree->SetBranchAddress("time",&time);
  tree->SetBranchAddress("x",&x);
  tree->SetBranchAddress("y",&y);
  tree->SetBranchAddress("z",&z);
  tree->SetBranchAddress("maxPrevEnergy1600",&maxEnergy1600);
  tree->SetBranchAddress("maxPrevEnergy3200",&maxEnergy3200);
  tree->SetBranchAddress("maxPrevEnergy4800",&maxEnergy4800);
  tree->SetBranchAddress("maxPrevEnergy",&maxEnergy);
  tree->SetBranchAddress("maxPrevEnergyTime1600",&maxEnergyTime1600);
  tree->SetBranchAddress("maxPrevEnergyTime3200",&maxEnergyTime3200);
  tree->SetBranchAddress("maxPrevEnergyTime4800",&maxEnergyTime4800);
  tree->SetBranchAddress("maxPrevEnergyTime",&maxEnergyTime);
  tree->SetBranchAddress("maxPrevEnergyLength1600",&maxEnergyLength1600);
  tree->SetBranchAddress("maxPrevEnergyLength3200",&maxEnergyLength3200);
  tree->SetBranchAddress("maxPrevEnergyLength4800",&maxEnergyLength4800);
  tree->SetBranchAddress("maxPrevEnergyLength",&maxEnergyLength);
  for (size_t loc = 0; loc < fCutNames.size()-1; ++loc) {
    //std::vector<bool>::reference proxyForBit = passedCut.at(loc);
    //auto && proxyForBit = passedCut.at(loc);
    tree->SetBranchAddress(Form("passed%sCut",fCutNames.at(loc+1).c_str()),&passedCut.at(loc));
  }
  //tree->SetBranchAddress("hits",&hits);
  //tree->SetBranchAddress("numPMTs",&numPMTs);
  //tree->SetBranchAddress("vetoPromptTop",&vetoPromptValues.at(0));
  //tree->SetBranchAddress("vetoPromptBottom",&vetoPromptValues.at(1));
  //tree->SetBranchAddress("vetoPromptCLeft",&vetoPromptValues.at(2));
  //tree->SetBranchAddress("vetoPromptCRight",&vetoPromptValues.at(3));
  //tree->SetBranchAddress("vetoPromptCFront",&vetoPromptValues.at(4));
  //tree->SetBranchAddress("vetoPromptCBack",&vetoPromptValues.at(5));
  //tree->SetBranchAddress("vetoTop",&vetoValues.at(0));
  //tree->SetBranchAddress("vetoBottom",&vetoValues.at(1));
  //tree->SetBranchAddress("vetoCLeft",&vetoValues.at(2));
  //tree->SetBranchAddress("vetoCRight",&vetoValues.at(3));
  //tree->SetBranchAddress("vetoCFront",&vetoValues.at(4));
  //tree->SetBranchAddress("vetoCBack",&vetoValues.at(5));
  //tree->SetBranchAddress("weight",&weight);
  //tree->SetBranchAddress("largestPMTFraction",&largestPMTFraction);

  MsgInfo(MsgLog::Form("Size of passedCut %zu",passedCut.size()));
  std::vector<std::pair<unsigned int, unsigned int>> usedTimes;

  //const size_t kTimeLocation = std::distance(fCutNames.begin(),std::find(fCutNames.begin(),fCutNames.end(),"Time"));
  const size_t kBCMLocation = std::distance(fCutNames.begin(),std::find(fCutNames.begin(),fCutNames.end(),"BCM"));
  const size_t kPrevLocation = std::distance(fCutNames.begin(),std::find(fCutNames.begin(),fCutNames.end(),"Prev"));
  const size_t kLengthLocation = std::distance(fCutNames.begin(),std::find(fCutNames.begin(),fCutNames.end(),"Length"));
  //const size_t kNumPrev = std::count(fCutNames.begin(),fCutNames.end(),"Prev");
  const size_t kNumLength = std::count(fCutNames.begin(),fCutNames.end(),"Length");

  const long kEntries = tree->GetEntries();
  MsgInfo(MsgLog::Form("Number of events: %ld",kEntries));
  for (long e = 0; e < kEntries; ++e) {
    int nbytes = tree->GetEntry(e);
    if (nbytes == 0) {
      break;
    }

    if (e % 10000 == 0) {
      MsgInfo(MsgLog::Form("At %ld out of %ld",e+1,kEntries));
    }

    // want to convert from meters to centimeters
    x *= 100.0;
    y *= 100.0;
    z *= 100.0;
    double r = std::sqrt(x*x+y*y);

    size_t location = 0;
    if (time >= fTimeWindows.front()) {
      for (size_t region = 1; region < fNumTimeRegions-1; ++region) {
        if (time < fTimeWindows.at(region)) {
          location = region;
          break;
        }
      }
    }

    for (size_t loc = 0; loc < fCutNames.size(); ++loc) {

      if (loc >= kBCMLocation) {
        if (bcmIntegral < 17811.4 || bcmIntegral > 30978.4)
        //if (bcmIntegral < 24394.9-1316.8*3 || bcmIntegral > 24394.9+1316.8*3)
        {
          break;
        }
      }

      //if (loc >= kTimeLocation) {
      //  if (time > -1000) {
      //    break;
      //  }
      //}

      if (loc == kBCMLocation) {
        if (usedTimes.empty()) {
            usedTimes.emplace_back(std::make_pair(sec,ns));
            fTotalPOT += bcmIntegral*fgkBCMIntToPOT;
            ++fTotalTriggers;
        } else {
          auto p = usedTimes.back();
          if (p.first != sec && p.second != ns) {
            usedTimes.emplace_back(std::make_pair(sec,ns));
            fTotalPOT += bcmIntegral*fgkBCMIntToPOT;
            ++fTotalTriggers;
          }
        }
      }

      if (loc >= kPrevLocation) {
        if (maxEnergy != 0) {
          break;
        }
      }

      std::string extraName = "";
      if (loc >= kLengthLocation) {
        size_t number = loc - kLengthLocation;
        double threshold = 100;
        switch (number) {
          case 0: threshold = 200; extraName = "0"; break;
          case 1: threshold = 175; extraName = "1"; break;
          case 2: threshold = 150; extraName = "2"; break;
          case 3: threshold = 125; extraName = "3"; break;
          case 4: threshold = 100; extraName = "4"; break;
          case 5: threshold = 52; extraName = "5"; break;
          case 6: threshold = 50; extraName = "6"; break;
          case 7: threshold = 48; extraName = "7"; break;
          default: threshold = 52; extraName = "5"; break;
        }
        if (length > threshold) {
          break;
        }
      }
      if (loc != 0 && (loc < kLengthLocation || loc >= kLengthLocation+kNumLength)) {
        if (!passedCut.at(loc-1)) {
          break;
        }
          extraName = "";
      }

      Fill(Form("timeHist_%s%s",fCutNames.at(loc).c_str(),extraName.c_str()),time+1e-4,energy);
      Fill(Form("timeEnergyHist_%s%s",fCutNames.at(loc).c_str(),extraName.c_str()),time+1e-4,energy);
      Fill(Form("timeLengthHist_%s%s",fCutNames.at(loc).c_str(),extraName.c_str()),time+1e-4,length);
      
      if (time < -340) {

        Fill(Form("eventLengthDist_%s%s_%zu",fCutNames.at(loc).c_str(),extraName.c_str(),location),energy,length);
        Fill(Form("eventDist_%s%s_%zu",fCutNames.at(loc).c_str(),extraName.c_str(),location),energy);
        if (energy >= 3 && energy <= 1000 && time < -1000) {
          Fill(Form("positionXY_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),x,y);
          Fill(Form("positionRZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),r,z);
          Fill(Form("positionXZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),x,z);
          Fill(Form("positionYZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),y,z);
          Fill(Form("positionR2Z_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),std::pow(r/95.0,2.0),z);
        }
        if (energy >= 3 && energy <= 100 && time < -1000) {
          Fill(Form("positionXY_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),1),x,y);
          Fill(Form("positionRZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),1),r,z);
          Fill(Form("positionXZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),1),x,z);
          Fill(Form("positionYZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),1),y,z);
          Fill(Form("positionR2Z_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),1),std::pow(r/95.0,2.0),z);
        }
        if (energy >= 100 && energy <= 1000 && time < -1000) {
          Fill(Form("positionXY_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),2),x,y);
          Fill(Form("positionRZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),2),r,z);
          Fill(Form("positionXZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),2),x,z);
          Fill(Form("positionYZ_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),2),y,z);
          Fill(Form("positionR2Z_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),2),std::pow(r/95.0,2.0),z);
        }

        Fill(Form("EnergyMaxEnergy1600Hist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergyTime1600-time+1e-4,maxEnergy1600);
        Fill(Form("EnergyMaxEnergy3200Hist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergyTime3200-time+1e-4,maxEnergy3200);
        Fill(Form("EnergyMaxEnergy4800Hist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergyTime4800-time+1e-4,maxEnergy4800);
        Fill(Form("EnergyMaxEnergyHist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergyTime-time+1e-4,maxEnergy);

        Fill(Form("EnergyMaxLength1600Hist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergy1600,maxEnergyLength1600);
        Fill(Form("EnergyMaxLength3200Hist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergy3200,maxEnergyLength3200);
        Fill(Form("EnergyMaxLength4800Hist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergy4800,maxEnergyLength4800);
        Fill(Form("EnergyMaxLengthHist_%s%s_%d",fCutNames.at(loc).c_str(),extraName.c_str(),0),
            maxEnergy,maxEnergyLength);
      }

    } // end for loc < fCutNames.size()
  } // end for e
  MsgInfo(MsgLog::Form("Filled histograms from file %s now at POT %.0f (%.3g) and Triggers %.0f",
        fInFileName.c_str(),fTotalPOT,fTotalPOT,fTotalTriggers));
} // end of Get function

//-------------------------------------------------------------------------------------------------
void CCMNuDMPlots::Fill(std::string name, double par0, double par1, double par2)
{
  auto iter = fHists1D.find(name);
  if (iter == fHists1D.end()) {
    auto iter2 = fHists2D.find(name);
    if (iter2 == fHists2D.end()) {
      return;
    }

    if (iter2->second == nullptr) {
      return;
    }

    iter2->second->Fill(par0, par1, par2);

    return;
  }

  if (iter->second == nullptr) {
    return;
  }

  iter->second->Fill(par0,par1);

}// end CCMNuDMPlots::Fill

//-------------------------------------------------------------------------------------------------
void CCMNuDMPlots::DefineHists()
{

  std::vector<double> energyBins;
  energyBins.emplace_back(0.0);
  for (double value = 1; value < 10; value+=1.0) {
    energyBins.emplace_back(value);
  }
  for (double value = 10; value < 50; value+=2.0) {
    energyBins.emplace_back(value);
  }
  for (double value = 50; value < 100; value+=5) {
    energyBins.emplace_back(value);
  }
  for (double value = 100; value < 200; value+=10) {
    energyBins.emplace_back(value);
  }
  for (int exp = 0.0; exp < 5.0; ++exp) {
    for (int digit = 0; digit < 36; ++digit) {
      double value =(1.0+static_cast<double>(digit)/4.0)*std::pow(10.0,static_cast<double>(exp));
      if (value < 200) {
        continue;
      }
      energyBins.emplace_back(value);
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
  /*
  std::vector<double> energyBins;
  energyBins.push_back(0.0);
  //for (double  value = 1; value < 20; value+=0.25) {
  //  energyBins.push_back(value);
  //}
  //for (double  value = 20; value < 100; ++value) {
  //  energyBins.push_back(value);
  //}
  for (int exp = 0.0; exp < 5.0; ++exp) {
    for (int digit = 1; digit < 10; ++digit) {
      //double value =(1.0+static_cast<double>(digit)/4.0)*std::pow(10.0,static_cast<double>(exp));
      double value =static_cast<double>(digit)*std::pow(10.0,static_cast<double>(exp));
      //if (value < 100) {
      //  continue;
      //}
      energyBins.push_back(value);
      //if (value >= 500.0) {
      //  ++digit;
      //}
      //if (exp < 1) {
      //  digit += 0.9;
      //} else if (exp < 2) {
      //  digit += 0.15;
      //} else {
      //  digit += 0.9;
      //}
    }
  }
  */


  fHists1D.emplace("beamHist",
      std::make_shared<TH1D>("beamHist",
        ";Time (ns);arb.",
        Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime));
  fHists1D.emplace("piHist",
      std::make_shared<TH1D>("piHist",
        ";Time (ns);arb.",
        Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime));
  fHists1D.emplace("muonHist",
      std::make_shared<TH1D>("muonHist",
        ";Time (ns);arb.",
        Utility::fgkNumBins,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime));

  size_t prevNumber = 0;
  for (const auto & cutName : fCutNames) {
    std::string extraName = "";
    if (cutName == "Length") {
      extraName = std::to_string(prevNumber);
      ++prevNumber;
    }
    std::string name = "timeEnergyHist_"+cutName+extraName;
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";Time (ns);Energy (PE);Events/Day",
          Utility::fgkNumBins/5,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime,
          energyBins.size()-1,&energyBins.front()));

    name = "timeLengthHist_"+cutName+extraName;
    fHists2D.emplace(name.c_str(),std::make_shared<TH2D>(name.c_str(),
          ";Time (ns);Length (ns);Count",
          Utility::fgkNumBins/5,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime,
          Utility::fgkNumBins/5,0, std::fabs(Utility::fgkWindowEndTime-Utility::fgkWindowStartTime)));


    name = "EnergyMaxEnergy1600Hist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";Time Before (ns);1600 Max Energy (PE);Events/Day",
          4500,-9000,0,
          energyBins.size()-1,&energyBins.front()));

    name = "EnergyMaxEnergy3200Hist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";Time Before (ns);3200 Max Energy (PE);Events/Day",
          4500,-9000,0,
          energyBins.size()-1,&energyBins.front()));

    name = "EnergyMaxEnergy4800Hist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";Time Before (ns);4800 Max Energy (PE);Events/Day",
          4500,-9000,0,
          energyBins.size()-1,&energyBins.front()));

    name = "EnergyMaxEnergyHist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";Time Before (ns); Max Energy (PE);Events/Day",
          4500,-9000,0,
          energyBins.size()-1,&energyBins.front()));

    name = "EnergyMaxLength1600Hist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";1600 Max Energy (PE);Length (ns);Events/Day",
          energyBins.size()-1,&energyBins.front(),
            Utility::fgkNumBins/5,0, std::fabs(Utility::fgkWindowEndTime-Utility::fgkWindowStartTime)));

    name = "EnergyMaxLength3200Hist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";3200 Max Energy (PE);Length (ns);Events/Day",
          energyBins.size()-1,&energyBins.front(),
            Utility::fgkNumBins/5,0, std::fabs(Utility::fgkWindowEndTime-Utility::fgkWindowStartTime)));

    name = "EnergyMaxLength4800Hist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";4800 Max Energy (PE);Length (ns);Events/Day",
          energyBins.size()-1,&energyBins.front(),
            Utility::fgkNumBins/5,0, std::fabs(Utility::fgkWindowEndTime-Utility::fgkWindowStartTime)));

    name = "EnergyMaxLengthHist_"+cutName+extraName+"_0";
    fHists2D.emplace(name.c_str(),
        std::make_shared<TH2D>(name.c_str(),
          ";Max Energy (PE);Length (ns);Events/Day",
          energyBins.size()-1,&energyBins.front(),
            Utility::fgkNumBins/5,0, std::fabs(Utility::fgkWindowEndTime-Utility::fgkWindowStartTime)));

    name = "timeHist"+cutName+extraName;
    fHists1D.emplace(name.c_str(),
        std::make_shared<TH1D>(name.c_str(),
          ";Time (#mus);Accumulated Integral/Day",
          Utility::fgkNumBins/5,Utility::fgkWindowStartTime,Utility::fgkWindowEndTime));

    for (size_t i=0; i < fNumTimeRegions; ++i) {
      name = "eventDist_"+cutName+extraName+"_"+std::to_string(i);
      fHists1D.emplace(name.c_str(),
          std::make_shared<TH1D>(name.c_str(),
            ";Energy (PE);Events/PE/Day",
            energyBins.size()-1,&energyBins.front()));

      name = "positionXY_"+cutName+extraName+"_"+std::to_string(i);
      fHists2D.emplace(name.c_str(),
          std::make_shared<TH2D>(name.c_str(),
            ";X (cm);Y (cm);Events/Day",
            220,-120,120,220,-120,120));

      name = "positionRZ_"+cutName+extraName+"_"+std::to_string(i);
      fHists2D.emplace(name.c_str(),
          std::make_shared<TH2D>(name.c_str(),
            ";R (cm);Z (cm);Events/Day",
            220,-120,120,100,-50,50));

      name = "positionXZ_"+cutName+extraName+"_"+std::to_string(i);
      fHists2D.emplace(name.c_str(),
          std::make_shared<TH2D>(name.c_str(),
            ";R (cm);Z (cm);Events/Day",
            220,-120,120,100,-50,50));

      name = "positionYZ_"+cutName+extraName+"_"+std::to_string(i);
      fHists2D.emplace(name.c_str(),
          std::make_shared<TH2D>(name.c_str(),
            ";R (cm);Z (cm);Events/Day",
            220,-120,120,100,-50,50));

      name = "positionR2Z_"+cutName+extraName+"_"+std::to_string(i);
      fHists2D.emplace(name.c_str(),
          std::make_shared<TH2D>(name.c_str(),
            ";(R/95 cm)^{2};Z (cm);Events/Day",
            1440,0,std::pow(120.0/95.0,2.0),100,-50,50));

      name = "eventLengthDist_"+cutName+extraName+"_"+std::to_string(i);
      fHists2D.emplace(name.c_str(),
          std::make_shared<TH2D>(name.c_str(),
            ";Energy (PE);Length (ns);Events/PE/ns/Day",
            energyBins.size()-1,&energyBins.front(),
            Utility::fgkNumBins/5,0, std::fabs(Utility::fgkWindowEndTime-Utility::fgkWindowStartTime)));
    }
  } // end loop over fCutNames
}
