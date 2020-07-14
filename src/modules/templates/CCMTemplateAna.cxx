/*!**********************************************
 * \file CCMTemplateAna.cxx
 * \author H. Ray (UF)
 * \date June 30, 2020
 *
 * Analysis Template for users to make a module
 *
 * !!!! VIP !!!!
 * Users must also edit the CMakeLists.txt and LinkDef.h files
 * in the module subdir where your module resides (eg reco)
 ***********************************************/
// include your .h file here
#include "CCMTemplateAna.h"

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMModuleTable.h"

#include "Events.h"
#include "SimplifiedEvent.h"
#include "MsgLog.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"

#include "RawData.h"
#include "Pulses.h"
#include "Utility.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <limits>
#include <memory>
#include <array>
#include <map>
#include <locale>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TVector3.h"


//See CCMModuleTable for info
MODULE_DECL(CCMTemplateAna);

// global var used in code applied for every entry in the input file go here

//_______________________________________________________________________________________
// *********************************
// CAUTION CAUTION CAUTION CAUTION
// vars in the constructor & destructor must be in the SAME ORDER as in the header file!
// *********************************
CCMTemplateAna::CCMTemplateAna(const char* version) 
  : CCMModule("CCMTemplateAna"),
    fEvents(nullptr),
    fRawData(nullptr),
    fPulses(nullptr),
    fOutFileName(""),
    fTriggerType("BEAM"),
    fNumTriggers(0),

    
    // module-specific vars in .xml file
    fMyDouble(0.4),
    //   fThreshold(0.4),
    fMyString("MYSTRING"),
    //   fTriggerType("BEAM"),
    fMyBool(false),
    //   fPulseTimeCut(false)
    fMyLowValue(0.1),
    fMyHighValue(0.9),
    fTreeName("tree"),

    
    // output file & tree
    fOutfile(nullptr),
    fTree(nullptr),

    
    // module-specific vars for output ntuple tree
    fmyOutDouble(0)
    // fEnergy(0)
    
{
  //Default constructor
  this->SetCfgVersion(version);

}


//_______________________________________________________________________________________
CCMTemplateAna::CCMTemplateAna(const CCMTemplateAna& clufdr) 
  : CCMModule(clufdr),
    fEvents(clufdr.fEvents),
    fRawData(clufdr.fRawData),
    fPulses(clufdr.fPulses),
    fOutFileName(clufdr.fOutFileName),
    fTriggerType(clufdr.fTriggerType),
    fNumTriggers(clufdr.fNumTriggers),


    // module-specific vars in .xml file
    fMyDouble(clufdr.fMyDouble),
    //   fThreshold(clufdr.fThreshold),
    fMyString(clufdr.fMyString),
    //   fTriggerType(clufdr.fTriggerType),
    fMyBool(clufdr.fMyBool),
    //   fPulseTimeCut(clufdr.fPulseTimeCut),
    fMyLowValue(clufdr.fMyLowValue),
    fMyHighValue(clufdr.fMyHighValue),
    fTreeName(clufdr.fTreeName),

    
    // output file & tree
    fOutfile(clufdr.fOutfile),
    fTree(clufdr.fTree),
  

    // module-specific vars for output ntuple tree
    fmyOutDouble(clufdr.fmyOutDouble)
    // fEnergy(clufdr.fEnergy)

{
  // copy constructor
}


//_______________________________________________________________________________________
CCMTemplateAna::~CCMTemplateAna()
{ 
  // destructor
  fMyHists.clear();
  //   fTimeHists.clear();
}


//_______________________________________________________________________________________
void CCMTemplateAna::Configure(const CCMConfig& c ) 
{
  //Initialize any parameters here
  //by reading them from the CCMConfig object.
  MsgInfo("Inside Coonfiguration file");

  c("TriggerType").Get(fTriggerType);


  c("MyDouble").Get(fMyDouble);
  //   c("Threshold").Get(fThreshold);
  c("MyString").Get(fMyString);
  //   c("TriggerType").Get(fTriggerType);
  c("MyBool").Get(fMyBool);
  //   c("PulseTimeCut").Get(fPulseTimeCut);


  if (fMyBool) {
    c("MyLowValue").Get(fMyLowValue);
    c("MyHighValue").Get(fMyHighValue);
    
    fMyLowValue = std::max(fMyLowValue,Utility::fgkWindowStartTime);
    fMyHighValue = std::min(fMyHighValue,Utility::fgkWindowEndTime);
    
    //   c("PulseTimeLowValue").Get(fPulseTimeLowValue);
    //   c("PulseTimeHighValue").Get(fPulseTimeHighValue);

    //   fPulseTimeLowValue = std::max(fPulseTimeLowValue,Utility::fgkWindowStartTime);
    //   fPulseTimeHighValue = std::min(fPulseTimeHighValue,Utility::fgkWindowEndTime);

  }

  
  MsgInfo("Input parameter values");
  MsgInfo(MsgLog::Form("-TriggerType: %s",fTriggerType.c_str()));

  MsgInfo(MsgLog::Form("-MyDouble: %f",fMyDouble));
  //   MsgInfo(MsgLog::Form("-Threshold: %f",fThreshold));
  MsgInfo(MsgLog::Form("-MyString: %s",fMyString.c_str()));
  //   MsgInfo(MsgLog::Form("-TriggerType: %s",fTriggerType.c_str()));
  MsgInfo(MsgLog::Form("-MyBool: %d",fMyBool));
  //   MsgInfo(MsgLog::Form("-PulseTimeCut: %d",fPulseTimeCut));


  MsgInfo(MsgLog::Form("-MyLowValue: %f",fMyLowValue));
  //   MsgInfo(MsgLog::Form("-PulseTimeLowValue: %f",fPulseTimeLowValue));
  MsgInfo(MsgLog::Form("-MyHighValue: %f",fMyHighValue));
  //   MsgInfo(MsgLog::Form("-PulseTimeHighValue: %f",fPulseTimeHighValue));

  
  
  // did you do the initialization?
  fIsInit = true;

} // Configure


//_______________________________________________________________________________________
CCMResult_t CCMTemplateAna::ProcessTrigger()
{
  // this is where you get the event found by FindEvent module
  // and apply your analysis code

  // since the class does not own the TFile, always check to see if it is present
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  
  // Check which trigger occured in DAQ window
  // make sure it is the one we want based on input .xml file
  // should have already been done in FindEvent; kept in Template
  // as a cross-check
  if (!fRawData->IsTriggerPresent(fTriggerType)) {
    return kCCMDoNotWrite;
  }

  // Reset vectors and histograms used for each event
  ResetVectorsAndHistos();

  
  // global counter for number of triggers processed
  ++fNumTriggers;

  
  // This function is code that is run over every event
  // It could go in ProcessEvent but was separated out
  // so that ProcessEvent isn't one huge unweildy function
  TemplateAna();

 
  // at the end of every event fill the ntuple
  if (fOutfile != nullptr) {
    fOutfile->cd();
    fTree->Fill();
  }

  return kCCMSuccess;

} // ProcessTrigger


//_______________________________________________________________________________________
void CCMTemplateAna::SetupOutFile()
{
  // setup the output file that gets saved
  if (fOutfile) {
    if (fOutfile->GetName() != fOutFileName) {
      MsgWarning(MsgLog::Form("New outfile name %s is different than old %s, this should not happen, sticking with ld",
			      fOutFileName.c_str(),fOutfile->GetName()));
    }
    return;
  }

  // ntuple tree vars
  fmyOutDouble = 0;
  // fEnergy = 0;

  
  if (!fOutFileName.empty()) {
    fOutfile = gROOT->GetFile(fOutFileName.c_str());
    if (!fOutfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fOutFileName.c_str()));
      return;
    }
    fOutfile->cd();
    fTree = new TTree("tree","tree");

    fTree->Branch("myOutDouble",&fmyOutDouble);
    // fTree->Branch("energy",&fEnergy);
  
  } // if no output file name

} //SetupOutputFile


//_______________________________________________________________________________________
CCMResult_t CCMTemplateAna::EndOfJob() 
{
  // print out any global counters
  MsgInfo(MsgLog::Form("Num Triggers %ld passed trigger type cut",fNumTriggers));
  
  
  // write out the ntuple & histograms to output file
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  if (fOutfile != nullptr) {
    fOutfile->cd();
    fTree->Write();

    for (auto & hist : fMyHists) {
      //   for (auto & hist : fTimeHists) {
        hist->Write();
    } // hists

  } // outfile
    
  return kCCMSuccess;

} // EndOfJob


//_______________________________________________________________________________________
void CCMTemplateAna::DefineVectorsAndHistos()
{
  fMyVector.resize(Utility::fgkNumBins,0.f);
  //   fPulsesTime.resize(Utility::fgkNumBins,0.f);


  fMySingleHist = new TH1D("MySingleHist",";Energy (PE);Count",8000, -9.92, 6.08);
  //   fPreEventHist = new TH1D("PreEventHist",";Energy (PE);Count",energyBins.size()-1,&energyBins.front());

  
  for (int i=0; i < 3; ++i) {
    fMyHists.push_back(new TH1D(Form("MyHists%d",i),";Time (#mus);Count", 8000, -9.92, 6.08));
    //   fTimeHist.push_back(new TH1D(Form("TimeHist%d",i),";Time (#mus);Count",8000,-9.92,6.08));
  }


} // DefineVectorsAndHistos


//_______________________________________________________________________________________
void CCMTemplateAna::ResetVectorsAndHistos()
{
  auto itMyVectorBegin = fMyVector.begin();
  //   auto itPulsesTimeBegin = fPulsesTime.begin();
  auto itMyVectorEnd = fMyVector.end();
  //   auto itPulsesTimeEnd = fPulsesTime.end();
  
  std::fill(itMyVectorBegin,itMyVectorEnd,0);
  //   std::fill(itPulsesTimeBegin,itPulsesTimeEnd,0);


  // no histogram in this template to be reset for each trigger processed
  
} // ResetVectorsAndHistos


//_______________________________________________________________________________________
void CCMTemplateAna::TemplateAna()
{
  // *************************************************************
  // Get the number of events in this trigger for your ana loop
  // *************************************************************
  const size_t kNumEvents = fEvents -> GetNumEvents();

  // loop over all events in this trigger
  for (size_t i = 0; i < kNumEvents; i++) {
    // get the event & all information linked to it
    auto simplifiedEvent = fEvents->GetSimplifiedEvent(i);

    
    // get var from event we want to write out to ntuple
    fmyOutDouble = simplifiedEvent.GetIntegralTank(true);
    // fEnergy = simplifiedEvent.GetIntegralTank(true);
    
    

    // TYLER: WHAT DO I DO ABOUT INCLUDING THESE IN THE TEMPLATE??
    // vector example
    // FindEvent never uses TimeBegin & TimeEnd ??
    
    // auto itMyVectorBegin = fMyVector.begin();
    // auto itPulsesTimeBegin = fPulsesTime.begin();
    
    // fPulsesTime.at(bin) += hitFraction;
    
    
    // *************************************************************
    // Declare variables used locally in analysis
    // *************************************************************
    int key = 0;
    
    double time = 0;
    double pulseIntegral = 0;
    
    
    // *************************************************************
    // Example loop over all pulses in the event
    // *************************************************************
    // Required
    auto pmtInfo = PMTInfoMap::GetPMTInfo(0);
    
    const size_t kNumPulses = fPulses->GetNumPulses();
    for (size_t loc = 0; loc < kNumPulses; ++loc) {
      
      // Required: make sure the PMT is to be used for the analysis
      key = fPulses->GetKey(loc);
      if (!PMTInfoMap::IsActive(key)) {
	continue;
      }

      // Required: make sure the PMT info is present
      pmtInfo = PMTInfoMap::GetPMTInfo(key);
      if (!pmtInfo) {
	continue;
      }
      
      
      // check to see if there is a cut on which pulses to include (fMyBool, from xml file)
      // if so, apply cut based on values in .xml file
      time = fPulses->GetPulseTime(loc);
      
      if (fMyBool) {
	if (time < fMyLowValue || time > fMyHighValue) continue;
      } // end if bool cut
      
      
	// get the integral of the pulse and make sure it is a real number and above threshold
	// threshold (fMyDouble, set in xml file)
      pulseIntegral = fPulses->GetPulseIntegral(loc);
      
      if (std::isnan(pulseIntegral) || std::isinf(pulseIntegral) || pulseIntegral < fMyDouble) {
	continue;
      }

      
      // example of using member functions from PMTInfoMap
      if (pmtInfo->IsVeto()) continue;


      // fill histogram
      fMyHists.at(key)->Fill(pulseIntegral);
      //   fTimeHists.at(key)->Fill(pulseTimeF+1e-4);

      
    } // loop over each pulse in this event
    
  
      // *************************************************************
      // Examples of how to use the output logging function
      // Look at MsgLog: print frequency/situation depends on numbers (1, 2, 3, etc) 
      // *************************************************************
    if (MsgLog::GetGlobalLogLevel()) {
      MsgDebug(2,MsgLog::Form("Print one number %d and a second %d",fMyLowValue, fMyHighValue));
    }
    
    if (MsgLog::GetGlobalLogLevel()) {
      MsgDebug(1,"Veto prompt integral");
    }


  } // loop over all events in this trigger

  
} // TemplateAna

