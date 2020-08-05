/*!**********************************************
 * \file CCMSingletTriplet.cxx
 * \author H. Ray (UF)
 * \date August 5, 2020
 *
 * To perform singlet/triplet ratio analysis
 *
 * !!!! VIP !!!!
 * Users must also edit the CMakeLists.txt and LinkDef.h files
 * in the module subdir where your module resides (eg reco)
 ***********************************************/
// include your .h file here
#include "CCMSingletTriplet.h"

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
#include "SinglePulse.h"
#include "SimplifiedEvent.h"

#include <sstream>
#include <numeric>
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <memory>
#include <math.h>

#include <fstream>
#include <limits>
#include <locale>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TMultiGraph.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"


#include "TF1.h"
#include "TVector3.h"


//See CCMModuleTable for info
MODULE_DECL(CCMSingletTriplet);


// ****************************************************************
// global var used in code applied for every entry in the input file go here
// ****************************************************************

// array of histograms for each PMT (key)
TH1F *h_PMTintegral[160];
TH1F *h_PMTfiredTime[160];

// to print out every event integral passing event selection
TH1F *h_each_event_int[1000];

// for all events where we calculate s/t ratio
double myavg = 0;
double this_ratio = 0;
double num_ratios = 0.0;
 
// global counters to see how many events are cut along the way
int foundStrobe = 0;
int foundDAQtime = 0;
int foundNStime = 0;
int foundIsActive = 0;
int foundPMTinfo = 0;
int foundVeto = 0;
int foundADCIntegral = 0;
int foundADCthresh = 0;
int foundSinglet = 0;
int foundSingletRange = 0;
int foundMaxBin = 0;
int foundTripletRange = 0;


//_______________________________________________________________________________________
// *********************************
// CAUTION CAUTION CAUTION CAUTION
// vars in the constructor & destructor must be in the SAME ORDER as in the header file!
// *********************************
CCMSingletTriplet::CCMSingletTriplet(const char* version) 
  : CCMModule("CCMSingletTriplet"),
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
CCMSingletTriplet::CCMSingletTriplet(const CCMSingletTriplet& clufdr) 
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
CCMSingletTriplet::~CCMSingletTriplet()
{ 
  // destructor
  fMyHists.clear();
  //   fTimeHists.clear();
}


//_______________________________________________________________________________________
void CCMSingletTriplet::Configure(const CCMConfig& c ) 
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

  MsgInfo(MsgLog::Form("### Total Strobe events found in file: %i", foundStrobe ));

  myavg = this_ratio / num_ratios;
  
  MsgInfo(MsgLog::Form(" S/T Avg Ratio for this set of files = %ld, (ratio, num ratios) = %ld, %ld", myavg, this_ratio, num_ratios ));

  MsgInfo(MsgLog::Form("---------------------------------------"));
  // MsgInfo(MsgLog::Form("Total Events:        %ld", nEntries ));
  MsgInfo(MsgLog::Form("Strobe Events:       %ld", foundStrobe ));
  MsgInfo(MsgLog::Form("Valid DAQ Time:      %ld", foundDAQtime ));
  MsgInfo(MsgLog::Form("Valid ns Time:       %ld", foundNStime ));
  MsgInfo(MsgLog::Form("Active PMT:          %ld", foundIsActive ));
  MsgInfo(MsgLog::Form("PMT Info found:      %ld", foundPMTinfo ));
  MsgInfo(MsgLog::Form("Veto PMTs:           %ld", foundVeto ));
  MsgInfo(MsgLog::Form("Valid Integral:      %ld", foundADCIntegral ));
  MsgInfo(MsgLog::Form("Above ADC Threshold: %ld", foundADCthresh ));
  MsgInfo(MsgLog::Form("Candidate Singlet:   %ld", foundSinglet ));
  MsgInfo(MsgLog::Form("Valid Singlet Range: %ld", foundSingletRange ));
  MsgInfo(MsgLog::Form("Max Bin in Prompt:   %ld", foundMaxBin ));
  MsgInfo(MsgLog::Form("Valid Triplet Range: %ld", foundTripletRange ));
  
  // did you do the initialization?
  fIsInit = true;

} // Configure


//_______________________________________________________________________________________
CCMResult_t CCMSingletTriplet::ProcessTrigger()
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
  SingletTriplet();

 
  // at the end of every event fill the ntuple
  if (fOutfile != nullptr) {
    fOutfile->cd();
    fTree->Fill();
  }

  return kCCMSuccess;

} // ProcessTrigger


//_______________________________________________________________________________________
void CCMSingletTriplet::SetupOutFile()
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
CCMResult_t CCMSingletTriplet::EndOfJob() 
{
  // print out any global counters
  MsgInfo(MsgLog::Form("Num Triggers %ld passed trigger type cut",fNumTriggers));

  MsgInfo(MsgLog::Form("### Total Strobe events found in file: %ld", foundStrobe ));

  myavg = this_ratio / num_ratios;
  
  MsgInfo(MsgLog::Form(" S/T Avg Ratio for this set of files = %ld, (ratio, num ratios) = %ld, %ld", myavg, this_ratio, num_ratios ));

  MsgInfo(MsgLog::Form("---------------------------------------"));
  // MsgInfo(MsgLog::Form("Total Events:        %ld", nEntries ));
  MsgInfo(MsgLog::Form("Strobe Events:       %ld", foundStrobe ));
  MsgInfo(MsgLog::Form("Valid DAQ Time:      %ld", foundDAQtime ));
  MsgInfo(MsgLog::Form("Valid ns Time:       %ld", foundNStime ));
  MsgInfo(MsgLog::Form("Active PMT:          %ld", foundIsActive ));
  MsgInfo(MsgLog::Form("PMT Info found:      %ld", foundPMTinfo ));
  MsgInfo(MsgLog::Form("Veto PMTs:           %ld", foundVeto ));
  MsgInfo(MsgLog::Form("Valid Integral:      %ld", foundADCIntegral ));
  MsgInfo(MsgLog::Form("Above ADC Threshold: %ld", foundADCthresh ));
  MsgInfo(MsgLog::Form("Candidate Singlet:   %ld", foundSinglet ));
  MsgInfo(MsgLog::Form("Valid Singlet Range: %ld", foundSingletRange ));
  MsgInfo(MsgLog::Form("Max Bin in Prompt:   %ld", foundMaxBin ));
  MsgInfo(MsgLog::Form("Valid Triplet Range: %ld", foundTripletRange ));
    


  
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
void CCMSingletTriplet::DefineVectorsAndHistos()
{
  fMyVector.resize(Utility::fgkNumBins,0.f);
  //   fPulsesTime.resize(Utility::fgkNumBins,0.f);


  fMySingleHist = new TH1D("MySingleHist",";Energy (PE);Count",8000, -9.92, 6.08);
  //   fPreEventHist = new TH1D("PreEventHist",";Energy (PE);Count",energyBins.size()-1,&energyBins.front());

  
  for (int i=0; i < 3; ++i) {
    fMyHists.push_back(new TH1D(Form("MyHists%d",i),";Time (#mus);Count", 8000, -9.92, 6.08));
    //   fTimeHist.push_back(new TH1D(Form("TimeHist%d",i),";Time (#mus);Count",8000,-9.92,6.08));
  }


  // *******************************************************************************
  // declare histograms, arrays of histograms
  // these are for debug / low level info purposes
  // *******************************************************************************

  // plot pulse time b/c getting lots > 8000 samples, not sure what this means
  h_pulseSample = new TH1F("h_pulseSample", "Sample number of pulse; Sample Number; Occupancy/Sample", 8100, 0, 8100);
  
  // make these plots before I any quality cuts at all on pulses, check if info is present
  h_pmtID_all = new TH1F("h_pmtID_all", "PMT occupancy, all pulses in event; PMT ID; Occupancy per PMT", 160, 0, 160);
  h_pmtIDvsTime_all = new TH2F("h_pmtID_vs_time_all", "ns Time of Pulse (x) vs PMT occupancy (y) (all); Time of Pulse (ns) (100 ns bins, raw time); PMT ID", 18, -1000, 17000, 160, 0, 160);

  // after check if info is present for pmt
  h_pmtID_exist    = new TH1F("h_pmtID_exist", "PMT occupancy, PMTs exist", 160, 0, 160);
  h_existPMTLoc    = new TH2F("h_existPMTloc", "Disc. Board (x) vs Channel (y), PMTs exist", 13, -0.5, 12.5, 16, 0, 16);
  h_HVexistPMTLoc  = new TH2F("h_HVexistPMTloc", "HV Board (x) vs Channel (y), PMTs exist", 5, -0.5, 4.5, 48, -0.5, 47.5);

  // plot for veto PMTs; I skip those pulses
  h_pmtID_veto   = new TH1F("h_pmtID_veto", "PMT occupancy, All Veto PMTs", 160, 0, 160);
  h_vetoPMTLoc   = new TH2F("h_vetoPMTloc", "Disc. Board (x) vs Channel (y), ALL Veto PMTs", 13, -0.5, 12.5, 16, 0, 16);
  h_HVvetoPMTLoc = new TH2F("h_HVvetoPMTloc", "HV Board (x) vs Channel (y), ALL Veto PMTs", 5, -0.5, 4.5, 48, -0.5, 47.5);
 
  // after all event selection, still in pulse loop
  h_pmtID       = new TH1F("h_pmtID", "PMT occupancy, Analysis Pulses", 160, 0, 160);
  h_pmtIDvsTime = new TH2F("h_pmtID_vs_time", "ns Time of Analysis Pulse (x) vs PMT occupancy (y)", 18, -1000, 17000, 160, 0, 160);
  h_PMTLoc      = new TH2F("h_PMTloc", "Disc. Board (x) vs Channel (y), Analysis PMTs", 13, -0.5, 12.5, 16, 0, 16);
  h_HVPMTLoc    = new TH2F("h_HVPMTloc", "HV Board (x) vs Channel (y), Analysis PMTs", 5, -0.5, 4.5, 48, -0.5, 47.5);

  for (int i=0; i < 160; i++) {
    h_PMTintegral.push_back(new TH1D(Form("PMTintegral_%d", i), Form("Integral for PMT (p.e.), PMT %d", i), 1000, 0, 10000));
    h_PMTfiredTime.push_back(new TH1D(Form("PMTfiredTime_%d", i), Form("ns Time of Pulses for PMT, PMT %d", i), 18, -1000, 17000));
    //   fTimeHist.push_back(new TH1D(Form("TimeHist%d",i),";Time (#mus);Count",8000,-9.92,6.08));

  }
  
  // for pulses passing event selection cuts
  h_hasnstimePulse = new TH1F("h_hasnstimePulse", "Number of Pulses in 2 ns Window", 8100, 0, 16200);
  h_hasDAQnstimePulse = new TH1F("h_hasDAQnstimePulse", "Number of Pulses in 2 ns Window", 8000, -9920, 6080);
  
  h_singletpe = new TH1F("h_singletpe", "Singlet Candidate p.e.", 1000, 0, 1000);
  h_tripletpe = new TH1F("h_tripletpe", "Triplet Candidate p.e.", 1000, 0, 10000);
  
  // total event plots
  h_totalInt = new TH1F("h_totalInt", "Total Pulse Integral per Event (p.e.)", 100, 0, 100000);
  h_pmtsWithPulse = new TH1F("h_pmtsWithPulse", "Per Event Number PMTs with Pulses", 160, 0, 160);

  // histos to investigate events failing b/c max pe isn't in the prompt region
  h_maxpe = new TH1F("h_maxpe", "max pe when not in prompt", 2000, 0, 500);
  h_deltatime = new TH1F("h_deltatime", "delta time, max peak and singlet end (ns)", 11000, -6000, 16000);
  

  // *******************************************************************************
  // this is for the analysis
  // *******************************************************************************

  // event integral for all events passing event selection cuts
  h_event_int = new TH1F("h_event_int", "Total p.e. Integral per 2 ns window, Event", 8000, -9920, 6080);


} // DefineVectorsAndHistos


//_______________________________________________________________________________________
void CCMSingletTriplet::ResetVectorsAndHistos()
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
void CCMSingletTriplet::SingletTriplet()
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

  
} // SingletTriplet

