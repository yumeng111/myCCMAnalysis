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


// for all events where we calculate s/t ratio
double myavg = 0.0;
double this_ratio = 0.0;
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

    // ntuple vars to be written to output ntuple tree
    evt(0),
    singInt(0),
    tripInt(0),
    ratio(0),
    thistime(0),
    singletLowBin(0),
    singletHighBin(0),
    eventEndBin(0),
    tripletLowBin(0),
    tripletHighBin(0),
    notrip(0)
    
{
  //Default constructor
  this->SetCfgVersion(version);

  DefineVarsAndHistos();

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
  
    // ntuple vars to be written to output ntuple tree
    evt(clufdr.evt),
    singInt(clufdr.singInt),
    tripInt(clufdr.tripInt),
    ratio(clufdr.ratio),
    thistime(clufdr.thistime),
    singletLowBin(clufdr.singletLowBin),
    singletHighBin(clufdr.singletHighBin),
    eventEndBin(clufdr.eventEndBin),
    tripletLowBin(clufdr.tripletLowBin),
    tripletHighBin(clufdr.tripletHighBin),
    notrip(clufdr.notrip)

{
  // copy constructor
}


//_______________________________________________________________________________________
CCMSingletTriplet::~CCMSingletTriplet()
{ 
  // destructor

}


//_______________________________________________________________________________________
void CCMSingletTriplet::Configure(const CCMConfig& c ) 
{
  // Initialize any parameters here by reading them from the CCMConfig object.
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
CCMResult_t CCMSingletTriplet::ProcessTrigger()
{
  // this is where you get the event found by FindEvent module
  // and apply your analysis code

  // since the class does not own the TFile, always check to see if it is present
  fOutfile = gROOT->GetFile(fOutFileName.c_str());
  

  // Check which trigger occured in DAQ window
  // make sure it is the one we want based on input .xml file
  // should have already been done in FindEvent; kept in as a cross-check
  if (!fRawData->IsTriggerPresent(fTriggerType)) return kCCMDoNotWrite;
  else foundStrobe += 1;
  

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

  // ntuple vars to be written to output ntuple tree
  evt = 0;
  singInt = 0;
  tripInt = 0;
  ratio = 0;
  thistime = 0;
  singletLowBin = 0;
  singletHighBin = 0;
  eventEndBin = 0;
  tripletLowBin = 0;
  tripletHighBin = 0;
  notrip = 0;

  
  if (!fOutFileName.empty()) {
    fOutfile = gROOT->GetFile(fOutFileName.c_str());
    if (!fOutfile) {
      MsgWarning(MsgLog::Form("Could not find ROOT file with name %s already opened",fOutFileName.c_str()));
      return;
    }
    fOutfile->cd();
    fTree = new TTree("tree","tree");

    fTree -> Branch("evt", &evt);
    fTree -> Branch("singInt", &singInt);
    fTree -> Branch("tripInt", &tripInt);
    fTree -> Branch("ratio", &ratio);
    fTree -> Branch("thistime", &thistime);
    fTree -> Branch("singletLowBin", &singletLowBin);
    fTree -> Branch("singletHighBin", &singletHighBin);
    fTree -> Branch("eventEndBin", &eventEndBin);
    fTree -> Branch("tripletLowBin", &tripletLowBin);
    fTree -> Branch("tripletHighBin", &tripletHighBin);
    fTree -> Branch("notrip", &notrip);
    
  } // if no output file name

} //SetupOutputFile


//_______________________________________________________________________________________
CCMResult_t CCMSingletTriplet::EndOfJob() 
{
  // print out any global counters
  MsgInfo(MsgLog::Form("Num Triggers %ld passed trigger type cut",fNumTriggers));


  MsgInfo(" ");
  MsgInfo(MsgLog::Form("### Total Strobe events found in file: %ld", foundStrobe ));

  myavg = this_ratio / num_ratios;
  
  MsgInfo(MsgLog::Form(" S/T Avg Ratio for this set of files = %f, (ratio, num ratios) = %ld, %ld", myavg, this_ratio, num_ratios ));

  MsgInfo(MsgLog::Form("---------------------------------------"));
  // MsgInfo(MsgLog::Form("Total Events:        %ld", nEntries ));
  MsgInfo(MsgLog::Form("Strobe Events:       %ld", foundStrobe ));
  MsgInfo(MsgLog::Form("Active PMT:          %ld", foundIsActive ));
  MsgInfo(MsgLog::Form("PMT Info found:      %ld", foundPMTinfo ));
  MsgInfo(MsgLog::Form("Valid DAQ Time:      %ld", foundDAQtime ));
  MsgInfo(MsgLog::Form("Valid ns Time:       %ld", foundNStime ));
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
    fTree -> Write();

    // ########### ONLY IF DEBUG FLAG SET DO THESE ##################
    
    // *******************************************************************************
    // print all histograms to pdf file
    // *******************************************************************************

    // #### DEBUG FLAG 1
    // make these plots before I any quality cuts at all on pulses, check if info is present
    h_pmtID_all -> Write();
    h_pmtIDvsTime_all -> Write();
    
    // after check if info is present for pmt
    h_pmtID_exist -> Write();
    h_existPMTLoc -> Write();
    h_HVexistPMTLoc -> Write();
    
    // plot for veto PMTs; I skip those pulses
    h_pmtID_veto -> Write();
    h_vetoPMTLoc -> Write();
    h_HVvetoPMTLoc -> Write();
    
    // after all event selection, still in pulse loop
    h_pmtID -> Write();
    h_pmtIDvsTime -> Write();
    h_PMTLoc -> Write();
    h_HVPMTLoc -> Write();
    
    // per-PMT plots for pulses passing selection cuts
    for (auto & hist : h_PMTintegral) hist->Write();
    for (auto & hist : h_PMTfiredTime) hist->Write();
    
    // for pulses passing event selection cuts
    h_hasnstimePulse -> Write();
    h_hasDAQnstimePulse -> Write();
   
    // total event plots
    h_totalInt -> Write();
    h_pmtsWithPulse -> Write();

    
    // #### DEBUG FLAG 2
    // plot pulse time b/c getting lots > 8000 samples, not sure what this means
    // see comment by R. T. Thornton explaining these, 6/23/20
    h_pulseSample -> Write();

    
    // #### DEBUG FLAG 4
    // plot of singlet and triplet candidate integrals
    // use for threshold setting information
    h_singletpe -> Write();
    h_tripletpe -> Write();

    
    // #### DEBUG FLAG 5
    // plot of max pe peak relative to prompt singlet
    h_maxpe -> Write();
    h_deltatime -> Write();

    
  } // outfile
    
  return kCCMSuccess;

} // EndOfJob


//_______________________________________________________________________________________
void CCMSingletTriplet::DefineVarsAndHistos()
{
  fMyVector.resize(Utility::fgkNumBins,0.f);
  //   fPulsesTime.resize(Utility::fgkNumBins,0.f);


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
  h_pmtID_exist    = new TH1F("h_pmtID_exist", "PMT occupancy, PMTs exist; PMT ID; Occupancy per PMT", 160, 0, 160);
  h_existPMTLoc    = new TH2F("h_existPMTloc", "Disc. Board (x) vs Channel (y), PMTs exist; Disc. Board; Channel", 13, -0.5, 12.5, 16, 0, 16);
  h_HVexistPMTLoc  = new TH2F("h_HVexistPMTloc", "HV Board (x) vs Channel (y), PMTs exist; HV Board; Channel", 5, -0.5, 4.5, 48, -0.5, 47.5);

  // plot for veto PMTs; I skip those pulses
  h_pmtID_veto   = new TH1F("h_pmtID_veto", "PMT occupancy, All Veto PMTs; PMT ID; Occupancy per PMT", 160, 0, 160);
  h_vetoPMTLoc   = new TH2F("h_vetoPMTloc", "Disc. Board (x) vs Channel (y), ALL Veto PMTs; Disc. Board; Channel", 13, -0.5, 12.5, 16, 0, 16);
  h_HVvetoPMTLoc = new TH2F("h_HVvetoPMTloc", "HV Board (x) vs Channel (y), ALL Veto PMTs; HV Board; Channel", 5, -0.5, 4.5, 48, -0.5, 47.5);
 
  // after all event selection, still in pulse loop
  h_pmtID       = new TH1F("h_pmtID", "PMT occupancy, Analysis Pulses; PMT ID; Occupancy per PMT", 160, 0, 160);
  h_pmtIDvsTime = new TH2F("h_pmtID_vs_time", "ns Time of Analysis Pulse (x) vs PMT occupancy (y); Time of Pulse (ns) (100 ns bins, raw time); PMT ID", 18, -1000, 17000, 160, 0, 160);
  h_PMTLoc      = new TH2F("h_PMTloc", "Disc. Board (x) vs Channel (y), Analysis PMTs; Disc. Board; Channel", 13, -0.5, 12.5, 16, 0, 16);
  h_HVPMTLoc    = new TH2F("h_HVPMTloc", "HV Board (x) vs Channel (y), Analysis PMTs; HV Board; Channel", 5, -0.5, 4.5, 48, -0.5, 47.5);

  for (int i=0; i < 160; i++) {
    h_PMTintegral.push_back(new TH1D(Form("PMTintegral_%d", i), Form("Integral for PMT (p.e.), PMT %d; Integral (p.e.); Occupancy per 100 p.e.", i), 1000, 0, 10000));
    h_PMTfiredTime.push_back(new TH1D(Form("PMTfiredTime_%d", i), Form("ns Time of Pulses for PMT, PMT %d; Time (ns) (raw); Pulses per 100 ns", i), 18, -1000, 17000));
  }
  
  // for pulses passing event selection cuts
  h_hasnstimePulse = new TH1F("h_hasnstimePulse", "Number of Pulses in 2 ns Window; Time of Pulse (ns) (raw, 2 ns windows); Number of Pulses", 8100, 0, 16200);
  h_hasDAQnstimePulse = new TH1F("h_hasDAQnstimePulse", "Number of Pulses in 2 ns Window; Time of Pulse (ns) (DAQ, 2 ns windows); Number of Pulses", 8000, -9920, 6080);
  
  h_singletpe = new TH1F("h_singletpe", "Singlet Candidate p.e.; p.e. of singlet; Number of events", 1000, 0, 1000);
  h_tripletpe = new TH1F("h_tripletpe", "Triplet Candidate p.e.; p.e. of triplet; Number of events", 1000, 0, 10000);
  
  // total event plots
  h_totalInt = new TH1F("h_totalInt", "Total Pulse Integral per Event (p.e.); Total Event Pulse Integral (p.e.); Occupancy per 1000 p.e.", 100, 0, 100000);
  h_pmtsWithPulse = new TH1F("h_pmtsWithPulse", "Per Event Number PMTs with Pulses; Number PMTs per Event with Pulses; Number of Events", 160, 0, 160);

  // histos to investigate events failing b/c max pe isn't in the prompt region
  h_maxpe = new TH1F("h_maxpe", "max pe when not in prompt; Max p.e. when peak not in Prompt; Number of Events", 2000, 0, 500);
  h_deltatime = new TH1F("h_deltatime", "delta time, max peak and singlet end (ns); Delta Time (ns) max peak - singlet (start/end); Number of Events", 11000, -6000, 16000);
  

  // *******************************************************************************
  // this is for the analysis
  // *******************************************************************************

  // event integral for all events passing event selection cuts
  h_event_int = new TH1F("h_event_int", "Total p.e. Integral per 2 ns window, Event", 8000, -9920, 6080);


} // DefineVarsAndHistos


//_______________________________________________________________________________________
void CCMSingletTriplet::ResetVarsAndHistos()
{
  auto itMyVectorBegin = fMyVector.begin();
  //   auto itPulsesTimeBegin = fPulsesTime.begin();
  auto itMyVectorEnd = fMyVector.end();
  //   auto itPulsesTimeEnd = fPulsesTime.end();
  
  std::fill(itMyVectorBegin,itMyVectorEnd,0);
  //   std::fill(itPulsesTimeBegin,itPulsesTimeEnd,0);

  
} // ResetVarsAndHistos


//_______________________________________________________________________________________
void CCMSingletTriplet::SingletTriplet()
{

  // Reset vars, vectors, and histograms used for each trigger
  ResetVarsAndHistos();

  // should all of these go into ResetVarsAndHistos ?
  
  // clear flag vars for global counters
  passDAQtime = 0; passNStime = 0; passIsActive = 0; passPMTinfo= 0;
  passVeto = 0; passADCIntegral = 0; passADCthresh = 0;

  
  // clear the event integral histogram!
  h_event_int -> Reset();

  // for each entry, clear the counter of which PMTs fired a pulse
  for (int i=0; i < 160; i++) hasPulse[i] = 0;

  
  // Required
  auto pmtInfo = PMTInfoMap::GetPMTInfo(0);
  
  
  // *************************************************************
  // Declare variables used locally in analysis
  // ONLY do ones used in pulse loop; declare vars used outside pulse
  // loop when used for better readability/understanding of code
  // *************************************************************
  int key = 0;
  
  double time = 0;
  float myTime = -999;
  double pulseIntegral = 0;

  int myboard = -1;
  int mychan = -1;
  
  int hvboard = -999;
  int hvchan = -999;
  
  // counter of pmt pulse integral for each event
  float sumInt = 0;

  double threshold = -1;
  double pe = -1;

  
  // *************************************************************
  // Get the number of pulses in this trigger for your ana loop
  // *************************************************************
  const size_t kNumPulses = fPulses->GetNumPulses();
 
  
  // *******************************************************************************
  // START PULSE LOOP
  // *******************************************************************************
  
  for (size_t loc = 0; loc < kNumPulses; ++loc) {

    // -------------------------------------------------------
    // get PMT key
    // -------------------------------------------------------
    key = fPulses->GetKey(loc);
 
    
    // -------------------------------------------------------
    // check PMT info, is it in the database
    // -------------------------------------------------------
    // Required: make sure the PMT is to be used for the analysis
    if (!PMTInfoMap::IsActive(key)) {
      // if (debugFlag) std::cout << "*** ERROR: PMT IS NOT IN DB, event: " << evt << " key: " << key << " pulse: " << loc << std::endl;
      continue;
    }

    passIsActive = 1;

	
    // -------------------------------------------------------
    // does it exist for this PMT: channels in digitizer don't have PMTs connected to them,
    // they're not loaded into channel map.  just a sanity check
    // -------------------------------------------------------
    // Required: make sure the PMT info is present
    pmtInfo = PMTInfoMap::GetPMTInfo(key);
    if (!pmtInfo) {
      // std::cout << "*** ERROR: PMT INFO NOT PRESENT" << std::endl;
      continue;
    }

    passPMTinfo = 1;


    // -------------------------------------------------------
    // plot PMT info before any quality cuts on pulses (other than being active & info present)
    // -------------------------------------------------------
    // these are DISCRIMINATOR board & channel, not HV!
    myboard = pmtInfo->GetBoard();
    mychan = pmtInfo->GetBoardChan();
    
    hvboard = -999;
    hvchan = -999;
    PMTInfoMap::ConvertKeyToHVBoardChan(key, hvboard, hvchan);
    
    h_pmtID_exist -> Fill (key);
    h_existPMTLoc -> Fill( myboard, mychan );
    h_HVexistPMTLoc -> Fill( hvboard, hvchan);
    
    
    // -------------------------------------------------------
    // make sure the time is within the range of the DAQ window for analysis
    // does this time account for the correction due to 11th board?
    // -------------------------------------------------------
    time = fPulses->GetPulseTime(loc);
    
    // this error happens a LOT!  why?? values > 8000.  fill histos to explore
    h_pulseSample -> Fill(time);
    
    if (time > 8000) {
      // std::cout << fPulses->GetKey(loc) << std::endl;
      
      // R. T. Thornton 6/23/20
      // this should dominated by 1" PMTs because we shit them in time by 42 ns because of the
      // response time difference between the 8" tubes and the 1" tubes
      // looking at the output that does seem to be the case
    }

    
    // why am I using 7960 here?
    if (std::isnan(time) || std::isinf(time) || time > 7960.0 || time < 0) {
      // this debug statement was b/c we saw lots of PMTs with 8000 ns time - see comment above, this is now understood
      // if (debugFlag) if (time > 7999) std::cout << "*** ERROR: RAW PULSE TIME IS NOT IN RANGE: " << time << std::endl;
      continue;
    }
    
    passDAQtime = 1;
          
   
    // -------------------------------------------------------
    // convert time in samples to time in ns
    // -------------------------------------------------------
    myTime = -999;

    // converts samples into time in ns
    myTime = (time * ( (16 * pow(10,3)) / 8000));
    
    if (myTime == -999 || myTime < 0 || myTime > (16 * pow(10,3))) {
      // if (debugFlag) std::cout << "*** ERROR: BAD NS TIME: " << myTime << std::endl;
      continue;
    }
    
    passNStime = 1;
    

    // -------------------------------------------------------
    // plot PMT info
    // -------------------------------------------------------
    h_pmtID_all -> Fill (key);
    h_pmtIDvsTime_all -> Fill (myTime, key);

    
    // -------------------------------------------------------
    // don't include any Veto PMTs in the analysis
    // plot their ID, disc board/channel, HV board/channel as sanity check
    // -------------------------------------------------------
    if (pmtInfo->IsVeto()) {
      h_pmtID_veto -> Fill (key);
      h_vetoPMTLoc -> Fill( myboard, mychan );
      h_HVvetoPMTLoc -> Fill( hvboard, hvchan );
      continue;
    }
    
    passVeto = 1;
      
   
    // -------------------------------------------------------
    // check that the integrals are realistic (ADC counts)
    // -------------------------------------------------------
    pulseIntegral = fPulses->GetPulseIntegral(loc);
    if (pulseIntegral < 0 || std::isnan(pulseIntegral) || std::isinf(pulseIntegral)) {
      // std::cout << "*** ERROR: BAD PULSE INTEGRAL: " << pulseIntegral << std::endl;
      continue;
    }
    
    passADCIntegral = 1;
    
    
    // -------------------------------------------------------
    // convert integral in ADC to integral in pe
    // using this function returns automatically the calibrated pe
    //
    // threshold check that the pulse is not coming from noise on the digitizer
    // can only do this after SPE calibration is done
    //
    // R. T. Thornton - 6/23/2020
    // -------------------------------------------------------
    threshold = pmtInfo -> GetADCThreshold();
    if (pulseIntegral < threshold) {
      // if (debugFlag) std::cout << "*** ERROR: PULSE BELOW ADC THRESHOLD " << threshold << '\t' << pulseIntegral << std::endl;
      continue;
    }
    
    pe = pmtInfo -> GetADCToPE();
    
    if (pe <= 0) {
      // this happens a lot.  why??  values are -7.07503, no variation
      // if (debugFlag) std::cout << "*** ERROR: BAD PE " << pe << std::endl;
      continue;
    }
    
    pulseIntegral /= pe;
    
    passADCthresh = 1;
    

    // -------------------------------------------------------
    // plot PMT info after all event selection cuts on pulses
    // -------------------------------------------------------
    h_pmtID -> Fill (key);
    h_pmtIDvsTime -> Fill (myTime, key);
    h_PMTLoc -> Fill( myboard, mychan );
    h_HVPMTLoc -> Fill( hvboard, hvchan );
    
    // arrays of integrals for each PMT
    h_PMTintegral[key]->Fill(fPulses->GetPulseIntegral(loc));
    h_PMTfiredTime[key]->Fill(myTime);
    

    // -------------------------------------------------------
    // take individual pulse of each tube, create an accumulated waveform
    // count #PMTs that were on, number of pulses, then plot these two
    //
    // see which PMTs registered pulses this event
    // -------------------------------------------------------
    hasPulse[key] = 1;
    
    // see which bin registered pulses this event
    // casting as int should give the correct bin for my histos (0.5 / 2 = 0, 1.5 / 2 = 0 etc)
    // mytime = ns time of pulse
    int thisbin = h_hasnstimePulse -> FindBin(myTime);
    
    // want to increment bin with this time by 1
    h_hasnstimePulse -> AddBinContent(thisbin, 1); // time = 0 would be bin 1
    
    // now convert time to be DAQ compliant & fill same histo, with shifted axes
    // converting ns into DAQ time
    myTime -= ( (9.92 * pow(10, 3)) );
    
    thisbin = -999; // clear histo bin value of this time
    thisbin = h_hasDAQnstimePulse -> FindBin(myTime);
    
    if (thisbin == -999 || myTime < -9920 || myTime > 6080) {
      //std::cout << "**** WARNING BAD SHIFTED DAQ NS TIME: " << myTime << " , bin " << thisbin << std::endl;
    }
    
    h_hasDAQnstimePulse -> AddBinContent(thisbin, 1); 
    

    // -------------------------------------------------------
    // used for s/t analysis
    // -------------------------------------------------------
    // for each event plot the cumulative puse integral per DAQ time bin
    h_event_int-> AddBinContent(thisbin, pulseIntegral); // adding pulseIntegral from this PMT to the DAQ time bn
    
    // do sum of integral of all pulses in event
    sumInt += pulseIntegral;
    
    
  } // END PULSE LOOP 
    

  // *******************************************************************************
  // increment per-event counters
  // *******************************************************************************
  if (passIsActive == 1) foundIsActive += 1;
  if (passDAQtime == 1) foundDAQtime += 1;
  if (passNStime == 1) foundNStime += 1;
  if (passPMTinfo == 1) foundPMTinfo += 1;
  if (passVeto == 1) foundVeto += 1;
  if (passADCIntegral == 1) foundADCIntegral += 1;
  if (passADCthresh == 1) foundADCthresh += 1;
  
  
  // *******************************************************************************
  // plot total integral of all pulses this event, number of PMTs in event with pulse
  // *******************************************************************************
  h_totalInt -> Fill (sumInt); 
  
  // plot number of PMTs this event with a pulse
  int totalPMTwPulse = 0;
  
  for (int i = 0; i < 160; i++) {
    if (hasPulse[i] != 0) totalPMTwPulse += 1;
  }
  
  h_pmtsWithPulse -> Fill(totalPMTwPulse);
  

  // *******************************************************************************
  // START SINGLET TRIPLET RATIO
  //
  // look at current events' histogram
  // search for singlet/triplet events
  // IGNORE UNDERFLOW BIN
  // *******************************************************************************

  /*
  // check to see if there is a cut on which pulses to include (fMyBool, from xml file)
  // if so, apply cut based on values in .xml file
  time = fPulses->GetPulseTime(loc);
  
  if (fMyBool) {
  if (time < fMyLowValue || time > fMyHighValue) continue;
  } // end if bool cut
  */
    
  
  // -------------------------------------------------------
  // ********** SINGLET **********
  // look for singlet start: first bin with pe >= 5, valid time range
  // FYI cosmics are HUGE, >20 pe or more
  // -------------------------------------------------------
  int foundStart = -1;
  
  for (int i = 1; i <= h_event_int -> GetNbinsX(); i++) {
    if (h_event_int -> GetBinContent(i) >= 5) {
      foundStart = i;
      // if (debugFlag) std::cout << "EVENT: " << evt << " found singlet start bin " << foundStart << std::endl;
    }
    if (foundStart != -1) break;
  } 
  
  if (foundStart == -1) {
    // if (debugFlag) std::cout << "***** INTERESTING: STROBE BUT NO SINGLET FOUND, EVENT " << evt << std::endl;
    return;
  }
  
  foundSinglet += 1;
  
  
  // ------------------------------------------------------------
  // get bin range of singlet
  // 
  // DAQ time range of this bin = bin - 2 ns (low edge) to bin (high edge) ns
  // start with time as low bin and integrate up to 0.09 mu sec b/c 3 sigma of singlet light comes in the first 27 ns sec
  // ------------------------------------------------------------
  int starttime =  (h_event_int -> GetBinLowEdge(foundStart)) + 1; // DAQ time in ns, +1 to put it in middle of bin b/c 2 ns wide bins
  
  int singletwidth = 0.09 * pow(10, 3); // need to convert mu sec into ns, start t = 0 at time of first pe bin and go up 0.09 mu sec (90 ns!) 
  int triprange = (4.0 * pow(10, 3)); // need to convert mu sec into ns, start t = 0 at end of singlet peak window and go up 4000 ns 
  
  singletLowBin = h_event_int -> FindBin(starttime);
  singletHighBin = h_event_int -> FindBin(starttime + singletwidth);
  
  eventEndBin = h_event_int -> FindBin(starttime + singletwidth+ triprange);
  
  
  // ------------------------------------------------------------
  // get max of EVENT = not the same as max of the histo (aka entire DAQ range)!
  // ------------------------------------------------------------
  int maxbin;
  
  h_event_int -> GetXaxis() -> SetRange(singletLowBin,eventEndBin);
  maxbin = h_event_int -> GetMaximumBin();
  h_event_int -> GetXaxis() -> UnZoom();
  
  /*
  // ------------------------------------------------------------
  // **** DEBUG ****
  // print out singlet info before check the bin ranges
  // (counter foundStrobe set to 1 first time we find one, thus need to - 1 to get first histo in array)
  // ------------------------------------------------------------
  if (debugFlag) {
    
    std::cout << "DEBUG: event " << evt << " number of strobe event = " << foundStrobe - 1
	      << " singlet start bin = " << foundStart << " pe value = " << h_event_int -> GetBinContent(foundStart)
	      << " bin low edge (DAQ time, ns) = " << h_event_int -> GetBinLowEdge(foundStart) 
	      << " bin up edge (DAQ time, ns) = " << h_event_int -> GetBinLowEdge(foundStart+1)
	      << " evt max pe bin = " << maxbin << " pe value = " << h_event_int -> GetBinContent(maxbin)
	      << std::endl;
    
    // low edge bin of singlet = bin containing the low edge of the first pe above threshold
    // high bin edge of singlet = bin containing the high edge of the first pe above threshold + singletwidth
    std::cout << "DEBUG: singlet start time (DAQ, ns) = " << starttime << " lowest bin of singlet = " << singletLowBin
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(singletLowBin) << " to " 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(singletLowBin + 1) << std::endl;
    std::cout << "   upper bin of singlet = " << singletHighBin 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(singletHighBin) << " to " 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(singletHighBin + 1) << std::endl;
    
  } // debugFlag for singlet
  */
    
  // ------------------------------------------------------------
  // check bin ranges are acceptable
  // ------------------------------------------------------------
  if ( (h_event_int -> GetBinLowEdge(singletLowBin) < -9920) || (h_event_int -> GetBinLowEdge(singletHighBin + 1) > 6080) ) {
    // if (debugFlag)
    //std::cout << "*** ERROR: bin range for singlet peak outside of DAQ window, Event: " << evt << " Strobe Event: " << foundStrobe - 1
    //		<< " Lower Bin Value = " << h_event_int -> GetBinLowEdge(singletLowBin)
    //		<< " Upper Bin Value = " << h_event_int -> GetBinLowEdge(singletHighBin + 1) << std::endl;
    return;
  }
  
  foundSingletRange += 1;
    
    
  // ------------------------------------------------------------
  // check the max bin lies in the singlet (prompt) region
  // ------------------------------------------------------------
  if (maxbin < singletLowBin || maxbin > singletHighBin) {
    //if (debugFlag) std::cout << "*** ERROR: maxbin is outside of the prompt region, Event: " << evt
    //			     << " Strobe Event: " << foundStrobe - 1 << std::endl;
    // if events fail plot the max pe value
    h_maxpe -> Fill(h_event_int -> GetBinContent(maxbin));
    // also plot the delta time between maxbin and singlet start/end 
    if (maxbin < singletLowBin)  h_deltatime -> Fill(h_event_int -> GetBinLowEdge(maxbin) - starttime);
    if (maxbin > singletHighBin) h_deltatime -> Fill(h_event_int -> GetBinLowEdge(maxbin) - (starttime + singletwidth));					      
    return;
  }
    
  foundMaxBin += 1;
  
  
  // ------------------------------------------------------------
  // get integral for singlet
  // ------------------------------------------------------------
  singInt = h_event_int -> Integral(singletLowBin, singletHighBin);
  
  
  // -------------------------------------------------------
  // ********** TRIPLET **********
  // found valid time range for singlet, now get range for triplet
  // check triplet integral above some min, check time range
  // -------------------------------------------------------
  
  tripletLowBin = singletHighBin + 1; // start right after singlet ends
  
  int tripfullrange = starttime + singletwidth + triprange; // (singlet start time + singlet range = end bin of singlet) + triplet range
  tripletHighBin = h_event_int -> FindBin(tripfullrange);
  
  
  if (h_event_int -> FindBin(tripfullrange) > 8099) {
    // 8100  bins in this histo, 0 = underflow, 8101 = overflow
    // std::cout << "**** WARNING: upper bin of triplet integral outside reasonable bounds, only integrating to 6080 ns DAQ time.  Event "
    //		<< evt << " Strobe Event: " << foundStrobe - 1 << std::endl;
    tripletHighBin = 8100;
  }
  
  /*
  // ------------------------------------------------------------
  // **** DEBUG ****
  // print out our triplet info before we check the bin ranges
  // ------------------------------------------------------------
  if (debugFlag) {
    std::cout << "DEBUG: lowest bin of triplet = " << tripletLowBin 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(tripletLowBin)  << " to " 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(tripletLowBin + 1)  << std::endl;
    std::cout << "   upper bin of triplet = " << tripletHighBin 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(tripletHighBin)  << " to " 
	      << " DAQ time (ns) = " << h_event_int -> GetBinLowEdge(tripletHighBin + 1)  << std::endl;
  }
  */
  
  // ------------------------------------------------------------
  // check bin ranges are acceptable
  // ------------------------------------------------------------
  if ( (h_event_int -> GetBinLowEdge(tripletLowBin) < -9920) || (h_event_int -> GetBinLowEdge(tripletHighBin + 1) > 6080) ) {
    //if (debugFlag)
    //std::cout << "*** ERROR: bin range for triplet peak outside of DAQ window, Event: " << evt << " Event: " << foundStrobe - 1 
    //		<< " Lower Bin Value = " << h_event_int -> GetBinLowEdge(tripletLowBin)
    //		<< " Upper Bin Value = " << h_event_int -> GetBinLowEdge(tripletHighBin + 1) << std::endl;
    return;
  }
  
  foundTripletRange += 1;
  
  
  // ------------------------------------------------------------
  // get integral for triplet, check above threshold
  // if not we still want to keep it to write out to ntuple
  // b/c want to calc efficiency - when find singlet but not triplet
  // ------------------------------------------------------------
  notrip = 0;
  tripInt = h_event_int -> Integral(tripletLowBin, tripletHighBin);
  
  // if (debugFlag) std::cout << "Singlet Integral: " << singInt << " Triplet Integral: " << tripInt << std::endl;
  h_singletpe -> Fill(singInt);
  h_tripletpe -> Fill(tripInt);
  
  double tripThreshold = 5; // because why not, if singlet has min of 5 p.e.
  
  if (tripInt < tripThreshold) {
    notrip = 1;
    // std::cout << "*** ERROR: TRIPLET NOT ABOVE THRESHOLD, Event: "
    //		<< evt << " Strobe Event: " << foundStrobe - 1 << std::endl;
  }

        
  // ------------------------------------------------------------
  // calclulate ratio, set to 0 if no triplet
  // ------------------------------------------------------------
  if (tripInt > 0) ratio = singInt / tripInt;
  else ratio = 0;
  
  // if (debugFlag && ratio > 1) std::cout << " ----------------------- WTF ratio greater than 1!  Event: " << evt << " Ratio: " << ratio << std::endl;

  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // write out events of interest to ntuple
  // *** include pathologies: singlet but no apparent triplet, ratio is crazy high
  //
  // for each entry in the ntuple want: event number (not strobe evt!)
  //     singlet value, triplet value, ratio, time of event,
  //     singlet and triplet low and high bin = daq time for range of integral
  //     flag to tell us if this is an event with no apparent triplet
  //
  // save event integral histogram to ntuple
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // var for ntuple
  thistime = fRawData->GetGPSSecIntoDay(); // this gets sec
  
  h_event_int -> Write();
  
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // do an average ratio for all s/t events in input files
  // skip events with no apparent triplet
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (notrip) return;
  
  this_ratio += ratio;
  num_ratios += 1;
  
  

  // get var from event we want to write out to ntuple
  // fmyOutDouble = simplifiedEvent.GetIntegralTank(true);
  // fEnergy = simplifiedEvent.GetIntegralTank(true);
  
  
  // auto itMyVectorBegin = fMyVector.begin();
  // auto itPulsesTimeBegin = fPulsesTime.begin();
  
  // fPulsesTime.at(bin) += hitFraction;
  
  
  
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
  
  
} // SingletTriplet

