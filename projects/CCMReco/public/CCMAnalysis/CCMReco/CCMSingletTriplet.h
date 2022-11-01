#ifndef CCMSingletTriplet_h
#define CCMSingletTriplet_h

#include <vector>
#include <iterator>
#include <algorithm>

#include "CCMAnalysis/CCMFramework/CCMModule.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

class TFile;
class TH1;
class TH1D;
class TH1F;
class TH2F;
class TTree;

// ***************************************************************
// Module SingletTriplet, a template analysis module
// ***************************************************************
class CCMSingletTriplet : public CCMModule
{
 public:
  
  // the constructor
  CCMSingletTriplet(const char* version);
  
  // copy constructor
  CCMSingletTriplet(const CCMSingletTriplet& clufdr);
  
  // destructor
  ~CCMSingletTriplet();
  
  
  // ***************************************************************
  // Configures things that are hardware and event specific settings
  // Sets all input parameters based on values in .xml config file
  // If not present in .xml file, sets to default values in .cc file
  // ***************************************************************
  void Configure(const CCMConfig& c);
  
  
  // ***************************************************************
  // This is where the action takes place
  // Main body of code, run over every entry in the input .root file
  // Takes as input triggers found by FindEvent module
  // (can be multiple events per trigger, event = detector activity)
  // ***************************************************************
  CCMResult_t ProcessTrigger();
  
  
  // ***************************************************************  
  // Returns true of the job has ended
  // Finish any job related tasks like printing global counters
  // ***************************************************************
  CCMResult_t EndOfJob();
  
  
  // ***************************************************************  
  // Methods to pass the data to process into the module
  // These are generic for every module
  // DON'T EDIT THESE LINES!
  // ***************************************************************
  void ConnectEvents(std::shared_ptr<Events> evt) { fEvents = evt; }
  void ConnectRawData(std::shared_ptr<RawData> rawData)  { fRawData = rawData; }
  void ConnectPulses(std::shared_ptr<Pulses> pulses) {fPulses = pulses; }
  void ConnectOutFileName(std::string name) { fOutFileName = name; SetupOutFile();}
  // should have already run FindEvent module to get specific events/triggers
  // you need for your analysis, however you may need to check the trigger
  // type, etc.
  void SetTriggerType(std::string triggerType) { fTriggerType = triggerType; }
  
  
  // *************************************************************** 
  // Module-specific functions to set input parameters
  // ex: use double to set threshold, string to set trigger type, bool to
  //     set a flag to do/not do parts of code
  // ***************************************************************
  void SetSingletMinPE(double singletMinPE) {fSingletMinPE = singletMinPE; }
  void SetSingletTimeWidth(double singletTimeWidth) {fSingletTimeWidth = singletTimeWidth; }
  
  void SetTripletMinPE(double tripletMinPE) {fTripletMinPE = tripletMinPE; }
  void SetTripletTimeWidth(double tripletTimeWidth) {fTripletTimeWidth = tripletTimeWidth; }

  
  void SetTreeName(std::string treeName) { fTreeName = treeName; }
  
 private:
  // private methods
  void SetupOutFile();
  
  // ***************************************************************
  // Define vectors & histograms used & reset them & prettify them
  // ***************************************************************
  void DefineVarsAndHistos();
  void ResetVarsAndHistos();
  void PrettyHist(TH1 *hist, const char *myX, const char *myY, int statFlag);
  
  
  // ***************************************************************
  // Analysis code you want applied for every entry in the input file
  // run inside of ProcessTrigger, only to keep ProcessTrigger
  // concise.  Technically could put all code in SingletTriplet into
  // ProcessTrigger function
  // ***************************************************************
  void SingletTriplet();
  
  
 private:
  // private data members
  
  // ***************************************************************
  // Need these parameters to: access events, raw data, pulses and output to a file
  // Define trigger type, keep track of total triggers processed
  // Generic for every module
  // ***************************************************************
  std::shared_ptr<Events>  fEvents;
  std::shared_ptr<RawData> fRawData;
  std::shared_ptr<Pulses>  fPulses;
  
  std::string fOutFileName;
  
  std::string fTriggerType;
  
  unsigned long int fNumTriggers;
  
  
  // ***************************************************************  
  // Module-specific vars based on your input parameters in .xml
  // ***************************************************************
  double fSingletMinPE;
  double fSingletTimeWidth;

  double fTripletMinPE;
  double fTripletTimeWidth;

  int fSaveDebugPlots;
  int fPrintDebugStatements;
  
  // var needed to put tree in the output file
  std::string fTreeName;
 
  
  // ***************************************************************
  // Output file & histos written to file
  // ***************************************************************
  TFile * fOutfile;
  TTree * fTree;

  
  // *******************************************************************************
  // declare vars for ntuple and ntuple here
  // *******************************************************************************
  long evt;
  double singInt;
  double tripInt;
  double ratio;
  double thistime;
  int singletLowBin;
  int singletHighBin;
  int eventEndBin;
  int tripletLowBin;
  int tripletHighBin;
  int notrip;

  
  // *******************************************************************************
  // declare histograms, arrays of histograms
  // these are for debug / low level info purposes
  // *******************************************************************************

  // plot pulse time b/c getting lots > 8000 samples, not sure what this means
  TH1F *h_pulseSample;
  
  // make these plots before I any quality cuts at all on pulses, check if info is present
  TH1F *h_pmtID_all;
  TH2F *h_pmtIDvsTime_all;

  // after check if info is present for pmt
  TH1F *h_pmtID_exist;
  TH2F *h_existPMTLoc;
  TH2F *h_HVexistPMTLoc;

  // plot for veto PMTs; I skip those pulses
  TH1F *h_pmtID_veto;
  TH2F *h_vetoPMTLoc;
  TH2F *h_HVvetoPMTLoc;
 
  // after all event selection, still in pulse loop
  TH1F *h_pmtID;
  TH2F *h_pmtIDvsTime;
  TH2F *h_PMTLoc;
  TH2F *h_HVPMTLoc;
 
  // per-PMT plots for pulses passing selection cuts
  std::vector<TH1D*> h_PMTintegral;
  std::vector<TH1D*> h_PMTfiredTime;
  
  // for pulses passing event selection cuts
  TH1F *h_hasnstimePulse;
  TH1F *h_hasDAQnstimePulse;
  
  TH1F *h_singletpe;
  TH1F *h_tripletpe;
  
  // total event plots
  TH1F *h_totalInt;
  TH1F *h_pmtsWithPulse;

  // histos to investigate events failing b/c max pe isn't in the prompt region
  TH1F *h_maxpe;
  TH1F *h_deltatime;

  
  // *******************************************************************************
  // this is for the analysis
  // *******************************************************************************

  // event integral for all events passing event selection cuts
  TH1F *h_event_int;

  
  // *******************************************************************************
  // global vars
  // *******************************************************************************
  int key = 0;
  int hasPulse[160];
  int hasnstimePulse[8000];
  
  double avgChannelNoise = 0;
  long numTriggers = 0;
  
  int board = 0;
  int channel = 0;

  std::vector<std::vector<double>> pmtWaveform;
  
  // for per-pulse loops, for global event counters
  int passDAQtime, passNStime, passIsActive, passPMTinfo, passVeto, passADCIntegral, passADCthresh;

  
};

#endif // CCMSingletTriplet_h

