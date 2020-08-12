#ifndef CCMTemplateAna_h
#define CCMTemplateAna_h

#include "Utility.h"
#include "CCMModule.h"

#include <vector>
#include <iterator>
#include <algorithm>

class TFile;
class TH1D;
class TTree;

// ***************************************************************
// Module TemplateAna, a template analysis module
// ***************************************************************
class CCMTemplateAna : public CCMModule
{
 public:
  
  // the constructor
  CCMTemplateAna(const char* version);
  
  // copy constructor
  CCMTemplateAna(const CCMTemplateAna& clufdr);
  
  // destructor
  ~CCMTemplateAna();
  
  
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
  void SetMyDouble(double myDouble) { fMyDouble = myDouble; }
  //    void SetThreshold(double threshold) { fThreshold = threshold; }
  
  void SetMyString(std::string myString) { fMyString = myString; }
  //    void SetTriggerType(std::string triggerType) { fTriggerType = triggerType; }
  
  void SetMyBool(bool myBool) { fMyBool = myBool; }
  //    void SetPulseTimeCut(bool flag) { fPulseTimeCut = flag; }
  
  
  // parameters that are set if your Bool is true
  void SetMyLowValue(double myValue) { fMyLowValue = myValue; }
  //    void SetPulseTimeLowValue(double value) { fPulseTimeLowValue = value; }
  
  void SetMyHighValue(double myValue) { fMyHighValue = myValue; }
  //    void SetPulseTimeHighValue(double value) { fPulseTimeHighValue = value; }
  
  
  void SetTreeName(std::string treeName) { fTreeName = treeName; }
  
 private:
  // private methods

  // ***************************************************************
  // !!!! VITAL !!!! for making sure histos/ntuple written to file
  // ***************************************************************
  void SetupOutFile();
  
  
  // ***************************************************************
  // Define vectors & histograms used & reset them
  // ***************************************************************
  void DefineVectorsAndHistos();
  void ResetVectorsAndHistos();
  
  
  // ***************************************************************
  // Analysis code you want applied for every entry in the input file
  // run inside of ProcessTrigger, only to keep ProcessTrigger
  // concise.  Technically could put all code in TemplateAna into
  // ProcessTrigger function
  // ***************************************************************
  void TemplateAna();
  
  
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
  double fMyDouble;
  //  double fThreshold;
  
  std::string fMyString;
  //  std::string fTriggerType;
  
  // vars used to apply or not apply analysis cuts
  int fMyBool;
  //  int fPulseTimeCut;
  
  
  // vars that are set/used if fMyBool == TRUE
  double fMyLowValue;
  //  double fPulseTimeLowValue;
  double fMyHighValue;
  //  double fPulseTimeHighValue;
  
  
  // var needed to put tree in the output file
  std::string fTreeName;
  
  
  // ***************************************************************
  // Example vector for body of code (eg TemplateAna function)
  // ***************************************************************
  std::vector<float> fMyVector;
  //   std::vector<float> fPulsesTime;
  
  
  // ***************************************************************
  // Output file & histos written to file
  // ***************************************************************
  TFile * fOutfile;
  TTree * fTree;
  

  // histos that go in the output file
  TH1D* fMySingleHist;
  //   TH1D* fPreEventHist;
  
  std::vector<TH1D*> fMyHists;
  //   std::vector<TH1D> fTimeHists;
  
  
  // vars that go in the output ntuple tree
  double fmyOutDouble;
  // double fEnergy;

  
};

#endif // CCMTemplateAna_h

