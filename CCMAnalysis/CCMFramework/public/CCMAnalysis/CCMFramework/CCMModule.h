/*------------------------------------------------------

  CCMModule

  Based on TPCModule from the NIFFTE Experiment 

  Base class for SciBath modules


  Adapter: R. T. Thornton (LANL)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMMODULE_H
#define CCMMODULE_H

#include <memory>
#include <string>
#include <vector>
#include <stdint.h>

#include "CCMAnalysis/utils/Utility.h"

class CCMModule; //needed for typedef
class CCMConfig;
class Events;
class AccumWaveform;
class RawData;
class Pulses;
class MCTruth;

typedef CCMModule* (*ModuleMaker_t)(const char* version);

class CCMModule 
{

public:

  CCMModule(const char* name); //constructor
  CCMModule(const CCMModule& mod); //copy constructor

  virtual ~CCMModule() {}; //destructor

  virtual CCMResult_t ProcessTrigger(); // Process an individual event

  virtual void Configure(const CCMConfig& c); // Beginning of job

  // Start new run
  virtual CCMResult_t NewRun(uint32_t run, uint32_t subrun);

  //Finish any job related tasks
  virtual CCMResult_t EndOfJob();

  virtual void CheckInit();

  //Methods to pass the data to process into the module
  virtual void ConnectAccumWaveform(std::shared_ptr<AccumWaveform> accumWaveform);
  virtual void ConnectMCTruth(std::shared_ptr<MCTruth> mcTruth);
  virtual void ConnectEvents(std::shared_ptr<Events> evt);
  virtual void ConnectBinaryRawData(std::shared_ptr<RawData> rawData);
  virtual void ConnectRawData(std::shared_ptr<RawData> rawData);
  virtual void ConnectPulses(std::shared_ptr<Pulses> pulses);
  virtual void ConnectInFileName(std::string name);
  virtual void ConnectOutFileName(std::string name);

  const char* Name()    const { return fName.c_str();       }
  const char* Version() const { return fCfgVersion.c_str(); }
  uint32_t CurrentRun() const { return fCurrentRun; }
  uint32_t CurrentSubRun() const { return fCurrentSubRun; }

  void SetCfgVersion(const char* cfgv);

protected:
  std::string  fName;       ///< Name of module
  std::string  fCfgVersion; ///< Version of configuration in use
  std::string  fInFileName; ///< Name of the current input file
  std::string  fOutFileName; ///< Name of the current output file
  bool         fIsInit;     ///< Init to false constructor
  uint32_t     fCurrentRun; ///< Current run number
  uint32_t     fCurrentSubRun; ///< Current subrun number

};

#endif //CCMMODULE_H
