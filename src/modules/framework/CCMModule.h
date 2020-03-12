/*------------------------------------------------------

  CCMModule

  Based on TPCModule from the NIFFTE Experiment 

  Base class for SciBath modules


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMMODULE_H
#define CCMMODULE_H

#include <stdint.h>
#include <string>

#include "Utility.h"

class CCMModule; //needed for typedef
class CCMConfig;

#include <vector>

typedef CCMModule* (*ModuleMaker_t)(const char* version);

class CCMModule 
{

public:

  CCMModule(const char* name); //constructor
  CCMModule(const CCMModule& mod); //copy constructor

  virtual ~CCMModule() {}; //destructor

  virtual CCMResult_t ProcessEvent(); // Process an individual event

  virtual void Configure(const CCMConfig& c); // Beginning of job

  // Start new run
  virtual CCMResult_t NewRun(uint32_t run);

  //Finish any job related tasks
  virtual CCMResult_t EndOfJob();

  virtual void CheckInit();

  //Methods to pass the data to process into the module
  //virtual void ConnectEvent(CCMEvent* evt);
  //virtual void ConnectHits(HitVec* hits);
  //virtual void ConnectBIBHits(HitVec* hits);
  //virtual void ConnectEventTimeInfo(CCMEventTimeInfo * evt);
  //virtual void ConnectPosTopology(CCMPosTopology * evt);
  //virtual void ConnectWaveInfo(WaveInfoVec* hits);
  //virtual void ConnectBIBWaveInfo(WaveInfoVec* hits);

  const char* Name()    const { return fName.c_str();       }
  const char* Version() const { return fCfgVersion.c_str(); }
  const uint32_t CurrentRun() const { return fCurrentRun; }

  void SetCfgVersion(const char* cfgv);

protected:
  std::string  fName;       ///< Name of module
  std::string  fCfgVersion; ///< Version of configuration in use
  bool         fIsInit;     ///< Init to false constructor
  uint32_t     fCurrentRun; ///< Current run number

};

#endif //CCMMODULE_H
