/*------------------------------------------------------

  CCMTaskManager

  Based on the TPCTaskManager and NiffteTaskMager from
  the NIFFTE Experiment


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMTASKMANAGER_H
#define CCMTASKMANAGER_H

#include <list>
#include <memory>
#include <string>
#include <vector>

// Include these here so we don't have to in the derived classes
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMModule.h"
#include "CCMAnalysis/CCMFramework/CCMTaskConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMConfigTable.h"
#include "CCMAnalysis/CCMFramework/CCMModuleTable.h"

#include "CCMAnalysis/CCMDataStructures/Events.h"
#include "CCMAnalysis/CCMDataStructures/Pulses.h"
#include "CCMAnalysis/CCMDataStructures/MCTruth.h"
#include "CCMAnalysis/CCMDataStructures/RawData.h"
#include "CCMAnalysis/CCMDataStructures/AccumWaveform.h"

#include "CCMAnalysis/CCMIO/IOUtils.h"

class CCMRootIO;
class CCMRawIO;

class CCMTaskManager {

public:

  CCMTaskManager(); //default constructor
  CCMTaskManager(const CCMTaskManager& task); //copy constructor
  CCMTaskManager(std::string configfile,
      std::vector<std::string> infileList,
      std::string outfile); // primary constructor
  virtual ~CCMTaskManager(); //destructor

  CCMResult_t Execute(int32_t nevt);
  CCMResult_t Terminate();

  static const CCMTaskConfig& GetTaskConfig() { return *fgkTaskConfig; }
protected:

  CCMResult_t RegisterModules(); //these are called internally

private:

  void ConnectDataToModules();
  void NewRun(int run, int subrun);

  static void SetTaskConfig(CCMTaskConfig* tskCfg) { fgkTaskConfig = std::unique_ptr<CCMTaskConfig>(tskCfg);}

  std::list<std::shared_ptr<CCMModule>> GetModuleList() const { return fModuleList;}

  CCMResult_t Execute(int32_t nevt, const std::vector<std::string> & fileList);
  CCMResult_t ExecuteTask();
  CCMResult_t FinishTask();

  void SetNextFile(CCMFileType file_type, std::string const & fname);
  void NextEvent(CCMFileType file_type);
  bool ReadOK(CCMFileType file_type) const;
  void ReadTrigger(CCMFileType file_type);
  uint32_t GetTriggerNumber(CCMFileType file_type) const;


  void ClearDataVectors();

  static std::unique_ptr<CCMTaskConfig> fgkTaskConfig;  //configuration object
					      //(stores all task
					      //parameters)

  std::list<std::shared_ptr<CCMModule>> fModuleList; //list of registered modules

  std::shared_ptr<CCMRootIO> fRootIO;
  std::shared_ptr<CCMRawIO> fRawIO;

  std::shared_ptr<AccumWaveform>  fAccumWaveform;
  std::shared_ptr<MCTruth>  fMCTruth;
  std::shared_ptr<Events>  fEvents;
  std::shared_ptr<RawData>  fRawData;
  std::shared_ptr<Pulses>  fPulses;

  std::string fCurrentInFileName;
  std::string fCurrentOutFileName;

  int fCurrentRunNum;
  int fCurrentSubRunNum;

};

#endif //CCMTASKMANAGER_H
