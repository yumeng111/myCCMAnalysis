/*------------------------------------------------------

  CCMTaskManager

  Based on the TPCTaskManager and NiffteTaskMager from
  the NIFFTE Experiment


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMTASKMANAGER_H
#define CCMTASKMANAGER_H

#include <string>
#include <vector>
#include <list>
#include <memory>

//Include these her so we don't have to in the derived classes
#include "CCMTaskConfig.h"
#include "CCMModule.h"
#include "CCMModuleTable.h"
#include "MsgLog.h"
#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMConfigTable.h"

#include "Utility.h"

#include "MsgLog.h"

#include "AccumWaveform.h"
#include "MCTruth.h"
#include "Events.h"
#include "RawData.h"
#include "Pulses.h"

class CCMRootIO;
class CCMRawIO;

class CCMTaskManager {

public:

  CCMTaskManager(); //default constructor
  CCMTaskManager(const CCMTaskManager& task); //copy constructor
  CCMTaskManager(std::string configfile, 
      std::vector<std::string> rootInfileList,
      std::vector<std::string> rootOutfileList,
      std::vector<std::string> rawInfileList,
      std::vector<std::string> rawOutfileList); // primary constructor
  virtual ~CCMTaskManager(); //destructor

  CCMResult_t Execute(int32_t nevt);
  CCMResult_t Terminate();

  static const CCMTaskConfig& GetTaskConfig() { return *fgkTaskConfig; }
protected:

  CCMResult_t RegisterModules(); //these are called internally

private:

  void ConnectDataToModules();
  void NewRun();

  static void SetTaskConfig(CCMTaskConfig* tskCfg) { fgkTaskConfig = std::unique_ptr<CCMTaskConfig>(tskCfg);}

  std::list<std::shared_ptr<CCMModule>> GetModuleList() const { return fModuleList;}
  
  CCMResult_t ExecuteRaw(int32_t nevt, const std::vector<std::string> & fileList, bool outRoot = true);
  CCMResult_t ExecuteRoot(int32_t nevt, const std::vector<std::string> & fileList, bool outRoot = true);

  CCMResult_t ExecuteTask();
  CCMResult_t FinishTask();

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
  std::shared_ptr<RawData>  fBinaryRawData;
  std::shared_ptr<RawData>  fRawData;
  std::shared_ptr<Pulses>  fPulses;

  std::string fCurrentInFileName;
  std::string fCurrentOutFileName;
  
  int fCurrentRunNum;
  int fCurrentSubRunNum;

};

#endif //CCMTASKMANAGER_H
