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

//class CCMAnalysisManager;
//class CCMHitInformation;
//class CCMEventTimeInfo;
//class CCMPosTopology;
//class CCMWaveformInformation;

class CCMTaskManager {

public:

  CCMTaskManager(); //default constructor
  CCMTaskManager(const CCMTaskManager& task); //copy constructor
  CCMTaskManager(std::string configfile, std::vector<std::string> infileList,
      std::vector<std::string> outfileList); // primary constructor
  virtual ~CCMTaskManager(); //destructor

  CCMResult_t Execute(int32_t nevt);
  CCMResult_t Terminate();

  static const CCMTaskConfig* GetTaskConfig() { return fgkTaskConfig; }
protected:

  CCMResult_t RegisterModules(); //these are called internally

private:

  void ConnectDataToModules();

  static void SetTaskConfig(CCMTaskConfig* tskCfg) { fgkTaskConfig = tskCfg;}

  std::list<CCMModule*> GetModuleList() const { return fModuleList;}

  CCMResult_t ExecuteTask();
  CCMResult_t FinishTask();

  void ClearDataVectors();

  static const CCMTaskConfig* fgkTaskConfig;  //configuration object
					      //(stores all task
					      //parameters)

  std::list<CCMModule*> fModuleList; //list of registered modules

  //CCMAnalysisManager * fCCMAnaMan;
  //CCMEvent * fEvent;
  //CCMEventTimeInfo * fEvtTimeInfo;
  //CCMPosTopology * fPosTop;
  //HitVec fHitVec;
  //HitVec fBIBHitVec;
  //WaveInfoVec fWaveInfoVec;
  //WaveInfoVec fBIBWaveInfoVec;

};

#endif //CCMTASKMANAGER_H
