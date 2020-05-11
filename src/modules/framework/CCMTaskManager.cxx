/*------------------------------------------------------

  CCMTaskManager

  Based on the TPCTaskManager and NiffteTaskMager from
  the NIFFTE Experiment


Adapter: R. T. Thornton (IU)
Date: 14-May-2013

-------------------------------------------------------*/
#include "CCMTaskManager.h"
#include <cstdlib>
#include <ctime>

//------------------------------------------------------
const char* tstamp()
{
  //======================================================================
  // Provide a nicely formatted, current, time stamp string
  //======================================================================
  static char tbuff[32];
  time_t t;
  t = time(0);
  strcpy(tbuff, ctime(&t));
  tbuff[24] = '\0';
  return tbuff;
}

const CCMTaskConfig* CCMTaskManager::fgkTaskConfig = 0;  //Static instance of config object
//------------------------------------------------------
CCMTaskManager::CCMTaskManager()
{
  //default constructor
}

//------------------------------------------------------
CCMTaskManager::CCMTaskManager(const CCMTaskManager& task)
{
  //copy constructor
}

//------------------------------------------------------
CCMTaskManager::~CCMTaskManager() 
{ 
  //destructor
}

//------------------------------------------------------
CCMTaskManager::CCMTaskManager(std::string configfile, std::vector<std::string> infileList,
    std::vector<std::string> outfileList)
{
  // currently not using an analysis manager because each module is assumed seperate
  //fCCMAnaMan = CCMAnalysisManager::GetInstance(infileList);
  SetTaskConfig(new CCMTaskConfig(configfile,infileList,outfileList));
  //if(fgkTaskConfig->ProcessType() != "EventDisplay")
  //  fCCMAnaMan->SetOutFile(outfileList.at(0));

  CCMResult_t val = RegisterModules();
  if(val != kCCMSuccess)
  {
    MsgError("Error registering modules: ");
    for(unsigned int i = 0; i < GetTaskConfig()->ModuleList().size(); i++) {
      MsgError(MsgLog::Form("\t %s",GetTaskConfig()->ModuleList().at(i).c_str()));
    }
    MsgError("Exiting gracefully");
    exit(0);
  }

}

//------------------------------------------------------
CCMResult_t CCMTaskManager::Execute(int32_t nevt)
{
  CCMResult_t status = kCCMSuccess;

  if(nevt == 0) {
    return status;
  }

  //if(fgkTaskConfig->ProcessType() == "EventDisplay")
  return ExecuteTask();

  /*
  unsigned int fileListSize = fgkTaskConfig->InputFileList().size();
  MsgInfo(MsgLog::Form("FileListSize %u",fileListSize));
  int32_t count = 0;
  for(unsigned int fileCount = 0; fileCount < fileListSize; ++fileCount, fCCMAnaMan->AdvanceFile())
  {
    MsgInfo(MsgLog::Form("File Name: %s",fCCMAnaMan->CurrentFileName()));

    for(unsigned int c=0; fCCMAnaMan->ReadOK(kBIBTreeID); ++c,fCCMAnaMan->NextEvent(kBIBTreeID))
    {
      CCMHitInformation * hit = fCCMAnaMan->GetHit(0,kBIBTreeID);
      if(hit)
        fBIBHitVec.push_back(new CCMHitInformation(*hit));
    }

    for(fCCMAnaMan->NextEvent(kEventTreeID); // grab first event that passes cuts
        fCCMAnaMan->ReadOK(kEventTreeID); // check to see if tree read is still ok
        fCCMAnaMan->NextEvent(kEventTreeID)) // advance to next event that passes cuts
    {
      if(count == nevt && nevt != -1)
        break;

      if (count%1000 == 0)
        MsgInfo(MsgLog::Form("[%d] %s /%d:%d/ %s",count,tstamp(), fCCMAnaMan->RunNumber(),
              fCCMAnaMan->EventNumber(kEventTreeID), fCCMAnaMan->CurrentFileName()));

      fEvent = fCCMAnaMan->GetEvent();
      fEvtTimeInfo = fCCMAnaMan->GetEventTimeInfo();
      fPosTop = fCCMAnaMan->GetPosTopology();


      unsigned int numhits = fEvent->getnumhits();
      for(unsigned int i=0; i< numhits; ++i)
      {
        CCMHitInformation * hit = fCCMAnaMan->GetHit(i);
        if(hit)
          fHitVec.push_back(new CCMHitInformation(*hit));
      }

      ConnectDataToModules();

      status = ExecuteTask();

      fCCMAnaMan->SetEvent(fEvent);
      fCCMAnaMan->SetHitVec(&fHitVec);
      fCCMAnaMan->SetEventTimeInfo(fEvtTimeInfo);
      fCCMAnaMan->SetPosTopology(fPosTop);
      fCCMAnaMan->Write(kEventTreeID);
      ++count;

      ClearDataVectors();
    }

    fCCMAnaMan->SetBIBHitVec(&fBIBHitVec);
    fCCMAnaMan->Write(kBIBTreeID);

    for(unsigned int i=0; i <fBIBHitVec.size(); ++i) 
      delete fBIBHitVec.at(i);
    if(fBIBHitVec.size() > 0)
      fBIBHitVec.clear();

    fCCMAnaMan->SaveOtherTrees();

    if(count == nevt && nevt != -1)
      break;

  }
  */

  return status;
}

//------------------------------------------------------
CCMResult_t CCMTaskManager::Terminate()
{
  CCMResult_t status = FinishTask();

  //fCCMAnaMan->Close();
  return status;
}

//------------------------------------------------------
CCMResult_t CCMTaskManager::RegisterModules()
{
  //Register and configure requested modules
  if(!fgkTaskConfig)
    MsgFatal("No configuration object available.  Aborting");

  for(unsigned int i = 0; i < fgkTaskConfig->ModuleList().size(); i++) {
    std::string name = fgkTaskConfig->ModuleList().at(i).c_str();
    ModuleMaker_t mm = CCMModuleTable::Instance().Lookup(name.c_str());
    if(mm == 0) {
      MsgFatal(MsgLog::Form("Attempt to request a non-existing module %s.  Aborting.",name.c_str()));
    }
    fModuleList.push_back((*mm)("default")); // build the module

  }

  //Configure modules
  MsgInfo("Configuring modules");
  for (auto & module : fModuleList) {
    module->Configure(*(CCMConfigTable::Instance().GetConfig(module->Name(),module->Version())));
  }

  return kCCMSuccess;

}

//------------------------------------------------------
void CCMTaskManager::ConnectDataToModules()
{
  //Connect data vectors to the registered modules
  //for (auto & module : fModuleList) {
    //module->ConnectEvent(fEvent);
    //module->ConnectHits(&fHitVec);
    //module->ConnectBIBHits(&fBIBHitVec);
    //module->ConnectWaveInfo(&fWaveInfoVec);
    //module->ConnectBIBWaveInfo(&fBIBWaveInfoVec);
    //module->ConnectEventTimeInfo(fEvtTimeInfo);
    //module->ConnectPosTopology(fPosTop);
  //}

}

//------------------------------------------------------
CCMResult_t CCMTaskManager::ExecuteTask()
{
  //loop over registered modules and execute their ProcessEvent() methods
  CCMResult_t status = kCCMSuccess;
  int ctr = 0;
  for (auto & module : fModuleList) {
    status = module->ProcessEvent();
    ++ctr;

    if(status == kCCMFailure) {
      break;
    }
  }

  MsgDebug(2,MsgLog::Form("Processed %d modules",ctr));

  return status;

}


//------------------------------------------------------
CCMResult_t CCMTaskManager::FinishTask()
{
  //Do post-processing
  //loop over registered modules and execute their EndOfJob() methods
  CCMResult_t passed = kCCMSuccess;
  int ctr = 0;
  for (auto & module : fModuleList) {
    passed = module->EndOfJob();
    ++ctr;
    if(passed == kCCMFailure) {
      break;
    }
  }

  MsgDebug(2,MsgLog::Form("Finished job on %d modules",ctr));

  return passed;

}

//------------------------------------------------------
void CCMTaskManager::ClearDataVectors()
{
  //for(unsigned int i=0; i<fHitVec.size(); ++i) delete fHitVec.at(i);
  //for(unsigned int i=0; i<fWaveInfoVec.size(); ++i) delete fWaveInfoVec.at(i);
  //for(unsigned int i=0; i<fBIBWaveInfoVec.size(); ++i) delete fBIBWaveInfoVec.at(i);

  //if(fHitVec.size() > 0) fHitVec.clear();
  //if(fWaveInfoVec.size() > 0) fWaveInfoVec.clear();
  //if(fBIBWaveInfoVec.size() > 0) fBIBWaveInfoVec.clear();

  //if(fEvent)
  //  delete fEvent;
}

