/*------------------------------------------------------

  CCMTaskManager

  Based on the TPCTaskManager and NiffteTaskMager from
  the NIFFTE Experiment


Adapter: R. T. Thornton (IU)
Date: 14-May-2013

-------------------------------------------------------*/
#include "CCMTaskManager.h"
#include "CCMRawIO.h"
#include "CCMRootIO.h"


#include <cstdlib>
#include <ctime>

std::unique_ptr<CCMTaskConfig> CCMTaskManager::fgkTaskConfig = nullptr;  //Static instance of config object
//------------------------------------------------------
CCMTaskManager::CCMTaskManager()
{
  //default constructor
  fRootIO = std::make_shared<CCMRootIO>();//CCMRootIO::GetInstance();
  fRawIO = std::make_shared<CCMRawIO>();//CCMRawIO::GetInstance();
  
  fAccumWaveform = std::make_shared<AccumWaveform>();
  fMCTruth = std::make_shared<MCTruth>();
  fEvents = std::make_shared<Events>();
  fPulses = std::make_shared<Pulses>();
  fRawData = std::make_shared<RawData>();
  fBinaryRawData = std::make_shared<RawData>();

  fCurrentInFileName = "";
  fCurrentOutFileName = "";
  fCurrentRunNum = 0;
  fCurrentSubRunNum = 0;
}

//------------------------------------------------------
CCMTaskManager::CCMTaskManager(const CCMTaskManager& task)
{
  //copy constructor
  fRootIO = std::make_shared<CCMRootIO>();//CCMRootIO::GetInstance();
  fRawIO = std::make_shared<CCMRawIO>();//CCMRawIO::GetInstance();
  
  fAccumWaveform = std::make_shared<AccumWaveform>();
  fMCTruth = std::make_shared<MCTruth>();
  fEvents = std::make_shared<Events>();
  fPulses = std::make_shared<Pulses>();
  fRawData = std::make_shared<RawData>();
  fBinaryRawData = std::make_shared<RawData>();

  fCurrentInFileName = "";
  fCurrentOutFileName = "";
  fCurrentRunNum = 0;
  fCurrentSubRunNum = 0;
}

//------------------------------------------------------
CCMTaskManager::~CCMTaskManager() 
{ 
  //destructor
}

//------------------------------------------------------
CCMTaskManager::CCMTaskManager(std::string configfile, 
    std::vector<std::string> rootInfileList,
    std::vector<std::string> rootOutfileList,
    std::vector<std::string> rawInfileList,
    std::vector<std::string> rawOutfileList)
{

  fRootIO = std::make_shared<CCMRootIO>();//CCMRootIO::GetInstance();
  fRawIO = std::make_shared<CCMRawIO>();//CCMRootIO::GetInstance();

  SetTaskConfig(new CCMTaskConfig(configfile,rootInfileList,rootOutfileList,
        rawInfileList,rawOutfileList,fRootIO,fRawIO));

  CCMResult_t val = RegisterModules();
  if(val != kCCMSuccess)
  {
    MsgError("Error registering modules: ");
    auto mList = GetTaskConfig().ModuleList();
    for (auto & p : mList) {
      MsgError(MsgLog::Form("\t %s %s",p.first.c_str(),p.second.c_str()));
    }
    MsgError("Exiting gracefully");
    exit(0);
  }

  // default set the objects, they may be just
  // the default constructor values, but it makes
  // it easier so that there are no condition statements
  // in the ConnectDataToModules() function
  fAccumWaveform = std::make_shared<AccumWaveform>();
  fMCTruth = std::make_shared<MCTruth>();
  fEvents = std::make_shared<Events>(fRootIO->GetEvents());
  fPulses = std::make_shared<Pulses>(fRootIO->GetPulses());
  fRawData = std::make_shared<RawData>(fRootIO->GetRawData());
  fBinaryRawData = std::make_shared<RawData>(fRawIO->GetRawData());

  fCurrentInFileName = "";
  fCurrentOutFileName = "";
  fCurrentRunNum = 0;
  fCurrentSubRunNum = 0;
}

//------------------------------------------------------
CCMResult_t CCMTaskManager::Execute(int32_t nevt)
{
  CCMResult_t status = kCCMSuccess;

  if(nevt == 0) {
    return status;
  }

  //if(fgkTaskConfig->ProcessType() == "EventDisplay")
  //return ExecuteTask();

  auto rootFileList = fgkTaskConfig->InputFileList();
  auto rawFileList = fgkTaskConfig->RawInputFileList();
  auto rootFileListOut = fgkTaskConfig->OutputFileList();
  auto rawFileListOut = fgkTaskConfig->RawOutputFileList();

  bool outRoot = true;
  if (rootFileListOut.empty() || !rawFileListOut.empty()) {
    outRoot = false;
  }

  if (!rootFileListOut.empty()) {
    fCurrentOutFileName = rootFileListOut.front();
  } else {
    fCurrentOutFileName = rawFileListOut.front();
  }

  if (!rootFileList.empty()) {
    return ExecuteRoot(nevt,rootFileList,outRoot);
  } else if (!rawFileList.empty()) {
    return ExecuteRaw(nevt,rawFileList,outRoot);
  }

  return status;
}

//------------------------------------------------------
CCMResult_t CCMTaskManager::ExecuteRaw(int32_t nevt, const std::vector<std::string> & fileList, bool outRoot)
{
  CCMResult_t status = kCCMSuccess;

  MsgInfo(MsgLog::Form("FileListSize %zu",fileList.size()));
  int32_t count = 1;
  int digit = 1;
  int exp = 0;
  int run = 0;
  int subRun = 0;
  for (auto & file : fileList) {
    MsgDebug(2,MsgLog::Form("File Name: %s",file.c_str()));
    fCurrentInFileName = file;

    Utility::ParseStringForRunNumber(file,run,subRun);
    if (run != fCurrentRunNum || subRun != fCurrentSubRunNum) {
      fCurrentRunNum = run;
      fCurrentSubRunNum = subRun;
      NewRun();
    }

    for(; // no initiation
        fRawIO->ReadOK(); // check to see if tree read is still ok
        fRawIO->Advance()) // advance to next event that passes cuts
    {
      if(count == nevt && nevt != -1) {
        break;
      }

      int mod = digit*std::pow(10,exp);
      if (count%mod == 0) {
        MsgInfo(MsgLog::Form("[%d] %s /%d/ %s",count,Utility::tstamp(),
              fRawIO->GetTriggerNumber(), fRawIO->CurrentFileName()));
        ++digit;
        if (digit == 10) {
          digit = 1;
          ++exp;
        }
      } // end count module

      fBinaryRawData->operator=(fRawIO->GetRawData());

      ConnectDataToModules();

      status = ExecuteTask();

      if (status != kCCMFailure && status != kCCMDoNotWrite) {
        if (outRoot) {
          fRootIO->SetAccumWaveform(*fAccumWaveform);
          fRootIO->SetMCTruth(*fMCTruth);
          fRootIO->SetEvents(*fEvents);
          fRootIO->SetRawData(*fRawData);
          fRootIO->SetPulses(*fPulses);
          fRootIO->WriteTrigger();
        } else {
          fRawIO->SetRawData(*fRawData);
          fRawIO->WriteTrigger();
        }
      }
      ++count;

      // Don't think I need to do these
      //ClearDataVectors();
    }

    if(count == nevt && nevt != -1) {
      break;
    }

    fRawIO->AdvanceFile();
  }

  return status;
}

//------------------------------------------------------
CCMResult_t CCMTaskManager::ExecuteRoot(int32_t nevt, const std::vector<std::string> & fileList, bool outRoot)
{
  CCMResult_t status = kCCMSuccess;

  MsgInfo(MsgLog::Form("FileListSize %zu",fileList.size()));
  int32_t count = 0;
  int digit = 1;
  int exp = 0;
  int run = 0;
  int subRun = 0;
  for (auto & file : fileList) {
    MsgDebug(2,MsgLog::Form("File Name: %s",file.c_str()));
    fCurrentInFileName = file;

    Utility::ParseStringForRunNumber(file,run,subRun);
    if (run != fCurrentRunNum || subRun != fCurrentSubRunNum) {
      fCurrentRunNum = run;
      fCurrentSubRunNum = subRun;
      NewRun();
    }

    for(; // no initiation
        fRootIO->ReadOK(); // check to see if tree read is still ok
        fRootIO->Advance()) // advance to next event that passes cuts
    {
      if(count == nevt && nevt != -1) {
        break;
      }

      int mod = digit*std::pow(10,exp);
      if ((count+1)%mod == 0) {
        MsgInfo(MsgLog::Form("[%d] %s /%d/ %s",count+1,Utility::tstamp(),
              fRootIO->GetTriggerNumber(), fRootIO->CurrentFileName()));
        ++digit;
        if (digit == 10) {
          digit = 1;
          ++exp;
        }
      } // end count module

      fAccumWaveform->operator=(fRootIO->GetAccumWaveform());
      fMCTruth->operator=(fRootIO->GetMCTruth());
      fEvents->operator=(fRootIO->GetEvents());
      fRawData->operator=(fRootIO->GetRawData());
      fPulses->operator=(fRootIO->GetPulses());

      ConnectDataToModules();

      status = ExecuteTask();

      if (status != kCCMFailure && status != kCCMDoNotWrite) {
        if (outRoot) {
          fRootIO->SetAccumWaveform(*fAccumWaveform);
          fRootIO->SetMCTruth(*fMCTruth);
          fRootIO->SetEvents(*fEvents);
          fRootIO->SetRawData(*fRawData);
          fRootIO->SetPulses(*fPulses);
          fRootIO->WriteTrigger();
        } else {
          fRawIO->SetRawData(*fRawData);
          fRawIO->WriteTrigger();
        }
      }
      ++count;

      // Don't think I need to do these
      //ClearDataVectors();
    }

    if(count == nevt && nevt != -1) {
      break;
    }

    fRootIO->AdvanceFile();
  }

  return status;
} // end ExecuteRoot

//------------------------------------------------------
CCMResult_t CCMTaskManager::Terminate()
{
  CCMResult_t status = FinishTask();

  fRootIO->Close();
  fRawIO->Close();

  fRootIO.reset();
  fRawIO.reset();

  return status;
}

//------------------------------------------------------
CCMResult_t CCMTaskManager::RegisterModules()
{
  //Register and configure requested modules
  if(!fgkTaskConfig) {
    MsgFatal("No configuration object available.  Aborting");
  }

  auto mList = fgkTaskConfig->ModuleList();
  for(size_t i = 0; i < mList.size(); ++i) {
    std::string name = mList.at(i).first;
    std::string version = mList.at(i).second;
    ModuleMaker_t mm = CCMModuleTable::Instance().Lookup(name.c_str());
    if(mm == 0) {
      MsgFatal(MsgLog::Form("Attempt to request a non-existing module %s.  Aborting.",name.c_str()));
    }
    fModuleList.push_back(std::shared_ptr<CCMModule>((*mm)(version.c_str()))); // build the module
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
  for (auto & module : fModuleList) {
    module->ConnectBinaryRawData(fBinaryRawData);
    module->ConnectAccumWaveform(fAccumWaveform);
    module->ConnectMCTruth(fMCTruth);
    module->ConnectEvents(fEvents);
    module->ConnectRawData(fRawData);
    module->ConnectPulses(fPulses);
    module->ConnectInFileName(fCurrentInFileName);
    module->ConnectOutFileName(fCurrentOutFileName);
  }

}

//------------------------------------------------------
void CCMTaskManager::NewRun()
{
  //Call the NewRun function for the modules as some of them
  //do special things at the start/end of each run/subrun
  //combination
  for (auto & module : fModuleList) {
    module->NewRun(fCurrentRunNum,fCurrentSubRunNum);
  }

}

//------------------------------------------------------
CCMResult_t CCMTaskManager::ExecuteTask()
{
  //loop over registered modules and execute their ProcessTrigger() methods
  CCMResult_t status = kCCMSuccess;
  int ctr = 0;
  for (auto & module : fModuleList) {
    status = module->ProcessTrigger();
    ++ctr;

    if(status == kCCMFailure) {
      break;
    }
  } // end for range-based for loop for fModuleList

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

