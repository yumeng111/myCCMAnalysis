/*------------------------------------------------------

  CCMTaskManager

  Based on the TPCTaskManager and NiffteTaskMager from
  the NIFFTE Experiment


Adapter: R. T. Thornton (IU)
Date: 14-May-2013

-------------------------------------------------------*/

#include <ctime>
#include <cstdlib>
#include <stdexcept>

#include "CCMAnalysis/CCMFramework/CCMTaskManager.h"

#include "CCMAnalysis/CCMIO/CCMRawIO.h"
#include "CCMAnalysis/CCMIO/CCMRootIO.h"
#include "CCMAnalysis/CCMIO/IOUtils.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

std::unique_ptr<CCMTaskConfig> CCMTaskManager::fgkTaskConfig = nullptr;  //Static instance of config object
//------------------------------------------------------
CCMTaskManager::CCMTaskManager()
{
  //default constructor
  fRootIO = std::make_shared<CCMRootIO>();
  fRawIO = std::make_shared<CCMRawIO>();

  fAccumWaveform = std::make_shared<AccumWaveform>();
  fMCTruth = std::make_shared<MCTruth>();
  fEvents = std::make_shared<Events>();
  fPulses = std::make_shared<Pulses>();
  fRawData = std::make_shared<RawData>();

  fCurrentInFileName = "";
  fCurrentOutFileName = "";
  fCurrentRunNum = 0;
  fCurrentSubRunNum = 0;
}

//------------------------------------------------------
CCMTaskManager::CCMTaskManager(const CCMTaskManager& task)
{
  //copy constructor
  fRootIO = std::make_shared<CCMRootIO>();
  fRawIO = std::make_shared<CCMRawIO>();

  fAccumWaveform = std::make_shared<AccumWaveform>();
  fMCTruth = std::make_shared<MCTruth>();
  fEvents = std::make_shared<Events>();
  fPulses = std::make_shared<Pulses>();
  fRawData = std::make_shared<RawData>();

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
    std::vector<std::string> infileList,
    std::string outfile)
{

  fRootIO = std::make_shared<CCMRootIO>();
  fRawIO = std::make_shared<CCMRawIO>();

  SetTaskConfig(new CCMTaskConfig(configfile, infileList, outfile, fRootIO,fRawIO));

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
  // TODO switch on which IO we are using
  fRawData = std::make_shared<RawData>(fRootIO->GetRawData());

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

  std::vector<std::string> input_file_list = fgkTaskConfig->InputFileList();
  std::string output_file = fgkTaskConfig->OutputFile();

  if(output_file == "") {
    throw std::runtime_error("No output file specified!");
  } else {
    fCurrentOutFileName = output_file;
  }

  return Execute(nevt, input_file_list);
}

CCMResult_t CCMTaskManager::Execute(int32_t n_events, std::vector<std::string> const & file_list) {
  CCMResult_t status = kCCMSuccess;

  MsgInfo(MsgLog::Form("FileListSize %zu",file_list.size()));
  int run = 0;
  int subRun = 0;

  // Setup output file
  fRootIO->SetOutFileName(fCurrentOutFileName);
  fRootIO->SetupOutputFile();

  Utility::ExponentialCounter event_counter;

  for(auto & file_name : file_list) {
    MsgDebug(2,MsgLog::Form("File Name: %s", file_name.c_str()));

    if(not FileExists(file_name)) {
      std::stringstream ss;
      ss << "File \"" << file_name << "\" does not exist!";
      MsgError(ss.str());
    }


    Utility::ParseStringForRunNumber(file_name, run, subRun);
    if (run != fCurrentRunNum || subRun != fCurrentSubRunNum) {
      fCurrentRunNum = run;
      fCurrentSubRunNum = subRun;
      NewRun(run,subRun);
    }

    CCMFileType file_type = DetermineFileType(file_name);
    SetNextFile(file_type, file_name);

    // Loop until we have processed n_events, or loop indefinitely if n_events is negative
    while(n_events < 0 or (int)(event_counter.Count()) <= n_events) {
      // Continue only if the reader for this file_type has more events
      if(not ReadOK(file_type))
        break;

      // Read in the next event
      ReadTrigger(file_type);

      // Increment the counter and print out periodically
      event_counter.Increment();
      if(event_counter) {
        MsgInfo(MsgLog::Form("[%d] %s /%d/ %s",event_counter.Count(),Utility::tstamp(),
              GetTriggerNumber(file_type), fCurrentInFileName.c_str()));
      }

      // Assign the pointers to each module
      ConnectDataToModules();

      // Run all the modules
      status = ExecuteTask();

      if (status != kCCMFailure && status != kCCMDoNotWrite) {
        // Set the pointers correctly for output
        fRootIO->SetAccumWaveform(*fAccumWaveform);
        fRootIO->SetMCTruth(*fMCTruth);
        fRootIO->SetEvents(*fEvents);
        fRootIO->SetRawData(*fRawData);
        fRootIO->SetPulses(*fPulses);

        // Write the trigger to a root file
        fRootIO->WriteTrigger();
      }

      // Load the next event
      NextEvent(file_type);
    }
  }
  return kCCMSuccess;
}

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
void CCMTaskManager::NewRun(int run, int subrun)
{
  //Call the NewRun function for the modules as some of them
  //do special things at the start/end of each run/subrun
  //combination
  for (auto & module : fModuleList) {
    module->NewRun(run,subrun);
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

void CCMTaskManager::SetNextFile(CCMFileType file_type, std::string const & fname) {
  fCurrentInFileName = fname;
  switch(file_type) {
    case CCMFileType::ROOT:
      fRootIO->SetInFileName(fname);
      fRootIO->SetupInputFile();
      break;
    case CCMFileType::RawBinary:
      fRawIO->SetInFileName(fname);
      fRawIO->SetupInputFile();
      break;
    default:
      std::stringstream ss;
      ss << "Received unknown file type: " << fname;
      MsgError(ss.str());
      break;
  }
}

void CCMTaskManager::NextEvent(CCMFileType file_type) {
  switch(file_type) {
    case CCMFileType::ROOT:
      fRootIO->Advance();
      break;
    case CCMFileType::RawBinary:
      fRawIO->Advance();
      break;
    default:
      std::stringstream ss;
      ss << "Received unknown file type!";
      MsgError(ss.str());
      break;
  }
}

bool CCMTaskManager::ReadOK(CCMFileType file_type) const {
  switch(file_type) {
    case CCMFileType::ROOT:
      return fRootIO->ReadOK();
      break;
    case CCMFileType::RawBinary:
      return fRawIO->ReadOK();
      break;
    default:
      std::stringstream ss;
      ss << "Received unknown file type!";
      MsgError(ss.str());
      return false;
      break;
  }
}

void CCMTaskManager::ReadTrigger(CCMFileType file_type) {
  switch(file_type) {
    case CCMFileType::ROOT:
      fAccumWaveform->operator=(fRootIO->GetAccumWaveform());
      fMCTruth->operator=(fRootIO->GetMCTruth());
      fEvents->operator=(fRootIO->GetEvents());
      fRawData->operator=(fRootIO->GetRawData());
      fPulses->operator=(fRootIO->GetPulses());
      break;
    case CCMFileType::RawBinary:
      fRawData->operator=(fRawIO->GetRawData());
      break;
    default:
      std::stringstream ss;
      ss << "Received unknown file type!";
      MsgError(ss.str());
      break;
  }
}

uint32_t CCMTaskManager::GetTriggerNumber(CCMFileType file_type) const {
  switch(file_type) {
    case CCMFileType::ROOT:
      return fRootIO->GetTriggerNumber();
      break;
    case CCMFileType::RawBinary:
      return fRawIO->GetTriggerNumber();
      break;
    default:
      std::stringstream ss;
      ss << "Received unknown file type!";
      MsgError(ss.str());
      return 0;
      break;
  }
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

