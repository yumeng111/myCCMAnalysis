/*------------------------------------------------------

  CCMModule

  Based on TPCModule from the NIFFTE Experiment 

  Base class for SciBath modules


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#include "CCMModule.h"
#include "MsgLog.h"

//------------------------------------------------------
CCMModule::CCMModule(const char* name)
  : fName(name), fCfgVersion(""), fIsInit(false), fCurrentRun(-1)
{
  //Default constructor
}

//------------------------------------------------------
CCMModule::CCMModule(const CCMModule& mod)
  : fName(mod.fName), fCfgVersion(mod.fCfgVersion), fIsInit(mod.fIsInit), fCurrentRun(mod.fCurrentRun)
{
  //Copy constructor

}
  
//------------------------------------------------------
CCMResult_t CCMModule::ProcessEvent()
{
  //Default Event method for all modules
  MsgError(MsgLog::Form("Calling default %s::Event().  Please overload",fName.c_str()));
  
  return kCCMSuccess;

}

//------------------------------------------------------
CCMResult_t CCMModule::NewRun(uint32_t run)
{
  //Default NewRun method for all modules
  //doesn't do anything except set internal run number
  MsgDebug(2,MsgLog::Form("Calling default %s::NewRun().",fName.c_str()));

  if(fCurrentRun != run) {
    MsgDebug(2,MsgLog::Form("\tChanging run number from %d to %d",fCurrentRun,run));
    fCurrentRun = run;
  }
  
  return kCCMSuccess;

}

//------------------------------------------------------
CCMResult_t CCMModule::EndOfJob()
{
  //Default EndOfJob method for all modules
  MsgError(MsgLog::Form("Calling default %s::EndOfJob().  Please overload.",fName.c_str()));
  
  return kCCMSuccess;

}

//------------------------------------------------------
void CCMModule::Configure(const CCMConfig& /* cfg */)
{
  //Default configure method
  MsgError(MsgLog::Form("Calling default %s::Configure().  Please overload.",fName.c_str()));

  //use CCMConfig object to configure version

}

//------------------------------------------------------  
void CCMModule::SetCfgVersion(const char* cfgv)
{
  if (fCfgVersion == cfgv) return;  //Nothing to do...
  fCfgVersion = cfgv;

}

//------------------------------------------------------  
void CCMModule::CheckInit()
{
  if (fIsInit) return;
  MsgFatal(MsgLog::Form("Module %s/%s is not initialized.  Stopping execution.",fName.c_str(),
			fCfgVersion.c_str()));

}

