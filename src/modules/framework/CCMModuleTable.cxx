/*------------------------------------------------------

  CCMModuleTable

  Based on the TPCModuleTable from the NIFFTE Experiment

  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#include "CCMModuleTable.h"
#include "MsgLog.h"

CCMModuleTable* CCMModuleTable::fInstance = 0;

//------------------------------------------------------
CCMModuleTable& CCMModuleTable::Instance()
{
  //return the singleton
  if (fInstance == 0) fInstance = new CCMModuleTable();
  return *fInstance;

}

//------------------------------------------------------
void CCMModuleTable::Insert(const char* name, ModuleMaker_t p)
{
  //Insert a module maker into the table
  std::string n(name);
  ModMakerMap_t::iterator itr = fMakerMap.find(n);
  if (itr == fMakerMap.end()) {
    fMakerMap[n] = p;
    return;
  }
  MsgError(MsgLog::Form("Proxy %s already in table",name));

}
  
//------------------------------------------------------
ModuleMaker_t CCMModuleTable::Lookup(const char* name)
{
  //Find a module maker in the table
  std::string n(name);
  ModMakerMap_t::iterator itr = fMakerMap.find(n);
  if (itr != fMakerMap.end()) return *(itr->second);
  
  MsgError(MsgLog::Form("Proxy %s is not in table",name));
  MsgError("*== Proxies in list are:");
  itr = fMakerMap.begin();
  ModMakerMap_t::iterator itrEnd(fMakerMap.end());
  for(; itr!=itrEnd; itr++) {
    MsgError(MsgLog::Form("\t %s",(itr->first).c_str()));
  }
  MsgFatal("Fatal error.  Exiting.");

  return 0;

}
