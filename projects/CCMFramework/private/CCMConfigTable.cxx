/*------------------------------------------------------

  CCMConfigTable

  Based on TPCConfigTable from the NIFFTE Experiment which is
  based on CfgTable from the E907/MIPP Experiment

  A collection of configuration objects

  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/

#include "CCMAnalysis/CCMFramework/CCMConfigTable.h"

#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"

CCMConfigTable* CCMConfigTable::fInstance = 0;

//------------------------------------------------------
CCMConfigTable& CCMConfigTable::Instance()
{
  //return the sole instance of this class
  if (fInstance == 0) fInstance = new CCMConfigTable;
  return (*fInstance);
}

//------------------------------------------------------
CCMConfig* CCMConfigTable::GetConfig(const char* name, const char* version) const
{
  //Find and return the configuration that matches the requested name
  //and version
  std::string n(name);
  std::string v(version);
  
  // Loop over configurations looking for one that matches
  ConfigList::iterator itr   (fConfigList.begin());
  ConfigList::iterator itrEnd(fConfigList.end());
  for (; itr!=itrEnd; ++itr) {
    CCMConfig* cfg = *itr;
    if (n == cfg->Name() && v == cfg->Version()) return cfg;
  }

  // Not found
  MsgError(MsgLog::Form("Module %s.%s not found in configuration table!",
            name,version));
  return 0;
}

//------------------------------------------------------
void CCMConfigTable::Print() const 
{
  // Print all the configurations
  MsgInfo("Printing the Configuration table");

  ConfigList::const_iterator itr(fConfigList.begin());
  ConfigList::const_iterator itrEnd(fConfigList.end());
  for (; itr!=itrEnd; ++itr) {
    CCMConfig* c = *itr;
    MsgInfo(MsgLog::Form("%s.%s",c->Name(),c->Version()));
  }

}

//------------------------------------------------------
void CCMConfigTable::AdoptConfig(CCMConfig* cfg) 
{
  // Place the configuration into the list

  std::string n(cfg->Name());
  std::string v(cfg->Version());
  
  // Check if a configuration matching this name and version exists.
  // If yes, replace it
  ConfigList::iterator itr   (fConfigList.begin());
  ConfigList::iterator itrEnd(fConfigList.end());
  for (; itr!=itrEnd; ++itr) {
    CCMConfig* inlist = *itr;
    if (n == inlist->Name() && v == inlist->Version()) {
      // replace existing config with new
      delete inlist;
      *itr = cfg;
      return;
    }
  }
  
  // Insert the configuration into the list
  fConfigList.push_back(cfg);
}
