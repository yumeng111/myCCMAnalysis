/*------------------------------------------------------

  CCMConfig

  Based on TPCConfig from the NIFFTE experiment which is
  based on CfgConfig from the E907/MIPP Experiment

  A collection of configuration parameters

  Configurations are known by their name and version. Their source
  records where the data used to file them comes from (eg. XML file
  name

  Configurations once created are not modifiable within a job

  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/

#include <sstream>

#include "CCMAnalysis/CCMFramework/CCMConfig.h"

#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"

//------------------------------------------------------
CCMConfig::CCMConfig(const char* name, const char* version, const char* source) 
  : fName    ( name    ),
    fVersion ( version ),
    fSource  ( source  )
{ 
  //default constructor

}

//------------------------------------------------------
CCMConfig::~CCMConfig()
{
  //destructor
  this->RemoveParam("*");

}

//------------------------------------------------------
const CCMConfigParam& CCMConfig::operator()(const char* pname) const
{
  // Allows read-only access to parameters

  CCMConfigParamMap::iterator itr = fParam.find(pname);
  if (itr==fParam.end()) {
    MsgError(MsgLog::Form("Parameter %s not found in configuration %s:%s",
			  pname,fName.c_str(),fVersion.c_str()));
    throw 1;
  }
  return (*itr->second);

}

//------------------------------------------------------
const CCMConfigParam& CCMConfig::operator()(const std::string& pname) const
{
  // Allows read-only access to parameters
  CCMConfigParamMap::iterator itr = fParam.find(pname);
  if (itr==fParam.end()) {
    MsgError(MsgLog::Form("Parameter %s not found in configuration %s:%s",
			  pname.c_str(),fName.c_str(),fVersion.c_str()));
    throw 1;
  }

  return (*itr->second);

}

//------------------------------------------------------
CCMConfigParam& CCMConfig::Param(const char* pname) 
{
  //Allows non-const access to parameters
  std::string p(pname);
  CCMConfigParamMap::iterator itr = fParam.find(p);
  if (itr==fParam.end()) {
    MsgError(MsgLog::Form("Parameter %s not found in configuration %s:%s",
			  pname,fName.c_str(),fVersion.c_str()));
    throw 1;
  }

  return (*itr->second);

}

//------------------------------------------------------
CCMConfigParam& CCMConfig::Param(const std::string& pname)
{
  //Allows non-const access to parameters
  CCMConfigParamMap::iterator itr = fParam.find(pname);
  if (itr==fParam.end()) {
    MsgError(MsgLog::Form("Parameter %s not found in configuration %s:%s",
                          pname.c_str(),fName.c_str(),fVersion.c_str()));
    throw 1;
  }

  return (*itr->second);

}

//------------------------------------------------------
void CCMConfig::AdoptParam(CCMConfigParam* p)
{
  // Insert the parameter p into a configuration. If the parameter
  // already exists in the configuration then replace it. Ownership of
  // the parameter stored at p is transferred to the configuration.

  std::string s = p->Name();
  
  // Check the parameter set for one of this name. If one is found
  // replace it
  CCMConfigParam* param = fParam[s];
  if (param==0) delete param;
  
  fParam[s] = p;
}

//------------------------------------------------------
void CCMConfig::RemoveParam(const char* which) 
{
  //Remove a selected parameter from the configuration. "*" removes
  //all parameters

  std::string s(which);
  // If which is "*" remove all parameters
  if (s=="*") {
    CCMConfigParamMap::iterator    itr(fParam.begin());
    CCMConfigParamMap::iterator itrEnd(fParam.end());
    for (; itr!=itrEnd; ++itr) {
      if (itr->second) { delete itr->second; itr->second = 0; }
    }
    return;
  }
  
  // Otherwise just remove the named set
  CCMConfigParam* param = fParam[s];
  if (param==0) { delete param; fParam[s] = 0; }

}

//------------------------------------------------------
void CCMConfig::Copy(CCMConfig& dest, const CCMConfig& src,
                     const char* destVersion, const char* destSource) 
{
  //Copy one configuration to another
  //  dest - configuration to copy to
  //  src  - configuration to copy from
  //  destVersion - version tag appied to dest (default appends '.copy')
  //  destSource  - source tag appied to dest (default appends '(copy)' )

  dest.fName = src.fName;
  if (destVersion==0) { 
    dest.fVersion  = src.fVersion;
    dest.fVersion += ".copy";
  }
  else {
    dest.fVersion = destVersion;
  }
  if (destSource==0) {
    dest.fSource  = src.fSource;
    dest.fSource += "(copy)";
  }
  else {
    dest.fSource = destSource;
  }

  // Wipe the parameter map clean and build a new one using the
  // CCMConfigParam copy constructor on each element in the source map
  dest.RemoveParam("*");
  CCMConfigParamMap::const_iterator    itr(src.fParam.begin());
  CCMConfigParamMap::const_iterator itrEnd(src.fParam.end());
  for (; itr!=itrEnd; ++itr) {
    CCMConfigParam* p = new CCMConfigParam(*itr->second);
    dest.AdoptParam(p);
  }

}

