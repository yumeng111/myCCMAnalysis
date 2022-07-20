/*------------------------------------------------------

  CCMConfigTable

  Based on TPCConfigTable from the NIFFTE Experiment which is
  based on CfgTable from the E907/MIPP Experiment

  A collection of configuration objects

  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMCONFIGTABLE_H
#define CCMCONFIGTABLE_H

#include <map>
#include <list>
#include <string>
#include <iostream>

class CCMConfig;

class CCMConfigTable
{

public:

  //! List of configurations
  typedef std::list<CCMConfig*> ConfigList;

  //! Configuration name/version pair
  typedef std::pair<std::string,std::string> NVPair;

  //! List of configuration name/version pair
  typedef std::list<NVPair> NVPairList;

public:
  static CCMConfigTable& Instance();

  CCMConfig* GetConfig   (const char* name, const char* version) const;

  void       Print()                   const;
  
  void        AdoptConfig(CCMConfig* cfg);
  ConfigList& GetConfigList() { return fConfigList; }

private:

  static  CCMConfigTable*  fInstance;   //!< Sole instance of the table class
  mutable ConfigList       fConfigList; //!< List of configurations defined

};

#endif //CCMCONFIGTABLE_H
