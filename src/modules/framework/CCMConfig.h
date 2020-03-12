/*------------------------------------------------------

  CCMConfig

  Based on CCMConfig from the NIFFTE Experiment which is
  based on CfgConfig from the E907/MIPP Experiment

  A collection of configuration parameters

  Configurations are known by their name and version. Their source
  records where the data used to file them comes from (eg. XML file
  name

  Configurations once created are not modifiable within a job
  
  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMCONFIG_H
#define CCMCONFIG_H

#include <map>
#include <string>

class CCMConfigParam;

class CCMConfig
{

public:

  //These methods are normally all you'll ever need
  const char* Name() const { return fName.c_str(); }
  const char* Version() const { return fVersion.c_str(); }
  const char* Source() const { return fSource.c_str(); }

  const CCMConfigParam& operator()(const char*        pname) const;
  const CCMConfigParam& operator()(const std::string& pname) const;

public:
  
  // These methods should *only* be used by
  // the configuration builder
  CCMConfig(const char* name, 
            const char* version, 
            const char* source);
  ~CCMConfig();

  void      AdoptParam(CCMConfigParam* p);
  void      RemoveParam(const char* which = "*");
  CCMConfigParam& Param(const char*        pname);
  CCMConfigParam& Param(const std::string& pname);

  void      SetSource(const char* filename) { fSource = filename; }

  static void Copy(CCMConfig& dest, 
                   const CCMConfig& src,
                   const char* destVersion=0,
                   const char* destSource =0);

public:

  typedef std::map<std::string, CCMConfigParam*> CCMConfigParamMap;
  CCMConfigParamMap& ParamMap() { return fParam; }

private:

  std::string               fName;    ///< Name of configuration
  std::string               fVersion; ///< Version tag of configuration
  std::string               fSource;  ///< Source of config. data (file name eg.)
  mutable CCMConfigParamMap fParam;   ///< Collection of parameters

};

#endif //CCMCONFIG_H
