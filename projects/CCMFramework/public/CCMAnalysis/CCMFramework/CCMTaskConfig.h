/*------------------------------------------------------

  CCMTaskConfig

  Based on the TPCTaskConfig from the NIFFTE Experiment


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMTASKCONFIG_H
#define CCMTASKCONFIG_H

#include <memory>
#include <string>
#include <vector>

class CCMRootIO;
class CCMRawIO;

class CCMTaskConfig {

public:

  CCMTaskConfig(); //default constructor
  CCMTaskConfig(const CCMTaskConfig& reco); //copy constructor
  virtual ~CCMTaskConfig(); //destructor

  CCMTaskConfig(std::string configfile, 
      std::vector<std::string> infileList, std::string outfile,
      std::shared_ptr<CCMRootIO> rootIO, std::shared_ptr<CCMRawIO> rawIO);
  CCMTaskConfig(std::string configfile,std::shared_ptr<CCMRootIO> rootIO, std::shared_ptr<CCMRawIO> rawIO); //primary constructor

  int ReadConfigFile(std::shared_ptr<CCMRootIO> rootIO, std::shared_ptr<CCMRawIO> rawIO);
  void Print();
  int Split(const char* line, const char* tok, std::vector<std::string>& fields);

  std::string ConfigFile() const;

  std::vector<std::string> InputFileList() const;
  std::string InputFileName(unsigned int i) const;
  std::string OutputFile() const;

  std::string ProcessType() const;
  std::vector<std::pair<std::string,std::string>> ModuleList() const;
  std::string ConfigInfo() const;

protected:

private:

  std::string fConfigFile;                   //name of config file
  std::vector<std::string> fInputFileList;   //list of input file names
  std::string fOutputFile;                   //name of output file
  std::string fProcessType;                  //type of run to process:
					     //DetSim, Reco, Ana, etc.
  std::vector<std::pair<std::string,std::string>> fModuleList;      //list of modules to be run
  std::string fConfigInfo;                   //string containing all
					     //config info for this job
};

#endif //CCMTASKCONFIG_H
