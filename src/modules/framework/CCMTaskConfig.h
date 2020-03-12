/*------------------------------------------------------

  CCMTaskConfig

  Based on the TPCTaskConfig from the NIFFTE Experiment


  Adapter: R. T. Thornton (IU)
  Date: 14-May-2013

-------------------------------------------------------*/
#ifndef CCMTASKCONFIG_H
#define CCMTASKCONFIG_H

#include <vector>
#include <string>

class CCMTaskConfig {

public:

  CCMTaskConfig(); //default constructor
  CCMTaskConfig(const CCMTaskConfig& reco); //copy constructor
  virtual ~CCMTaskConfig(); //destructor

  CCMTaskConfig(std::string configfile, std::vector<std::string> infileList, std::vector<std::string> outfileList); //primary constructor
  CCMTaskConfig(std::string configfile); //primary constructor

  int ReadConfigFile();
  void Print();
  int Split(const char* line, const char* tok, std::vector<std::string>& fields);

  std::string ConfigFile() const { return fConfigFile; }
  std::vector<std::string> InputFileList() const { return fInputFileList; }
  std::vector<std::string> OutputFileList() const { return fOutputFileList; }
  std::string InputFileName(unsigned int i) const { if (fInputFileList.size() > i) return fInputFileList[i]; else return "";}
  std::string OutputFile(unsigned int i) const { if (fOutputFileList.size() > i) return fOutputFileList[i]; else return "";}
  std::string ProcessType() const { return fProcessType; } 
  std::vector<std::string> ModuleList() const { return fModuleList; }
  std::string ConfigInfo() const { return fConfigInfo;}

protected:

private:

  std::string fConfigFile;                   //name of config file
  std::vector<std::string> fInputFileList;   //list of input file names
  std::vector<std::string> fOutputFileList;   //list of output file names
  std::string fProcessType;                  //type of run to process:
					     //DetSim, Reco, Ana, etc.
  std::vector<std::string> fModuleList;      //list of modules to be run
  std::string fConfigInfo;                   //string containing all
					     //config info for this job
};

#endif //CCMTASKCONFIG_H
