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
#include <memory>

class CCMRootIO;
class CCMRawIO;

class CCMTaskConfig {

public:

  CCMTaskConfig(); //default constructor
  CCMTaskConfig(const CCMTaskConfig& reco); //copy constructor
  virtual ~CCMTaskConfig(); //destructor

  CCMTaskConfig(std::string configfile, 
      std::vector<std::string> rootInfileList, std::vector<std::string> rootOutfileList,
      std::vector<std::string> rawInfileList, std::vector<std::string> rawOutfileList,
      std::shared_ptr<CCMRootIO> rootIO, std::shared_ptr<CCMRawIO> rawIO);
  CCMTaskConfig(std::string configfile,std::shared_ptr<CCMRootIO> rootIO, std::shared_ptr<CCMRawIO> rawIO); //primary constructor

  int ReadConfigFile(std::shared_ptr<CCMRootIO> rootIO, std::shared_ptr<CCMRawIO> rawIO);
  void Print();
  int Split(const char* line, const char* tok, std::vector<std::string>& fields);

  std::string ConfigFile() const { return fConfigFile; }

  std::vector<std::string> InputFileList() const { return fInputFileList; }
  std::vector<std::string> OutputFileList() const { return fOutputFileList; }
  std::string InputFileName(unsigned int i) const { if (fInputFileList.size() > i) return fInputFileList[i]; else return "";}
  std::string OutputFile(unsigned int i) const { if (fOutputFileList.size() > i) return fOutputFileList[i]; else return "";}
  std::vector<std::string> RawInputFileList() const { return fRawInputFileList; }
  std::vector<std::string> RawOutputFileList() const { return fRawOutputFileList; }
  std::string RawInputFileName(unsigned int i) const { if (fRawInputFileList.size() > i) return fRawInputFileList[i]; else return "";}
  std::string RawOutputFile(unsigned int i) const { if (fRawOutputFileList.size() > i) return fRawOutputFileList[i]; else return "";}

  std::string ProcessType() const { return fProcessType; } 
  std::vector<std::pair<std::string,std::string>> ModuleList() const { return fModuleList; }
  std::string ConfigInfo() const { return fConfigInfo;}

protected:

private:

  std::string fConfigFile;                   //name of config file
  std::vector<std::string> fInputFileList;   //list of input root file names
  std::vector<std::string> fOutputFileList;   //list of output root file names
  std::vector<std::string> fRawInputFileList;   //list of input binary file names
  std::vector<std::string> fRawOutputFileList;   //list of output binary file names
  std::string fProcessType;                  //type of run to process:
					     //DetSim, Reco, Ana, etc.
  std::vector<std::pair<std::string,std::string>> fModuleList;      //list of modules to be run
  std::string fConfigInfo;                   //string containing all
					     //config info for this job
};

#endif //CCMTASKCONFIG_H
