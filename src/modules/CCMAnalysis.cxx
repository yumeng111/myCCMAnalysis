#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <stdint.h>

#include "TH1.h"
#include "TApplication.h"

#include "Utility.h"
#include "MsgLog.h"
#include "CCMTaskManager.h"

//
extern "C" {
#include <unistd.h>
#include <sys/time.h>
#include "getopt.h"
}
// ROOT
//#include "TRint.h"
//#include "TSystem.h"
//#include "TROOT.h"
//#include "TDatabasePDG.h"
//

void Usage() {

  std::stringstream* usage = new std::stringstream;
  (*usage) << "Usage: CCMAnalysis [options] \n"
    << "options are:\n"
    << "  -c configfile.xml     : Name of XML config file for job processing\n"
    << "  -i file.root          : Add input root data file\n"
    << "  -I file_list.txt      : Indirect root file list\n"
    << "  -o file.root          : Set output root data file\n"
    << "  -r file.bin           : Add input binary data file\n"
    << "  -R file_list.txt      : Indirect binary file list\n"
    << "  -b file.bin           : Set output binary data file\n"
    << "  -l outputLogFileName  : The name of the log file to save the output (default is none)\n"
    << "  -d #                  : Debugging level of output\n"
    << "  -n #                  : Set number of events to process (default is all)\n"
    << "must include -c and either (-i/-I) or  (-r/-R) but not both\n";

  MsgError((char*)(usage->str().c_str()));
  delete usage;
  return;
}

//main program
int main (int argc, char** argv) 
{

  // set it so square of the sum of the weights
  // is used to calculate the error on any histograms
  // that are generated
  TH1::SetDefaultSumw2(true);
  
  // set so ROOT will not manage any memory of ROOT objects
  // that are created
  TH1::AddDirectory(0);

  bool removeFirstSubRun = false;
  std::string cfgfile;
  std::string outputLogFile = "";
  std::vector<std::string> outfileList;
  std::vector<std::string> infileList;
  std::vector<std::string> rawOutfileList;
  std::vector<std::string> rawInfileList;
  int debug = 0;        //debugging level
  int32_t nevents = -1; //process all events in file unless otherwise directed

  static const int kConfigOpt   = 'c';
  static const int kInputOpt    = 'i';
  static const int kRawInputOpt = 'r';
  static const int kRawIndirectOpt = 'R';
  static const int kIndirectOpt = 'I';
  static const int kRemoveFirstSubRun = 's';
  static const int kOutputOpt   = 'o';
  static const int kRawOutputOpt = 'b';
  static const int kOutputLogOpt  = 'l';
  static const int kDebugOpt    = 'd';
  static const int kNevtOpt     = 'n';
  static const int kHelpOpt     = 'h';
  static struct option long_options[] = {
    {"config",       required_argument, 0, kConfigOpt},
    {"input",        required_argument, 0, kInputOpt},
    {"indirect",     required_argument, 0, kIndirectOpt},
    {"rawinput",     required_argument, 0, kRawInputOpt},
    {"rawindirect",  required_argument, 0, kRawIndirectOpt},
    {"removefirstsubrun", no_argument,  0, kRemoveFirstSubRun},
    {"output",       required_argument, 0, kOutputOpt},
    {"rawoutput",    required_argument, 0, kRawOutputOpt},
    {"logout",       required_argument, 0, kOutputLogOpt},
    {"debug",        required_argument, 0, kDebugOpt},
    {"nevt",         required_argument, 0, kNevtOpt},
    {"help",         no_argument,       0, kHelpOpt},
    { NULL, 0, 0, 0} // This is a filler for -1 
  };

  while (1) {
    int c;
    int optindx = 0;
    c = getopt_long(argc, argv, "c:i:o:l:I:d:n:r:R:b:sh", long_options, &optindx);

    if (c==-1) break;
    std::string fname;
    switch (c) {
      case kConfigOpt:    cfgfile           = std::string(optarg);            break;
      case kInputOpt:     fname             = std::string(optarg); 
                          infileList = Utility::GetListOfFiles(optarg);       break;
      case kRawInputOpt:  fname             = std::string(optarg); 
                          rawInfileList.push_back(fname);                     break;
      case kIndirectOpt:  Utility::IndirectFileList(optarg,infileList);                break;
      case kRawIndirectOpt:  Utility::IndirectFileList(optarg,rawInfileList);          break;
      case kOutputOpt:    fname             = std::string(optarg); 
                          outfileList.push_back(fname);                       break;
      case kRawOutputOpt: fname             = std::string(optarg); 
                          rawOutfileList.push_back(fname);                    break;
      case kOutputLogOpt: outputLogFile     = std::string(optarg);            break;
      case kDebugOpt:     debug             = std::atoi(optarg);              break;
      case kNevtOpt:      nevents           = std::atoi(optarg);              break;
      case kRemoveFirstSubRun: removeFirstSubRun = true;                      break;
      case kHelpOpt:      Usage(); exit(0);                                   break;
      default:
                          MsgError(MsgLog::Form("Unknown option %d %s",optind,argv[optind]));
                          Usage();
                          exit(1);
    }
  }
  for (; optind<argc; ++optind) {
    infileList.push_back(std::string(argv[optind]));
  }

  if (cfgfile.empty() || (infileList.empty() && rawInfileList.empty())) {
    MsgError("Problem reading input parameters or no input file supplied");
    Usage();
    return EXIT_SUCCESS;
  }

  if (removeFirstSubRun) {
    MsgInfo("Going to remove the first subrun");
    auto it = infileList.begin();
    while ((it = std::find_if(infileList.begin(),infileList.end(),
            [](const std::string& s) { return s.find("000000") != std::string::npos;})) != infileList.end()) {
      infileList.erase(it);
    }

    if (infileList.empty() && rawInfileList.empty()) {
      MsgFatal("No more input files after removing those that match the 0th sub run");
    }
  }

  MsgLog::SetGlobalDebugLevel(debug);
  MsgLog::SetPrintRepetitions(false);
  if (!outputLogFile.empty()) {
    MsgInfo(MsgLog::Form("Setting output log file %s",outputLogFile.c_str()));
    MsgLog::SetFileOutput(outputLogFile.c_str());
  }

  if (!infileList.empty() && !rawInfileList.empty()) {
    MsgError("Both binary and root input file(s) are set, only one type should be set");
    Usage();
    return EXIT_SUCCESS;
  }

  if (!outfileList.empty() && !rawOutfileList.empty()) {
    MsgError("Both binary and root output file(s) are set, only one type should be set");
    Usage();
    return EXIT_SUCCESS;
  }

  // Default is to not save anything but only print to screen
  //if(outfileList.empty()) {
  //  std::string outfile = "ccm-output.root";
  //  outfileList.push_back(outfile);
  //  //MsgInfo(MsgLog::Form("Setting up default output file: %s",outfile.c_str()));
  //}
  
  std::unique_ptr<CCMTaskManager> taskMan = 
    std::make_unique<CCMTaskManager>(cfgfile,infileList,outfileList,rawInfileList,rawOutfileList);

  // check to see if process type is event display
  // if so start the TApplication in order to show
  // the event display
  std::unique_ptr<TApplication> app = nullptr;
  if(taskMan->GetTaskConfig().ProcessType() == "EventDisplay") {
    app = std::unique_ptr<TApplication>(new TApplication("EventDisplay",&argc,argv));
  }

  taskMan->Execute(nevents);

  taskMan->Terminate(); //writes data to output file

  if(app != nullptr) {
    app->Run();
  }

  delete MsgLog::Instance();

  return EXIT_SUCCESS;

}

