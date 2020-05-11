#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <fstream>
#include <stdint.h>

#include "TApplication.h"

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
    << "  -c configfile.xml     : Name of XML config file for job processing (REQUIRED)\n"
    << "  -l outputLogFileName  : The name of the log file to save the output default is none"
    << "  -i file.root          : Add input data file\n"
    << "  -I file_list.txt      : Indirect file list\n"
    << "  -o file.root          : Set output data file\n"
    << "  -d #                  : Debugging level of output\n"
    << "  -n #                  : Set number of events to process\n";

  MsgError((char*)(usage->str().c_str()));
  delete usage;
  return;
}

//Read a file containing a list of files 
void IndirectFileList(const char* file, std::vector<std::string>& infileList)
{
  // Open the file and pull the file names in
  std::ifstream infile;
  infile.open(file);
  if (!infile) {
    MsgError(MsgLog::Form("File %s not found.",file));
    exit(1);
  }

  std::string fname;
  std::string line;
  while (infile >> line) {
    // from the indirect file
    fname = line;
    infileList.push_back(fname);
  }
  infile.close();
}

//main program
int main (int argc, char** argv) 
{

  std::string cfgfile;
  std::string outputLogFile = "";
  std::vector<std::string> outfileList;
  std::vector<std::string> infileList;
  int debug = 0;        //debugging level
  int32_t nevents = -1; //process all events in file unless otherwise directed

  static const int kConfigOpt   = 'c';
  static const int kInputOpt    = 'i';
  static const int kIndirectOpt = 'I';
  static const int kOutputOpt   = 'o';
  static const int kOutputLogOpt   = 'l';
  static const int kDebugOpt    = 'd';
  static const int kNevtOpt     = 'n';
  static const int kHelpOpt     = 'h';
  static struct option long_options[] = {
    {"config",       required_argument, 0, kConfigOpt},
    {"input",        required_argument, 0, kInputOpt},
    {"indirect",     required_argument, 0, kIndirectOpt},
    {"output",       required_argument, 0, kOutputOpt},
    {"logout",       required_argument, 0, kOutputLogOpt},
    {"debug",        required_argument, 0, kDebugOpt},
    {"nevt",         required_argument, 0, kNevtOpt},
    {"help",         no_argument,       0, kHelpOpt},
    { NULL, 0, 0, 0} // This is a filler for -1 
  };

  while (1) {
    int c;
    int optindx = 0;
    c = getopt_long(argc, argv, "c:i:o:l:I:d:n:h", long_options, &optindx);

    if (c==-1) break;
    std::string fname;
    switch (c) {
      case kConfigOpt:    cfgfile           = std::string(optarg);            break;
      case kInputOpt:     fname             = std::string(optarg); 
                          infileList.push_back(fname);                        break;
      case kIndirectOpt:  IndirectFileList(optarg,infileList);                break;
      case kOutputOpt:    fname             = std::string(optarg); 
                          outfileList.push_back(fname);                       break;
      case kOutputLogOpt: outputLogFile     = std::string(optarg);            break;
      case kDebugOpt:     debug             = std::atoi(optarg);              break;
      case kNevtOpt:      nevents           = std::atoi(optarg);              break;
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

  if(cfgfile.empty()) {// || infileList.empty()) {
    MsgError("Problem reading input parameters");
    Usage();
    return -1;
  }

  MsgLog::SetGlobalDebugLevel(debug);
  MsgLog::SetPrintRepetitions(false);
  if (!outputLogFile.empty()) {
    MsgLog::SetFileOutput(outputLogFile.c_str());
  }

  if(outfileList.empty()) {
    std::string outfile = "ccm-output.root";
    outfileList.push_back(outfile);
    //MsgInfo(MsgLog::Form("Setting up default output file: %s",outfile.c_str()));
  }

  CCMTaskManager* reco = new CCMTaskManager(cfgfile,infileList,outfileList);

  // check to see if process type is event display
  // if so start the TApplication in order to show
  // the event display
  TApplication * app = 0;
  if(reco->GetTaskConfig()->ProcessType() == "EventDisplay")
    app = new TApplication("EventDisplay",&argc,argv);

  reco->Execute(nevents);

  reco->Terminate(); //writes data to output file

  if(app)
    app->Run();

  delete reco;
  delete MsgLog::Instance();
  if(app)
    delete app;

  return 0;

}

