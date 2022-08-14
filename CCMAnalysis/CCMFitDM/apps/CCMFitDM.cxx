#include <map>
#include <mutex>
#include <chrono>
#include <future>
#include <random>
#include <thread>
#include <vector>
#include <cstring>
#include <fstream>
#include <numeric>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMTaskConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMConfigTable.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"

#include "TCanvas.h"
#include "TColor.h"
#include "TDecompChol.h"
#include "TF1.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraphSmooth.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TPad.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

extern "C" {
#include <unistd.h>
#include <sys/time.h>
#include "getopt.h"
}

//////////////////////////////////////////////////////////////////
// Global Variables
//////////////////////////////////////////////////////////////////
static int gStartBin = 0;
static int gEndBin = 0;
static int gNumRecoBins = 0;
static size_t gNumNullFakeData = 5000;
static size_t gNumSigFakeData = 1000;
static double gDefaultAlphaD = 0.5;
static double gDefaultEpsilon = 1e-3;
static std::vector<std::tuple<double,double,double,double>> gResults;
static std::mutex gMTX;

static const int gkLightGray = TColor::GetColor(220.f/255.f,220.f/255.f,220.f/255.f);
static const int gkBlack     = TColor::GetColor(0.f,0.f,0.f);
static const int gkOrange    = TColor::GetColor(230.f/255.f,159.f/255.f,0.f);
static const int gkSkyBlue   = TColor::GetColor(86.f/255.f,180.f/255.f,233.f/255.f);
static const int gkBlueGreen = TColor::GetColor(0.f,158.f/255.f,115.f/255.f);
static const int gkBlue      = TColor::GetColor(0.f,114.f/255.f,178/255.f);
static const int gkVermilion = TColor::GetColor(213.f/255.f,94.f/255.f,0.f);
static const int gkRedPurple = TColor::GetColor(204.f/255.f,121.f/255.f,167.f/255.f);
static const int gkYellow    = TColor::GetColor(240.f/255.f,228.f/255.f,66.f/255.f);

//////////////////////////////////////////////////////////////////
// Global Functions
//////////////////////////////////////////////////////////////////
std::shared_ptr<TH1D> GetTotalHist(float mv, float mchi, int axis, std::string histName, 
    std::string masterDir, std::string cvDir, std::vector<std::string> sysDir,
    std::vector<std::shared_ptr<TH1D>> * sys = nullptr);

void CalculateSignalVariance(std::shared_ptr<TH1D> hist, const double omFrac, const double potFrac, const double pi0Frac,
    const double qfFrac) ;

void FindBest(std::ifstream & infile, TH1D *  best);

void AddOld(std::ifstream & infile, std::vector<TGraphAsymmErrors*> & grOld);

void ScaleDMSys(const TMatrixD & cv, const std::vector<TMatrixD> & fracError, TMatrixD & sysError, bool print = false);

double CalculateChi2(const TMatrixD & data, const TMatrixD & pred, const TMatrixD & sys);

void GetToyMC(std::vector<TMatrixD> & fakeDataVec, const TMatrixD & cv, const TMatrixD & sys, std::mt19937 & mt);

void CalculateSensCL(const TMatrixD dm, const TMatrixD & data, const TMatrixD & bkg, 
    const std::vector<TMatrixD> & sysFracError, const TMatrixD & bkgError, const int & index);

struct tokens: std::ctype<char> 
{
  tokens(): std::ctype<char>(get_table()) {}

  static std::ctype_base::mask const* get_table()
  {
    typedef std::ctype<char> cctype;
    static const cctype::mask *const_rc= cctype::classic_table();

    static cctype::mask rc[cctype::table_size];
    std::memcpy(rc, const_rc, cctype::table_size * sizeof(cctype::mask));

    rc[','] = std::ctype_base::space; 
    rc[' '] = std::ctype_base::alpha; 
    return &rc[0];
  }
};

//////////////////////////////////////////////////////////////////
// Main Function
//////////////////////////////////////////////////////////////////
void Usage() {

  std::stringstream* usage = new std::stringstream;
  (*usage) << "Usage: CCMFitDM [options] \n"
    << "options are:\n"
    << "  -c configfile.xml     : Name of XML config file for job processing\n"
    << "  -l outputLogFileName  : The name of the log file to save the output (default is none)\n"
    << "  -d #                  : Debugging level of output\n";

  MsgError((char*)(usage->str().c_str()));
  delete usage;
  return;
}

int main(int argc, char ** argv)
{

  // set it so square of the sum of the weights
  // is used to calculate the error on any histograms
  // that are generated
  TH1::SetDefaultSumw2(true);

  // set so ROOT will not manage any memory of ROOT objects
  // that are created
  TH1::AddDirectory(false);

  ////////////////////////////////////////////////
  // need to remove this block
  std::string rootFileNameCCM120 = "";
  std::string rootFileNameCCM200 = "";
  std::string rootFileNameCCM200Low = "";
  if (argc == 6) {
    rootFileNameCCM120 = std::string(argv[3]);
    rootFileNameCCM200 = std::string(argv[4]);
    rootFileNameCCM200Low = std::string(argv[5]);
  }
  ////////////////////////////////////////////////

  std::string cfgfile = "";
  std::string outputLogFile = "";
  int debug = 0;        //debugging level

  static const int kConfigOpt   = 'c';
  static const int kOutputLogOpt  = 'l';
  static const int kDebugOpt    = 'd';
  static const int kHelpOpt     = 'h';
  static struct option long_options[] = {
    {"config",       required_argument, 0, kConfigOpt},
    {"logout",       required_argument, 0, kOutputLogOpt},
    {"debug",        required_argument, 0, kDebugOpt},
    {"help",         no_argument,       0, kHelpOpt},
    { NULL, 0, 0, 0} // This is a filler for -1 
  };

  while (1) {
    int c;
    int optindx = 0;
    c = getopt_long(argc, argv, "c:l:d:h", long_options, &optindx);

    if (c==-1) break;
    std::string fname;
    switch (c) {
      case kConfigOpt:    cfgfile           = std::string(optarg);            break;
      case kOutputLogOpt: outputLogFile     = std::string(optarg);            break;
      case kDebugOpt:     debug             = std::atoi(optarg);              break;
      case kHelpOpt:      Usage(); exit(0);                                   break;
      default:
                          MsgError(MsgLog::Form("Unknown option %d %s",optind,argv[optind]));
                          Usage();
                          exit(1);
    }
  }

  MsgLog::SetGlobalDebugLevel(debug);
  MsgLog::SetPrintRepetitions(false);
  if (!outputLogFile.empty()) {
    MsgInfo(MsgLog::Form("Setting output log file %s",outputLogFile.c_str()));
    MsgLog::SetFileOutput(outputLogFile.c_str());
  }

  // variables that we will fill from 
  // reading the configuration file
  std::string dmMasterDir = "";
  std::string dmDir = "";
  std::string dmHistName = "";
  int dmHistAxis = 1;
  std::vector<std::string> dmSysDirs;
  std::string massComboFileName = "";
  int resultsFromFile = false;
  int applyEfficiency = true;
  int doFit = true;
  std::vector<std::string> previousResults;
  std::vector<std::string> previousResultsNames;
  std::string dmEffSmearFile = "";
  std::string dmEffHistName = "";
  std::string dmSmearHistName = "";
  std::string dataFileName = "";
  std::string dataHistName = "";
  std::string bkgHistName = "";
  std::string excessHistName = "";
  std::string outputDir = "";
  std::string otherLimitsDir = "";
  std::string postFix = "";

  std::tuple<std::string,std::string,double> inputOMSysInfo;
  std::tuple<std::string,std::string,double> inputQFSysInfo;
  std::tuple<std::string,std::string,double> inputPOTSysInfo;

  double defaultPOT = 1;
  double potScaling = 1;
  double scaleDataAndBkg = 1;
  double startEnergy = 3;
  double endEnergy = 15.5;
  double normalizeBkgValue = -1;

  double mvmChiRatio = 3;
  double newAlphaD = 0.5;

  if (cfgfile.empty()) {
    MsgError("Problem reading configuration file, or did not supply one");
    Usage();
    return EXIT_FAILURE;
  } else {
    MsgInfo(MsgLog::Form("Will read from %s",cfgfile.c_str()));
  }
  auto taskConfig = new CCMTaskConfig(cfgfile,nullptr,nullptr);

  std::string tempString = "";
  //we call the different parts of the fit modules although they are not true
  //modules as in CCMAnalysis uses modules, but since we are using the same 
  //CCMTaskConfig, we will use the same name
  auto mList = taskConfig->ModuleList();
  for(auto & module : mList) {
    std::string name = module.first;
    std::string version = module.second;
    const CCMConfig & c = *(CCMConfigTable::Instance().GetConfig(name.c_str(),version.c_str()));

    if (name == "DMTrueHists") {
      c("MasterDir").Get(dmMasterDir);
      c("Directory").Get(dmDir);
      c("HistName").Get(dmHistName);
      c("HistAxis").Get(dmHistAxis);
      c("DefaultAlphaD").Get(gDefaultAlphaD);
      c("DefaultEpsilon").Get(gDefaultEpsilon);
      c("SysDirs").Get(tempString);
      if (!tempString.empty()) {
        std::stringstream ss (tempString);
        ss.imbue(std::locale(std::locale(), new tokens));
        dmSysDirs.assign(
            (std::istream_iterator<std::string>(ss)),
            std::istream_iterator<std::string>());
      }
    } else if (name == "InputFiles") {
      c("MassComboFile").Get(massComboFileName);
      c("DMEfficiencySmear").Get(dmEffSmearFile);
      c("DMEfficiencyPlotName").Get(dmEffHistName);
      c("DMSmearPlotName").Get(dmSmearHistName);
      c("DataFile").Get(dataFileName);
      c("DataHistName").Get(dataHistName);
      c("ExcessHistName").Get(excessHistName);
      c("BackgroundHistName").Get(bkgHistName);
      c("ApplyDetEff").Get(applyEfficiency);
      c("ResultsFromFile").Get(resultsFromFile);
      if (resultsFromFile > 0) {
        c("PreviousResults").Get(tempString);
        std::stringstream ss (tempString);
        ss.imbue(std::locale(std::locale(), new tokens));
        previousResults.assign(
            (std::istream_iterator<std::string>(ss)),
            std::istream_iterator<std::string>());

        c("PreviousResultsNames").Get(tempString);
        std::stringstream ss2(tempString);
        ss2.imbue(std::locale(std::locale(), new tokens));
        previousResultsNames.assign(
            (std::istream_iterator<std::string>(ss2)),
            std::istream_iterator<std::string>());
      } // end if resultsFromFile
    } else if (name == "POTAndScaling") {
      c("DefaultPOT").Get(defaultPOT);
      c("POTScale").Get(potScaling);
      c("BkgScale").Get(scaleDataAndBkg);
      c("NormalizeBkg").Get(normalizeBkgValue);
    } else if (name == "FitParams") {
      c("DoFit").Get(doFit);
      c("StartEnergy").Get(startEnergy);
      c("EndEnergy").Get(endEnergy);
      int temp = 0;
      c("NumberNullThrows").Get(temp);
      gNumNullFakeData = temp;
      c("NumberSigThrows").Get(temp);
      gNumSigFakeData = temp;
    } else if (name == "OutputFiles") {
      c("Directory").Get(outputDir);
      c("PostFix").Get(postFix);
      c("OtherLimitsDir").Get(otherLimitsDir);
      c("MedPartMassRatio").Get(mvmChiRatio);
      c("OutputAlphaD").Get(newAlphaD);
    } else if (name == "Systematics") {
      std::string file = "";
      std::string histName = "";
      double fracError = 0;
      c("SysFileOM").Get(file);
      if (!file.empty()) {
        c("SysHistOM").Get(histName);
        fracError = 0;
      } else {
        c("SysNumberOM").Get(fracError);
        histName = "";
      }
      inputOMSysInfo = std::tuple<std::string,std::string,double>(file,histName,fracError);

      c("SysFilePOT").Get(file);
      if (!file.empty()) {
        c("SysHistPOT").Get(histName);
        fracError = 0;
      } else {
        c("SysNumberPOT").Get(fracError);
        histName = "";
      }
      inputPOTSysInfo = std::tuple<std::string,std::string,double>(file,histName,fracError);

      c("SysFileQF").Get(file);
      if (!file.empty()) {
        c("SysHistQF").Get(histName);
        fracError = 0;
      } else {
        c("SysNumberQF").Get(fracError);
        histName = "";
      }
      inputQFSysInfo = std::tuple<std::string,std::string,double>(file,histName,fracError);
    } // end if-else module names
  } // end for mList

  MsgInfo(MsgLog::Form("DMTrueHists: Master Dir = %s",dmMasterDir.c_str()));
  MsgInfo(MsgLog::Form("DMTrueHists: Dir = %s",dmDir.c_str()));
  MsgInfo(MsgLog::Form("DMTrueHists: Hist Name = %s",dmHistName.c_str()));
  MsgInfo(MsgLog::Form("DMTrueHists: Hist Axis = %d",dmHistAxis));
  MsgInfo(MsgLog::Form("DMTrueHists: Default Epsilon %.f",gDefaultEpsilon));
  MsgInfo(MsgLog::Form("DMTrueHists: Default Alpha D%.f",gDefaultAlphaD));
  for (auto & dir : dmSysDirs) {
    MsgInfo(MsgLog::Form("DMTrueHists: Adding Sys Dir = %s",dir.c_str()));
  }
  MsgInfo(MsgLog::Form("InputFiles: Mass Combo File %s",massComboFileName.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Take results from previous files %s", (resultsFromFile > 0) ? "true" : "false"));
  for (size_t i = 0; i < previousResults.size(); ++i) {
    MsgInfo(MsgLog::Form("InputFiles: Prev Result %s located %s",previousResultsNames.at(i).c_str(),previousResults.at(i).c_str()));
  }
  MsgInfo(MsgLog::Form("InputFiles: Apply Detector Efficiency %s (smearing will always be applied)", (applyEfficiency > 0) ? "true" : "false"));
  MsgInfo(MsgLog::Form("InputFiles: Dark Matter Efficiency and Smear Hist File %s",dmEffSmearFile.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Dark Matter Efficiency Hist Name %s",dmEffHistName.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Dark Matter Smear Hist Name %s",dmSmearHistName.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Data Hist File %s",dataFileName.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Data Hist Name %s",dataHistName.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Excess Hist Name %s",excessHistName.c_str()));
  MsgInfo(MsgLog::Form("InputFiles: Background Hist Name %s",bkgHistName.c_str()));
  MsgInfo(MsgLog::Form("POTAndScaling: Default POT %g",defaultPOT));
  MsgInfo(MsgLog::Form("POTAndScaling: POT Scale %g",potScaling));
  MsgInfo(MsgLog::Form("POTAndScaling: Background Scale %g",scaleDataAndBkg));
  MsgInfo(MsgLog::Form("POTAndScaling: Normalize Background Value %g",normalizeBkgValue));
  MsgInfo(MsgLog::Form("FitParams: Do Fit %s",(doFit > 0) ? "true" : "false"));
  MsgInfo(MsgLog::Form("FitParams: Start Energy %g",startEnergy));
  MsgInfo(MsgLog::Form("FitParams: End Energy %g",endEnergy));
  MsgInfo(MsgLog::Form("FitParams: Number Null Fake Data Throws %zu",gNumNullFakeData));
  MsgInfo(MsgLog::Form("FitParams: Number Sig Fake Data Throws %zu",gNumSigFakeData));
  MsgInfo(MsgLog::Form("OutputFiles: Output Directory %s",outputDir.c_str()));
  MsgInfo(MsgLog::Form("OutputFiles: Output Postfix %s",postFix.c_str()));
  MsgInfo(MsgLog::Form("OutputFiles: Other Limtis Dir %s",otherLimitsDir.c_str()));
  MsgInfo(MsgLog::Form("OutputFiles: Output ALphaD %f",newAlphaD));
  MsgInfo(MsgLog::Form("OutputFiles: Mediator/Particle Mass Ratio %f",mvmChiRatio));
  MsgInfo(MsgLog::Form("Systematics: POT Info file (%s) hist (%s) fracError (%f)",
        std::get<0>(inputPOTSysInfo).c_str(),std::get<1>(inputPOTSysInfo).c_str(),std::get<2>(inputPOTSysInfo)));
  MsgInfo(MsgLog::Form("Systematics: QF Info file (%s) hist (%s) fracError (%f)",
        std::get<0>(inputQFSysInfo).c_str(),std::get<1>(inputQFSysInfo).c_str(),std::get<2>(inputQFSysInfo)));
  MsgInfo(MsgLog::Form("Systematics: OM Info file (%s) hist (%s) fracError (%f)",
        std::get<0>(inputOMSysInfo).c_str(),std::get<1>(inputOMSysInfo).c_str(),std::get<2>(inputOMSysInfo)));

  //////////////////////////////////////////////////////////
  // Now that everything is set up we start 
  // getting all the necessary histograms/matrices and 
  // do the fit
  //////////////////////////////////////////////////////////

  TH2D * energyTrueVsReco = nullptr;
  TH1D * efficiencyTrueEnergy = nullptr;

  std::shared_ptr<TFile> rootFile = std::make_shared<TFile>(dmEffSmearFile.c_str(),"READ");
  rootFile->GetObject(dmSmearHistName.c_str(),energyTrueVsReco);
  if (!energyTrueVsReco) {
    MsgFatal(MsgLog::Form("Could not get smearing hist %s from %s",dmSmearHistName.c_str(),dmEffSmearFile.c_str()));
  }

  rootFile->GetObject(dmEffHistName.c_str(),efficiencyTrueEnergy);
  if (!efficiencyTrueEnergy) {
    MsgFatal(MsgLog::Form("Could not get efficiency hist %s from %s",dmEffHistName.c_str(),dmEffSmearFile.c_str()));
  }

  std::map<double,std::shared_ptr<TH1D>> mapBkgPlusSig;
  std::shared_ptr<TH1D> dataHist = std::shared_ptr<TH1D>(energyTrueVsReco->ProjectionY("dataHist"));
  std::shared_ptr<TH1D> bkgHist = std::shared_ptr<TH1D>(energyTrueVsReco->ProjectionY("bkgHist"));
  std::shared_ptr<TH1D> excessHist = std::shared_ptr<TH1D>(energyTrueVsReco->ProjectionY("excessHist"));
  dataHist->Reset("ICESM");
  bkgHist->Reset("ICESM");
  excessHist->Reset("ICESM");
  TH1D * tempHist = 0;
  rootFile = std::make_shared<TFile>(dataFileName.c_str(),"READ");
  rootFile->GetObject(dataHistName.c_str(),tempHist);
  if (!tempHist) {
    MsgFatal(MsgLog::Form("Could not get %s from %s",dataHistName.c_str(),dataFileName.c_str()));
  }
  for (int binX = 0; binX <= dataHist->GetNbinsX()+1; ++binX) {
    dataHist->SetBinContent(binX,tempHist->GetBinContent(binX)*scaleDataAndBkg*potScaling);
    dataHist->SetBinError(binX,tempHist->GetBinError(binX)*scaleDataAndBkg*potScaling);
  }
  delete tempHist;
  tempHist = nullptr;

  rootFile->GetObject(excessHistName.c_str(),tempHist);
  if (!tempHist) {
    MsgFatal(MsgLog::Form("Could not get %s from %s",excessHistName.c_str(),dataFileName.c_str()));
  }
  for (int binX = 0; binX <= excessHist->GetNbinsX()+1; ++binX) {
    excessHist->SetBinContent(binX,tempHist->GetBinContent(binX));
    excessHist->SetBinError(binX,tempHist->GetBinError(binX));
  }
  delete tempHist;
  tempHist = nullptr;

  rootFile->GetObject(bkgHistName.c_str(),tempHist);
  if (!tempHist) {
    MsgFatal(MsgLog::Form("Could not get %s from %s",bkgHistName.c_str(),dataFileName.c_str()));
  }
  for (int binX = 0; binX <= bkgHist->GetNbinsX()+1; ++binX) {
    bkgHist->SetBinContent(binX,tempHist->GetBinContent(binX)*scaleDataAndBkg*potScaling);
    bkgHist->SetBinError(binX,tempHist->GetBinError(binX)*scaleDataAndBkg*potScaling);
  }
  MsgDebug(2,MsgLog::Form("Bkg from file Integral %f After copy %f",tempHist->Integral(),bkgHist->Integral()));
  delete tempHist;
  tempHist = nullptr;

  gStartBin = bkgHist->FindBin(startEnergy);
  gEndBin = bkgHist->FindBin(endEnergy);
  gNumRecoBins = dataHist->GetNbinsX();
  MsgDebug(2,MsgLog::Form("Background Integral before normalization %.3f",std::sqrt(bkgHist->Integral(gStartBin,gEndBin))));
  if (normalizeBkgValue >= 0) {
    double bkgScale = normalizeBkgValue/bkgHist->Integral(gStartBin,gEndBin);
    bkgHist->Scale(bkgScale);
    dataHist->Scale(bkgScale);
  }
  MsgDebug(2,MsgLog::Form("Background Integral after normalization %.3f",std::sqrt(bkgHist->Integral(gStartBin,gEndBin))));

  TMatrixD dataMatrixBefore(gEndBin-gStartBin,1);
  TMatrixD dataMatrix(gEndBin-gStartBin,1);
  TMatrixD bkgMatrix(gEndBin-gStartBin,1);
  TMatrixD bkgErrorMatrix(gEndBin-gStartBin,gEndBin-gStartBin);
  TMatrixD dmMatrix(gEndBin-gStartBin,1);

  TH2D * qf_frac_error = nullptr;
  TH2D * om_frac_error = nullptr;
  TH2D * pot_frac_error = nullptr;
  if (!std::get<0>(inputQFSysInfo).empty()) {
    rootFile = std::make_shared<TFile>(std::get<0>(inputQFSysInfo).c_str(),"READ");
    if (!rootFile) {
      MsgFatal("QF fractional error file given but could not open");
    }
    rootFile->GetObject(std::get<1>(inputQFSysInfo).c_str(),qf_frac_error);
  }
  if (!std::get<0>(inputOMSysInfo).empty()) {
    rootFile = std::make_shared<TFile>(std::get<0>(inputOMSysInfo).c_str(),"READ");
    if (!rootFile) {
      MsgFatal("OM fractional error file given but could not open");
    }
    rootFile->GetObject(std::get<1>(inputOMSysInfo).c_str(),om_frac_error);
  }
  if (!std::get<0>(inputPOTSysInfo).empty()) {
    rootFile = std::make_shared<TFile>(std::get<0>(inputPOTSysInfo).c_str(),"READ");
    if (!rootFile) {
      MsgFatal("POT fractional error file given but could not open");
    }
    rootFile->GetObject(std::get<1>(inputPOTSysInfo).c_str(),pot_frac_error);
  }


  std::vector<TMatrixD> sysFracErrorMatrix;
  for (int i=0; i < 4; ++i) {
    sysFracErrorMatrix.emplace_back(TMatrixD(gEndBin-gStartBin,gEndBin-gStartBin));
    for (int bin = 0; bin < gEndBin-gStartBin; ++bin) {
      switch(i) {
        case 0: sysFracErrorMatrix.back()(bin,bin) = std::get<2>(inputOMSysInfo); break;//0.5; break;
        case 1: sysFracErrorMatrix.back()(bin,bin) = std::get<2>(inputPOTSysInfo); break;
        case 2: sysFracErrorMatrix.back()(bin,bin) = 0; break;// this is the place holder for any systematic that requires multiple BdNMC runs
        case 3: sysFracErrorMatrix.back()(bin,bin) = std::get<2>(inputQFSysInfo); break;//0.25; break;
      }
      for (int bin2 = 0; bin2 < gEndBin - gStartBin; ++bin2) {
        if (bin2 == bin) {
          continue;
        }
        sysFracErrorMatrix.back()(bin,bin2) = 0;
        sysFracErrorMatrix.back()(bin2,bin) = 0;
        switch(i) {
          case 1: sysFracErrorMatrix.back()(bin,bin2) = std::pow(std::get<2>(inputPOTSysInfo),2.0);
                  sysFracErrorMatrix.back()(bin2,bin) = sysFracErrorMatrix.back()(bin,bin2);
                  break;
        //  case 2: sysFracErrorMatrix.back()(bin,bin2) = 0.0049; // holder for systematic that requires multiple BdNMC runs
        //          sysFracErrorMatrix.back()(bin2,bin) = 0.0049; // holder for systematic that requires multiple BdNMC runs
        //          break;
          default: break;
        }
      }
    }
  }

  MsgDebug(2,MsgLog::Form("Bkg Hist before bin width %f",bkgHist->Integral()));
  for (int binX = 1; binX < gNumRecoBins+1; ++binX) {
    double binWidth = bkgHist->GetXaxis()->GetBinWidth(binX);
    dataHist->SetBinContent(binX,dataHist->GetBinContent(binX)*binWidth);
    bkgHist->SetBinContent(binX,bkgHist->GetBinContent(binX)*binWidth);
    bkgHist->SetBinError(binX,std::sqrt(bkgHist->GetBinContent(binX)));
    if (binX >= gStartBin && binX < gEndBin) {
      dataMatrixBefore(binX-gStartBin,0) = dataHist->GetBinContent(binX);
      dataMatrix(binX-gStartBin,0) = dataHist->GetBinContent(binX);
      bkgMatrix(binX-gStartBin,0) = bkgHist->GetBinContent(binX);
      bkgErrorMatrix(binX-gStartBin,binX-gStartBin) = bkgHist->GetBinContent(binX);
      if (bkgHist->GetBinContent(binX) < 1e-4) {
        bkgErrorMatrix(binX-gStartBin,binX-gStartBin) = 1e-4;
      } else {
        bkgErrorMatrix(binX-gStartBin,binX-gStartBin) = bkgHist->GetBinContent(binX);
      }
      for (int binY = 1; binY < gNumRecoBins+1; ++binY) {
        if (binY >= gStartBin && binY < gEndBin) {
          if (qf_frac_error != nullptr) {
            sysFracErrorMatrix.at(3)(binX-gStartBin,binY-gStartBin) = qf_frac_error->GetBinContent(binX,binY);
          }
          if (om_frac_error != nullptr) {
            sysFracErrorMatrix.at(0)(binX-gStartBin,binY-gStartBin) = om_frac_error->GetBinContent(binX,binY);
          }
          if (pot_frac_error != nullptr) {
            sysFracErrorMatrix.at(0)(binX-gStartBin,binY-gStartBin) = pot_frac_error->GetBinContent(binX,binY);
          }
        }
      }
    }
  }
  if (qf_frac_error) {
    delete qf_frac_error;
  }
  if (om_frac_error) {
    delete om_frac_error;
  }
  if (pot_frac_error) {
    delete pot_frac_error;
  }
  MsgDebug(2,MsgLog::Form("Bkg Hist after bin width %f",bkgHist->Integral()));

  //sysFracErrorMatrix.at(3).Print();
  //abort();

  auto rebinHist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(efficiencyTrueEnergy->Clone("rebinHist")));
  auto rebinHistReco = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(energyTrueVsReco->ProjectionY("rebinHistReco")));
  rebinHist->Reset("ICESM");
  rebinHistReco->Reset("ICESM");

  std::vector<std::shared_ptr<TH1D>> sysRebinHist;
  std::vector<std::shared_ptr<TH1D>> sysRebinHistReco;

  float mv = 0;
  float mchi = 0;

  double scaling = 0;
  std::vector<std::pair<double,double>> mVmChiVec;
  std::ifstream infile(Form("%sminiboone_full_nucleon_timing_vector_portal_cl_epsilon4alphaD.txt",otherLimitsDir.c_str()));
  while (infile >> mv >> mchi >> scaling) {
    //if (mv <= mvValues.back() && mchi < mchiValues.back()) {
      mVmChiVec.emplace_back(std::make_pair(mv,mchi));
    //}
  }
  infile.close();

  TLegend * leg = new TLegend(0.5,0.7,0.94,0.94);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(102);
  leg->SetHeader("m_{V} = 3m_{#chi} MeV");
  leg->SetNColumns(2);

  TGraphErrors * gr = new TGraphErrors();
  TGraph * grY = new TGraph();
  grY->SetLineColor(gkBlack);
  grY->SetLineStyle(1);
  grY->SetLineWidth(2);

  TGraph * grY44ns = new TGraph();
  grY44ns->SetLineColor(kBlack);
  grY44ns->SetLineStyle(2);
  grY44ns->SetLineWidth(4);

  gStyle->SetPalette(kBird);
  auto palette = TColor::GetPalette();
  gStyle->SetTitleFont(102,"t");
  gStyle->SetTitleFont(102,"x");
  gStyle->SetTitleFont(102,"y");
  gStyle->SetTitleFont(102,"z");
  gStyle->SetLabelFont(102,"x");
  gStyle->SetLabelFont(102,"y");
  gStyle->SetLabelFont(102,"z");

  std::vector<std::shared_ptr<TH1D>> keep;
  std::vector<std::shared_ptr<TH1D>> keepReco;
  std::vector<std::shared_ptr<TH1D>> keepRecoAll;


  std::vector<double> totalValue;
  std::vector<double> totalRMS;
  std::vector<double> totalMaxBin;

  int count = 0;
  std::vector<double> mchiValues;
  std::vector<double> mvValues;

  TGraph2D * trueNumEvents2D = new TGraph2D();
  trueNumEvents2D->SetName("trueNumEvents2D");
  TGraph2D * recoNumEvents2D = new TGraph2D();
  recoNumEvents2D->SetName("recoNumEvents2D");
  TGraph2D * trueMeanEnergy2D = new TGraph2D();
  trueMeanEnergy2D->SetName("trueMeanEnergy2D");
  TGraph2D * recoMeanEnergy2D = new TGraph2D();
  recoMeanEnergy2D->SetName("recoMeanEnergy2D");

  std::vector<std::shared_ptr<TH1D>> sysHistVec;

  TGraph2D * ySens2DGr = new TGraph2D();
  ySens2DGr->SetName("ySens2DGr");
  TGraph2D * yCL2DGr = new TGraph2D();
  yCL2DGr->SetName("yCL2DGr");

  TGraph2D * ySens2DGr200 = new TGraph2D();
  ySens2DGr200->SetName("ySens2DGr200");
  TGraph2D * yCL2DGr200 = new TGraph2D();
  yCL2DGr200->SetName("yCL2DGr200");

  TGraph2D * ySens2DGr200Low = new TGraph2D();
  ySens2DGr200Low->SetName("ySens2DGr200Low");
  TGraph2D * yCL2DGr200Low = new TGraph2D();
  yCL2DGr200Low->SetName("yCL2DGr200Low");

  TCanvas * canvas = new TCanvas("canvas","canvas",500,500);
  gStyle->SetOptStat(0);


  infile.open(massComboFileName);
  while (infile >> mv >> mchi) {
    //mv *= 1e3;
    //mchi *= 1e3;
    if (mv > 136) {
      continue;
    }

    gResults.emplace_back(mv,mchi,0,0);
    if (std::find(mchiValues.begin(),mchiValues.end(),mchi/1e3) == mchiValues.end()) {
      mchiValues.emplace_back(mchi/1e3);
    }
    if (std::find(mvValues.begin(),mvValues.end(),mv/1e3) == mvValues.end()) {
      mvValues.emplace_back(mv/1e3);
    }
  }
  infile.close();

  std::sort(mchiValues.begin(),mchiValues.end());
  std::sort(mvValues.begin(),mvValues.end());
  MsgInfo(MsgLog::Form("mchiValues: front %f back %f",mchiValues.front(),mchiValues.back()))
  MsgInfo(MsgLog::Form("mvValues: front %f back %f",mvValues.front(),mvValues.back()))

  if (resultsFromFile <= 0) {
    const auto kMaxNumberThreads = std::thread::hardware_concurrency();
    MsgInfo(MsgLog::Form("Maximum Number of threads = %zu",kMaxNumberThreads));
    std::vector<std::pair<size_t,std::future<void>>> futureVector;
    const size_t kNumMassCombos = gResults.size();
    for (size_t result = 0; result < kNumMassCombos; ++result) {
      gMTX.lock();
      mv = std::get<0>(gResults.at(result));
      mchi = std::get<1>(gResults.at(result));
      gMTX.unlock();
      MsgInfo(MsgLog::Form("Looking for mv %.1f mchi %.1f",mv,mchi));
      std::shared_ptr<TH1D> total = GetTotalHist(mv,mchi,dmHistAxis,dmHistName,dmMasterDir,dmDir,dmSysDirs,&sysHistVec);
      if (total == nullptr) {
        continue;
      }
      if (sysHistVec.size() > sysRebinHist.size()) {
        sysRebinHist.resize(sysHistVec.size());
        sysRebinHistReco.resize(sysHistVec.size());
        size_t count = 0;
        for (auto & hist : sysRebinHist) {
          if (hist != nullptr) {
            hist->Reset("ICESM");
          } else {
            hist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(rebinHist->Clone(Form("sysRebinHist_%zu",count))));
          }
          ++count;
        }
        count = 0;
        for (auto & hist : sysRebinHistReco) {
          if (hist != nullptr) {
            hist->Reset("ICESM");
          } else {
            hist = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(rebinHistReco->Clone(Form("sysRebinHistReco_%zu",count))));
          }
          ++count;
        }
      } else {
        for (auto & hist : sysRebinHist) {
          hist->Reset("ICESM");
        }
        for (auto & hist : sysRebinHistReco) {
          hist->Reset("ICESM");
        }
      }
      MsgDebug(2,MsgLog::Form("Got hist %s with Integral %f",total->GetName(),total->Integral()));
      rebinHist->Reset("ICESM");
      for (int i=0; i < total->GetNbinsX(); ++i) {
        double binCenter = total->GetXaxis()->GetBinLowEdge(i+1);
        double content = total->GetBinContent(i+1);
        int bin = rebinHist->FindBin(binCenter);
        rebinHist->AddBinContent(bin,content);
        for (size_t sys = 0; sys < sysHistVec.size(); ++sys) {
          sysRebinHist.at(sys)->AddBinContent(bin,sysHistVec.at(sys)->GetBinContent(i+1));
        }
      }
      keep.emplace_back(dynamic_cast<TH1D*>(rebinHist->Clone(Form("%s_keep_%.1f_%.1f",total->GetName(),mv,mchi))));
      keep.back()->Scale(defaultPOT*potScaling);
      for (auto & hist : sysRebinHist) {
        hist->Scale(defaultPOT*potScaling);
      }
      trueNumEvents2D->SetPoint(trueNumEvents2D->GetN(),mv/1e3,mchi/1e3,keep.back()->Integral());
      if (applyEfficiency > 0) {
        keep.back()->Multiply(efficiencyTrueEnergy);
        for (auto & hist : sysRebinHist) {
          hist->Multiply(efficiencyTrueEnergy);
        }
      }
      trueMeanEnergy2D->SetPoint(trueMeanEnergy2D->GetN(),mv/1e3,mchi/1e3,keep.back()->GetMean());
      keepReco.emplace_back(dynamic_cast<TH1D*>(rebinHistReco->Clone(Form("%s_keepReco_%.1f_%.1f",total->GetName(),mv,mchi))));
      keepReco.back()->Reset("ICESM");
      for (int binX = 1; binX < keep.back()->GetNbinsX()+1; ++binX) {
        double events = keep.back()->GetBinContent(binX);
        if (events == 0) {
          continue;
        }
        double sum = 0;
        for (int binY = 1; binY < energyTrueVsReco->GetNbinsY()+1; ++binY) {
          double scale = energyTrueVsReco->GetBinContent(binX,binY);
          if (scale == 0) {
            continue;
          }
          sum += scale;
          keepReco.back()->SetBinContent(binY,
              keepReco.back()->GetBinContent(binY) + scale*events);
        }
      }
      for (size_t sys = 0; sys < sysHistVec.size(); ++sys) {
        for (int binX = 1; binX < sysRebinHist.at(sys)->GetNbinsX()+1; ++binX) {
          double events = sysRebinHist.at(sys)->GetBinContent(binX);
          if (events == 0) {
            continue;
          }
          double sum = 0;
          for (int binY = 1; binY < energyTrueVsReco->GetNbinsY()+1; ++binY) {
            double scale = energyTrueVsReco->GetBinContent(binX,binY);
            if (scale == 0) {
              continue;
            }
            sum += scale;
            sysRebinHistReco.at(sys)->SetBinContent(binY,
                sysRebinHistReco.at(sys)->GetBinContent(binY) + scale*events);
          }
        }
      }
      ++count;
      if (count > 256) {
        count = 1;
      }
      keep.back()->SetLineColor(palette.At(count));
      keepReco.back()->SetLineColor(palette.At(count));
      recoNumEvents2D->SetPoint(recoNumEvents2D->GetN(),mv/1e3,mchi/1e3,keepReco.back()->Integral());
      recoMeanEnergy2D->SetPoint(recoMeanEnergy2D->GetN(),mv/1e3,mchi/1e3,keepReco.back()->GetMean());
      keep.back()->Scale(1.0,"width");
      //keepReco.back()->Scale(1.0/keepReco.back()->Integral(),"width");
      //

      if (sysRebinHistReco.size() > 0) {
        for (int binX = gStartBin; binX < gEndBin; ++binX) {
          double contentX = keepReco.back()->GetBinContent(binX);
          for (int binY = gStartBin; binY < gEndBin; ++binY) {
            double contentY = keepReco.back()->GetBinContent(binY);
            double sum = 0;
            for (auto & hist : sysRebinHistReco) {
              sum += (hist->GetBinContent(binX) - contentX)*(hist->GetBinContent(binY) - contentY);
            }
            sum /= static_cast<double>(sysRebinHistReco.size());
            sysFracErrorMatrix.at(2)(binX-gStartBin,binY-gStartBin) = sum/contentX/contentY;
          }
        }
      }

      if (doFit <= 0) {
        continue;
      }

      for (int bin = gStartBin; bin < gEndBin; ++bin) {
        dmMatrix(bin-gStartBin,0) = keepReco.back()->GetBinContent(bin);
      }

      if (futureVector.size() < kMaxNumberThreads-2) {
        futureVector.emplace_back(std::make_pair(result,std::async(CalculateSensCL,dmMatrix,std::ref(dataMatrix),std::ref(bkgMatrix),
              std::ref(sysFracErrorMatrix), std::ref(bkgErrorMatrix),std::ref(result))));
      } else {
        while (futureVector.size() == kMaxNumberThreads-2) {
          std::this_thread::sleep_for(std::chrono::seconds(10));
          int finished = 0;
          for (auto it = futureVector.begin(); it != futureVector.end();) {
            bool status = it->second.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
            if (status) {
              gMTX.lock();
              MsgInfo(MsgLog::Form("Finished and Saved mv %.1f mchi %.1f sens %g cl %g",
                    std::get<0>(gResults.at(it->first)), std::get<1>(gResults.at(it->first)),
                    std::get<2>(gResults.at(it->first)), std::get<3>(gResults.at(it->first))));
              gMTX.unlock();
              it = futureVector.erase(it);
              ++finished;
            } else {
              gMTX.lock();
              MsgInfo(MsgLog::Form("mV %.1f and mChi %.1f still running sens %g cl %g",
                    std::get<0>(gResults.at(it->first)), std::get<1>(gResults.at(it->first)),
                    std::get<2>(gResults.at(it->first)), std::get<3>(gResults.at(it->first))));
              gMTX.unlock();
              std::advance(it,1);
            }
          }
          MsgInfo(MsgLog::Form("%d finished since last check",finished));
        }
      }
    } // end for gResults

    //now we wait for all the futures to finish
    //once everyone is finished then we can plot the results
    while (!futureVector.empty()) {
      std::this_thread::sleep_for(std::chrono::seconds(10));
      int finished = 0;
      for (auto it = futureVector.begin(); it != futureVector.end();) {
        bool status = it->second.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
        if (status) {
          gMTX.lock();
          MsgInfo(MsgLog::Form("Finished and Saved mv %.1f mchi %.1f sens %g cl %g",
                std::get<0>(gResults.at(it->first)), std::get<1>(gResults.at(it->first)),
                std::get<2>(gResults.at(it->first)), std::get<3>(gResults.at(it->first))));
          gMTX.unlock();
          it = futureVector.erase(it);
          ++finished;
        } else {
          gMTX.lock();
          MsgInfo(MsgLog::Form("mV %.1f and mChi %.1f still running sens %g cl %g",
                std::get<0>(gResults.at(it->first)), std::get<1>(gResults.at(it->first)),
                std::get<2>(gResults.at(it->first)), std::get<3>(gResults.at(it->first))));
          gMTX.unlock();
          std::advance(it,1);
        }
      }
      MsgInfo(MsgLog::Form("%d finished since last check %zu still running",finished,futureVector.size()));
    }

    for (auto & hist : keepReco) {
      hist->Scale(1.0,"width");
    }

    TH2D * ySens2D = new TH2D("ySens2D",";m_{V} (GeV/c^{2});m_{#chi} (GeV/c^{2});90% #varepsilon^{4}#alpha_{D} Sens.",
        mvValues.size()-1,&mvValues.front(),mchiValues.size()-1,&mchiValues.front());
    TH2D * yCL2D = new TH2D("yCL2D",";m_{V} (GeV/c^{2});m_{#chi} (GeV/c^{2});90% #varepsilon^{4}#alpha_{D} CL",
        mvValues.size()-1,&mvValues.front(),mchiValues.size()-1,&mchiValues.front());
    TGraph2D * ySens2DGr1 = new TGraph2D();
    ySens2DGr1->SetName("ySens2DGr1");
    TGraph2D * yCL2DGr1 = new TGraph2D();
    yCL2DGr1->SetName("yCL2DGr1");

    double minSens2D = 9e9;
    double minCL2D = 9e9;
    double maxSens2D = 0;
    double maxCL2D = 0;
    for (auto & result : gResults) {
      mv = std::get<0>(result);
      mchi = std::get<1>(result);
      MsgInfo(MsgLog::Form("mv %.1f mchi %.1f sens %g cl %g",mv,mchi,std::get<2>(result),std::get<3>(result)));
      if (std::get<2>(result) > 0 && std::get<3>(result) > 0) {
        ySens2DGr1->SetPoint(ySens2DGr1->GetN(),mv/1e3,mchi/1e3,std::get<2>(result));
        yCL2DGr1->SetPoint(yCL2DGr1->GetN(),mv/1e3,mchi/1e3,std::get<3>(result));
        if (std::get<2>(result) != 0) {
          minSens2D = std::min(std::get<2>(result),minSens2D);
          maxSens2D = std::max(std::get<2>(result),maxSens2D);
        }
        if (std::get<3>(result) != 0) {
          minCL2D = std::min(std::get<3>(result),minCL2D);
          maxCL2D = std::max(std::get<3>(result),maxCL2D);
        }
      }
    }

    ySens2D->Print();
    ySens2D->Smooth();
    yCL2D->Smooth();

    for (auto & p : mVmChiVec) {
      ySens2DGr->SetPoint(ySens2DGr->GetN(),p.first,p.second,ySens2DGr1->Interpolate(p.first,p.second));
      yCL2DGr->SetPoint(yCL2DGr->GetN(),p.first,p.second,yCL2DGr1->Interpolate(p.first,p.second));
    }

    canvas->Divide(2,2);
    canvas->cd(1);
    gPad->SetMargin(0.12,0.2,0.1,0.05);
    trueNumEvents2D->Draw("colz");
    trueNumEvents2D->GetHistogram()->GetXaxis()->SetTitleFont(102);
    trueNumEvents2D->GetHistogram()->GetYaxis()->SetTitleFont(102);
    trueNumEvents2D->GetHistogram()->GetXaxis()->SetLabelFont(102);
    trueNumEvents2D->GetHistogram()->GetYaxis()->SetLabelFont(102);
    trueNumEvents2D->GetXaxis()->SetTitle(ySens2D->GetXaxis()->GetTitle());
    trueNumEvents2D->GetYaxis()->SetTitle(ySens2D->GetYaxis()->GetTitle());
    trueNumEvents2D->GetZaxis()->SetTitle("Number Events 100% Eff");
    trueNumEvents2D->SetMinimum(1e-1);
    gPad->SetLogy(true);
    gPad->SetLogx(true);
    gPad->SetLogz(true);
    canvas->cd(2);
    gPad->SetMargin(0.12,0.2,0.1,0.05);
    recoNumEvents2D->Draw("colz");
    recoNumEvents2D->GetHistogram()->GetXaxis()->SetTitleFont(102);
    recoNumEvents2D->GetHistogram()->GetYaxis()->SetTitleFont(102);
    recoNumEvents2D->GetHistogram()->GetXaxis()->SetLabelFont(102);
    recoNumEvents2D->GetHistogram()->GetYaxis()->SetLabelFont(102);
    recoNumEvents2D->GetXaxis()->SetTitle(ySens2D->GetXaxis()->GetTitle());
    recoNumEvents2D->GetYaxis()->SetTitle(ySens2D->GetYaxis()->GetTitle());
    recoNumEvents2D->GetZaxis()->SetTitle("Reco. Number Events");
    recoNumEvents2D->SetMinimum(1e-3);
    gPad->SetLogy(true);
    gPad->SetLogx(true);
    gPad->SetLogz(true);
    canvas->cd(3);
    gPad->SetMargin(0.12,0.2,0.1,0.05);
    trueMeanEnergy2D->Draw("colz");
    trueMeanEnergy2D->GetHistogram()->GetXaxis()->SetTitleFont(102);
    trueMeanEnergy2D->GetHistogram()->GetYaxis()->SetTitleFont(102);
    trueMeanEnergy2D->GetHistogram()->GetXaxis()->SetLabelFont(102);
    trueMeanEnergy2D->GetHistogram()->GetYaxis()->SetLabelFont(102);
    trueMeanEnergy2D->GetXaxis()->SetTitle(ySens2D->GetXaxis()->GetTitle());
    trueMeanEnergy2D->GetYaxis()->SetTitle(ySens2D->GetYaxis()->GetTitle());
    trueMeanEnergy2D->GetZaxis()->SetTitle("Mean True Energy after Reco. (keV)");
    trueMeanEnergy2D->SetMinimum(50);
    trueMeanEnergy2D->SetMaximum(1000);
    gPad->SetLogy(true);
    gPad->SetLogx(true);
    gPad->SetLogz(true);
    canvas->cd(4);
    gPad->SetMargin(0.12,0.2,0.1,0.05);
    recoMeanEnergy2D->Draw("colz");
    recoMeanEnergy2D->GetHistogram()->GetXaxis()->SetTitleFont(102);
    recoMeanEnergy2D->GetHistogram()->GetYaxis()->SetTitleFont(102);
    recoMeanEnergy2D->GetHistogram()->GetXaxis()->SetLabelFont(102);
    recoMeanEnergy2D->GetHistogram()->GetYaxis()->SetLabelFont(102);
    recoMeanEnergy2D->GetXaxis()->SetTitle(ySens2D->GetXaxis()->GetTitle());
    recoMeanEnergy2D->GetYaxis()->SetTitle(ySens2D->GetYaxis()->GetTitle());
    recoMeanEnergy2D->GetZaxis()->SetTitle("Mean Reco Energy (PE)");
    recoMeanEnergy2D->SetMinimum(3);
    recoMeanEnergy2D->SetMaximum(30);
    gPad->SetLogy(true);
    gPad->SetLogx(true);
    //gPad->SetLogz(true);
    canvas->Update();
    canvas->Print(Form("%sdm_energy_numEvents_mv_mchi_ccm%s.png",outputDir.c_str(),postFix.c_str()),"png");

    canvas->SetMargin(0.12,0.05,0.1,0.05);
    canvas->SetLogx(true);
    gr->SetLineStyle(1);
    //gr->SetLineColor(kBlack);
    //gr->SetLineWidth(3);
    gr->SetFillColor(4);
    gr->SetFillStyle(3010);
    gr->Draw("a3");
    gr->GetHistogram()->GetXaxis()->SetTitleFont(102);
    gr->GetHistogram()->GetYaxis()->SetTitleFont(102);
    gr->GetHistogram()->GetXaxis()->SetLabelFont(102);
    gr->GetHistogram()->GetYaxis()->SetLabelFont(102);
    gr->GetHistogram()->GetXaxis()->SetTitle("m_{#chi} (GeV c^{-2})");
    gr->GetHistogram()->GetYaxis()->SetTitle("Energy bin with most \"events\" (MeV)");
    //gr->GetHistogram()->GetYaxis()->SetTitle("Bunch Time with most \"events\" (ns)");
    //gr->GetHistogram()->GetYaxis()->SetTitle("Time since #pi generation with most \"events\" (ns)");
    gr->Draw("AL");
    canvas->Print(Form("%smax_dm_energy_bin_mv_3mchi_ccm%s.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

    dataMatrixBefore.Print();
    dataMatrix.Print();

    canvas->Clear();
    canvas->Divide(2,1);
    canvas->cd(1);
    gPad->SetLogx(true);
    gPad->SetLogy(true);
    gPad->SetMargin(0.12,0.05,0.1,0.05);
    double maximum = 0;
    double minimum = 9e9;
    for (auto & hist : keep) {
      maximum = std::max(maximum,hist->GetMaximum());
      minimum = std::min(minimum,hist->GetMinimum(1e-15));
    }
    count = 0;
    for (auto & hist : keep) {
      if (count == 0) {
        hist->SetMaximum(maximum*1.05);
        hist->SetMinimum(minimum*0.95);
        hist->GetYaxis()->SetTitleOffset(1.4);
        hist->GetXaxis()->SetRangeUser(50,5000);
        hist->Draw("hist");
      } else {
        hist->Draw("same hist");
      }
      ++count;
    }
    canvas->cd(2);
    gPad->SetMargin(0.12,0.05,0.1,0.05);
    gPad->SetLogy(true);
    maximum = 0;
    for (auto & hist : keepReco) {
      maximum = std::max(maximum,hist->GetMaximum());
      minimum = std::min(minimum,hist->GetMinimum(1e-15));
    }
    count = 0;
    for (auto & hist : keepReco) {
      if (count == 0) {
        hist->SetMaximum(maximum*1.05);
        hist->SetMinimum(minimum*0.95);
        hist->GetYaxis()->SetTitleOffset(1.4);
        hist->GetXaxis()->SetRangeUser(0,30);
        hist->Draw("hist");
      } else {
        hist->Draw("same hist");
      }
      ++count;
    }
    canvas->cd();
    canvas->Update();
    canvas->Print(Form("%sdm_energy_true_dist_mv_3mchi_ccm%s.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

    canvas->Clear();
    gPad->SetMargin(0.12,0.18,0.1,0.05);
    gPad->SetLogz(true);
    ySens2DGr->Draw("colz");
    ySens2DGr->GetZaxis()->SetTitleOffset(2.0);
    ySens2DGr->SetMinimum(minSens2D*0.5);
    ySens2DGr->SetMaximum(maxSens2D*1.1);
    ySens2DGr->GetXaxis()->SetTitle(ySens2D->GetXaxis()->GetTitle());
    ySens2DGr->GetYaxis()->SetTitle(ySens2D->GetYaxis()->GetTitle());
    ySens2DGr->GetZaxis()->SetTitle(ySens2D->GetZaxis()->GetTitle());
    canvas->Print(Form("%sdm_ep4ad_sens_2d%s.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

    canvas->Clear();
    gPad->SetMargin(0.12,0.18,0.1,0.05);
    gPad->SetLogz(true);
    yCL2DGr->Draw("colz");
    yCL2DGr->GetZaxis()->SetTitleOffset(2.0);
    yCL2DGr->SetMinimum(minCL2D*0.5);
    yCL2DGr->SetMaximum(maxCL2D*1.1);
    yCL2DGr->GetXaxis()->SetTitle(yCL2D->GetXaxis()->GetTitle());
    yCL2DGr->GetYaxis()->SetTitle(yCL2D->GetYaxis()->GetTitle());
    yCL2DGr->GetZaxis()->SetTitle(yCL2D->GetZaxis()->GetTitle());
    canvas->Print(Form("%sdm_ep4ad_cl_2d%s.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

    auto outFile = TFile::Open(Form("%sdm_true_reco_dists%s.root",outputDir.c_str(),postFix.c_str()),"RECREATE");
    outFile->cd();
    for (auto & hist : keep) {
      hist->Write();
    }
    for (auto & hist : keepReco) {
      hist->Write();
    }
    ySens2DGr->Write();
    yCL2DGr->Write();
    delete outFile;

    return EXIT_SUCCESS;


  } else {
    const long int kPrevSize = previousResultsNames.size();
    auto loc = std::distance(previousResultsNames.begin(),std::find(previousResultsNames.begin(),previousResultsNames.end(),"CCM120"));
    if (loc != kPrevSize) {
      rootFile = std::make_shared<TFile>(previousResults.at(loc).c_str(),"READ");
      rootFile->GetObject("ySens2DGr",ySens2DGr);
      rootFile->GetObject("yCL2DGr",yCL2DGr);

      ySens2DGr->SetNpx(500);
      ySens2DGr->SetNpy(500);
      yCL2DGr->SetNpx(500);
      yCL2DGr->SetNpy(500);
    }

    auto loc2 = std::distance(previousResultsNames.begin(),std::find(previousResultsNames.begin(),previousResultsNames.end(),"CCM200"));
    if (loc2 != kPrevSize) {
      rootFile = std::make_shared<TFile>(previousResults.at(loc2).c_str(),"READ");
      rootFile->GetObject("ySens2DGr",ySens2DGr200);
      rootFile->GetObject("yCL2DGr",yCL2DGr200);

      ySens2DGr200->SetNpx(500);
      ySens2DGr200->SetNpy(500);
      yCL2DGr200->SetNpx(500);
      yCL2DGr200->SetNpy(500);
    }

    loc = std::distance(previousResultsNames.begin(),std::find(previousResultsNames.begin(),previousResultsNames.end(),"CCM200 Low Thresh."));
    if (loc != kPrevSize && loc != loc2) {
      rootFile = std::make_shared<TFile>(previousResults.at(loc).c_str(),"READ");
      rootFile->GetObject("ySens2DGr",ySens2DGr200Low);
      rootFile->GetObject("yCL2DGr",yCL2DGr200Low);

      ySens2DGr200Low->SetNpx(500);
      ySens2DGr200Low->SetNpy(500);
      yCL2DGr200Low->SetNpx(500);
      yCL2DGr200Low->SetNpy(500);
    }

    TH1D * hist = 0;
    double bestCL = 9e9;
    double bestCLMV = 9e9;
    double bestCLMChi = 9e9;
    double bestSens = 9e9;
    double bestSensMV = 9e9;
    double bestSensMChi = 9e9;
    keepReco.resize(2);
    keep.resize(2);
    rootFile = std::make_shared<TFile>(previousResults.front().c_str(),"READ");
    for (auto & result : gResults) {
      mv = std::get<0>(result);
      mchi = std::get<1>(result);
      rootFile->GetObject(Form("%s_proj_1_keepReco_%.1f_%.1f",dmHistName.c_str(),mv,mchi),hist);
      if (!hist) {
        continue;
      }
      hist->GetXaxis()->UnZoom();
      keepRecoAll.emplace_back(std::shared_ptr<TH1D>(
            dynamic_cast<TH1D*>(hist->Clone(Form("%s_proj_1_keepReco_%.1f_%.1f_localAll",dmHistName.c_str(),mv,mchi)))));
      double sensEval = ySens2DGr->Interpolate(mv/1e3,mchi/1e3);
      if (sensEval < bestSens && sensEval >= 1e-14 && mchi > 1) {
        bestSens = sensEval;
        bestSensMV = mv;
        bestSensMChi = mchi;
        double scale = sensEval/std::pow(gDefaultEpsilon,4.0)/gDefaultAlphaD;
        keepReco.front() = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(hist->Clone(Form("%s_proj_1_keepReco_%.1f_%.1f_local",dmHistName.c_str(),mv,mchi))));
        keepReco.front()->Scale(scale);
        keepReco.front()->SetLineColor(gkBlue);
        keepReco.front()->SetMarkerColor(gkBlue);
        keepReco.front()->SetLineStyle(7);
        keepReco.front()->SetLineWidth(2);
        keep.front() = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(hist->Clone(Form("%s_proj_1_keepReco_%.1f_%.1f_local2",dmHistName.c_str(),mv,mchi))));
        keep.front()->Scale(scale);
        keep.front()->SetLineColor(gkBlue);
        keep.front()->SetMarkerColor(gkBlue);
        keep.front()->SetLineStyle(7);
        keep.front()->SetLineWidth(2);
      }
      double clEval = yCL2DGr->Interpolate(mv/1e3,mchi/1e3);
      if (clEval < bestCL && clEval >= 1e-14 && mchi > 10) {
        bestCL = clEval;
        bestCLMV = mv;
        bestCLMChi = mchi;
        double scale = clEval/std::pow(gDefaultEpsilon,4.0)/gDefaultAlphaD;
        MsgDebug(2,MsgLog::Form("CL Scale: %f",scale));
        keepReco.back() = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(hist->Clone(Form("%s_proj_1_keepReco_%.1f_%.1f_local",dmHistName.c_str(),mv,mchi))));
        keepReco.back()->Scale(scale);
        keepReco.back()->SetLineColor(gkVermilion);
        keepReco.back()->SetMarkerColor(gkVermilion);
        keepReco.back()->SetLineStyle(3);
        keepReco.back()->SetLineWidth(2);
        keep.back() = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(hist->Clone(Form("%s_proj_1_keepReco_%.1f_%.1f_local2",dmHistName.c_str(),mv,mchi))));
        keep.back()->Scale(scale);
        keep.back()->SetLineColor(gkVermilion);
        keep.back()->SetMarkerColor(gkVermilion);
        keep.back()->SetLineStyle(3);
        keep.back()->SetLineWidth(2);
      }
      delete hist;
      hist = nullptr;
    }
    if (hist) {
      delete hist;
    }

    MsgInfo(MsgLog::Form("Best Sens (mv,mchi) (%f,%f) = %g",bestSensMV,bestSensMChi,bestSens));
    MsgInfo(MsgLog::Form("Best CL (mv,mchi) (%f,%f) = %g",bestCLMV,bestCLMChi,bestCL));
  }

  canvas->Clear();
  gPad->SetMargin(0.12,0.05,0.1,0.05);
  gPad->SetLogy(true);
  double maximum = 0;
  double minimum = 9e9;
  for (auto & hist : keepRecoAll) {
    maximum = std::max(maximum,hist->GetMaximum());
    minimum = std::min(minimum,hist->GetMinimum(1e-15));
  }
  count = 0;
  for (auto & hist : keepRecoAll) {
    if (count == 0) {
      hist->SetMaximum(maximum*1.05);
      hist->SetMinimum(minimum*0.95);
      hist->GetYaxis()->SetTitleOffset(1.4);
      hist->GetXaxis()->SetRange(gStartBin,gEndBin-1);
      hist->Draw("hist");
    } else {
      hist->Draw("same hist");
    }
    ++count;
  }
  canvas->cd();
  canvas->Update();
  canvas->Print(Form("%sdm_energy_true_dist_mv_3mchi_ccm%s_recoOnly.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

  canvas->Clear();
  bkgHist->SetLineColor(kGray+2);
  bkgHist->SetFillColor(TColor::GetColorTransparent(kGray+2,0.25));

  auto pad1 = new TPad("pad1","pad1",0.0,0.25,1.0,1.0);
  auto pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.25);
  canvas->cd();
  pad1->Draw();
  pad1->cd();
  gPad->SetMargin(0.12,0.02,0.02,0.05);
  //canvas->cd(cut+1);
  gPad->SetLogx(false);
  gPad->SetLogy(false);
  dataHist->SetLineColor(1);
  dataHist->SetMarkerStyle(1);
  dataHist->SetMarkerSize(0.75);
  dataHist->SetLineWidth(2);
  //dataHist->SetYTitle("Events/10^{20} POT");
  dataHist->SetYTitle("Events/PE");
  dataHist->GetXaxis()->SetTitleOffset(999999);
  dataHist->GetXaxis()->SetLabelOffset(999999);
  dataHist->GetYaxis()->SetTitleFont(103);
  dataHist->GetYaxis()->SetLabelFont(103);
  dataHist->GetYaxis()->SetTitleSize(14);
  dataHist->GetYaxis()->SetLabelSize(10);
  dataHist->GetYaxis()->SetTitleOffset(1.8);

  auto bkgFracError = dynamic_cast<TH1D*>(bkgHist->Clone("bkgFracError"));

  TLatex latex;
  latex.SetTextFont(103);
  latex.SetTextSize(14);
  latex.SetTextColor(gkBlack);
  MsgDebug(2,MsgLog::Form("Bkg Hist before dm addition %f",bkgHist->Integral()));
  MsgDebug(2,MsgLog::Form("Number of bins data (%d) bkg (%d) excess (%d) keepFront (%d) keepBack(%d) keepRecoFront (%d) keepRecoBack (%d)",
        dataHist->GetNbinsX(),
        bkgHist->GetNbinsX(),
        excessHist->GetNbinsX(),
        keep.front()->GetNbinsX(),
        keep.back()->GetNbinsX(),
        keepReco.front()->GetNbinsX(),
        keepReco.back()->GetNbinsX()));
  for (int binX = 1; binX < gNumRecoBins+1; ++binX) {
    bkgFracError->SetBinContent(binX,1.0);
    bkgFracError->SetBinError(binX,bkgHist->GetBinError(binX)/bkgHist->GetBinContent(binX));
    excessHist->SetBinContent(binX,dataHist->GetBinContent(binX)-bkgHist->GetBinContent(binX));
    excessHist->SetBinError(binX,std::sqrt(bkgHist->GetBinContent(binX)));
    dataHist->SetBinError(binX,bkgHist->GetBinError(binX));
    if (binX < gEndBin && binX >= gStartBin) {
      for (size_t i = 0; i < keepReco.size(); ++i) {
        double binWidth = keepReco.at(i)->GetXaxis()->GetBinWidth(binX);
        double content = keepReco.at(i)->GetBinContent(binX);
        double error = 0;
        for (int sys = 0; sys < 4; ++sys) {
          error += sysFracErrorMatrix.at(sys)(binX-gStartBin,binX-gStartBin)*content*content*binWidth*binWidth;
        }
        keepReco.at(i)->SetBinError(binX,std::sqrt(error));
        keep.at(i)->SetBinError(binX,std::sqrt(error));
      }
    }
    for (auto & hist : keepReco) {
      double binWidth = hist->GetXaxis()->GetBinWidth(binX);
      double before = hist->GetBinContent(binX);
      hist->SetBinContent(binX,hist->GetBinContent(binX)*binWidth);
      MsgDebug(2,MsgLog::Form("Before %f After Scale %f",before,hist->GetBinContent(binX)));
    }
    for (auto & hist : keep) {
      double binWidth = hist->GetXaxis()->GetBinWidth(binX);
      double before = hist->GetBinContent(binX);
      double newValue = before*binWidth+bkgHist->GetBinContent(binX);
      hist->SetBinContent(binX,newValue);
      MsgDebug(2,MsgLog::Form("Before %f After Scale %f Bkg %f NewSum %f (%f)",before,before*binWidth,bkgHist->GetBinContent(binX),newValue,hist->GetBinContent(binX)));
    }
  }
  MsgDebug(2,MsgLog::Form("Bkg Hist after dm addition %f",bkgHist->Integral()));

  TLegend * legTime;
  legTime = new TLegend(0.40,0.40,0.94,0.94);
  legTime->SetBorderSize(0);
  legTime->SetFillStyle(0);
  legTime->SetTextFont(102);
  //double error = 0;
  //double integralBkg = bkgHist->IntegralAndError(gStartBin,gEndBin,error,"width");
  legTime->AddEntry(bkgHist.get(),"Measured Bkg.","f");
  legTime->AddEntry(dataHist.get(),"Beam ROI (stat. error)","lep");
      //Form("Background (%.2f Events)",integralBkg),"f");

  //double integral = dataHist->IntegralAndError(gStartBin,gEndBin,error,"width");
      //Form("Data (%.0f Events)",integral),"lep");
  //integral = excessHist->IntegralAndError(gStartBin,gEndBin,error,"width");
  //legTime->AddEntry((TObject*)0,
  //    Form("Excess (%.2f #pm %.2f Events)",integral,std::sqrt(integralBkg)),"");
  legTime->AddEntry(keep.front().get(),"Bkg. + DM1 (sys. error)","lep");
  legTime->AddEntry(keep.back().get(),"Bkg. + DM2 (sys. error)","lep");

  dataHist->Scale(1.0,"width");
  MsgDebug(2,MsgLog::Form("Bkg Hist before bin scaling %f",bkgHist->Integral()));
  bkgHist->Scale(1.0,"width");
  MsgDebug(2,MsgLog::Form("Bkg Hist after bin scaling %f",bkgHist->Integral()));
  excessHist->Scale(1.0,"width");

  dataHist->GetXaxis()->SetRange(gStartBin,gEndBin-1);
  dataHist->Draw("axis");
  bkgHist->Draw("same hist");
  //bkgHist->Draw("same E");
  for (auto & hist : keep) {
    hist->Scale(1.0,"width");
    hist->Draw("same E");
  }
  dataHist->Draw("same E");
  legTime->Draw();
  gPad->SetLogy(false);
  gPad->Update();

  canvas->cd();
  pad2->Draw();
  pad2->cd();
  gPad->SetMargin(0.12,0.02,0.2,0.05);
  gPad->SetLogx(false);
  excessHist->GetXaxis()->SetRange(gStartBin,gEndBin-1);
  excessHist->SetMarkerStyle(dataHist->GetMarkerStyle());
  excessHist->SetMarkerColor(dataHist->GetMarkerColor());
  excessHist->SetMarkerSize(dataHist->GetMarkerSize());
  excessHist->SetLineStyle(dataHist->GetLineStyle());
  excessHist->SetLineColor(dataHist->GetLineColor());
  excessHist->SetLineWidth(dataHist->GetLineWidth());

  excessHist->SetTitle("");
  excessHist->SetXTitle("Energy (PE)");
  excessHist->SetYTitle("Events/PE");
  excessHist->GetXaxis()->SetRange(gStartBin,gEndBin-1);
  excessHist->GetXaxis()->SetTitleFont(103);
  excessHist->GetXaxis()->SetLabelFont(103);
  excessHist->GetXaxis()->SetTitleSize(14);
  excessHist->GetXaxis()->SetLabelSize(10);
  excessHist->GetYaxis()->SetTitleFont(103);
  excessHist->GetYaxis()->SetLabelFont(103);
  excessHist->GetYaxis()->SetTitleSize(14);
  excessHist->GetYaxis()->SetLabelSize(10);

  excessHist->GetXaxis()->SetTitleOffset(3.0);
  excessHist->GetYaxis()->SetTitleOffset(1.8);
  excessHist->Draw("E");

  for (auto & hist : keepReco) {
    hist->Scale(1.0,"width");
    hist->Draw("same E");
  }
  //excessHist->SetMaximum(1.45);
  //excessHist->SetMinimum(0.8);
  //excessHist->Draw("same E");
  double error = 0;
  double integral = excessHist->IntegralAndError(gStartBin,gEndBin,error);
  MsgInfo(MsgLog::Form("Excess is %f +- %f",integral,error));

  latex.DrawLatex(9,400,"Background Subtracted");

  TLine line;
  line.SetLineColor(kGray+2);
  line.SetLineStyle(9);
  line.DrawLine(3,0,14,0);

  gPad->Update();

  canvas->Print(Form("%sdm_data_compare_mv_3mchi_ccm%s.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

  ySens2DGr->SetNpx(500);
  ySens2DGr->SetNpy(500);
  yCL2DGr->SetNpx(500);
  yCL2DGr->SetNpy(500);

  ySens2DGr200->SetNpx(500);
  ySens2DGr200->SetNpy(500);
  yCL2DGr200->SetNpx(500);
  yCL2DGr200->SetNpy(500);

  ySens2DGr200Low->SetNpx(500);
  ySens2DGr200Low->SetNpy(500);
  yCL2DGr200Low->SetNpx(500);
  yCL2DGr200Low->SetNpy(500);

  if (infile.is_open()) {
    infile.close();
  }

  //ySens2DGr->SetNpx(mVmChiVec.size()/2.0);
  //yCL2DGr->SetNpx(mVmChiVec.size()/2.0);
  //ySens2DGr->SetNpy(mVmChiVec.size()/2.0);
  //yCL2DGr->SetNpy(mVmChiVec.size()/2.0);
  //ySens2DGr->Print();
  //yCL2DGr->Print();


  //TGraphAsymmErrors * grBest = new TGraphAsymmErrors();
  //grBest->SetLineColor(TColor::GetColorTransparent(kGray+2,0.25));
  //grBest->SetFillColor(grBest->GetLineColor());

  std::vector<double> binsY;
  for (int i=0; i < 10; ++i) {
    binsY.push_back(static_cast<double>(i)/10000.0*1e3);
  }
  for (int i = 1; i < 2000; ++i) {
    binsY.push_back(static_cast<double>(i)/1000.0*1e3);
  }

  std::vector<double> binsY2;
  for (double x=0; x < 1; x+=1e-4) {
    binsY2.push_back(x*1e3);
  }

  std::sort(binsY.begin(),binsY.end());
  TH1D * grBest = new TH1D("grBest","",binsY.size()-1,&binsY.front());
  grBest->SetFillColor(TColor::GetColorTransparent(kGray+2,0.25));
  grBest->SetLineColor(TColor::GetColorTransparent(kGray+2,0.25));
  grBest->SetMarkerColor(TColor::GetColorTransparent(kGray+2,0.25));
  grBest->SetMarkerStyle(0);
  grBest->SetMarkerSize(0);

  std::vector<TGraphAsymmErrors*> grOld;


  TGraph * grYCCM200 = new TGraph();
  grYCCM200->SetLineColor(gkVermilion);
  grYCCM200->SetLineStyle(2);
  grYCCM200->SetLineWidth(3);

  TGraph * grYCCM200Low = new TGraph();
  grYCCM200Low->SetLineColor(gkRedPurple);
  grYCCM200Low->SetLineStyle(2);
  grYCCM200Low->SetLineWidth(3);

  TGraph * grYCCM200_1 = new TGraph();
  grYCCM200_1->SetLineColor(gkBlue);
  grYCCM200_1->SetLineStyle(3);
  grYCCM200_1->SetLineWidth(3);

  TGraph * grYCCM200_10 = new TGraph();
  grYCCM200_10->SetLineColor(gkBlue);
  grYCCM200_10->SetLineStyle(4);
  grYCCM200_10->SetLineWidth(3);

  TGraph * grYCCM200_1000 = new TGraph();
  grYCCM200_1000->SetLineColor(gkBlue);
  grYCCM200_1000->SetLineStyle(5);
  grYCCM200_1000->SetLineWidth(3);

  TGraph * grYMBN = new TGraph();
  grYMBN->SetLineColor(kGray+2);
  grYMBN->SetLineStyle(1);
  grYMBN->SetLineWidth(1);

  TGraph * grYMBe = new TGraph();
  grYMBe->SetLineColor(kGray+2);
  grYMBe->SetLineStyle(1);
  grYMBe->SetLineWidth(1);

  TGraph * grYLSND = new TGraph();
  grYLSND->SetLineColor(kGray+2);
  grYLSND->SetLineStyle(1);
  grYLSND->SetLineWidth(1);

  TGraph * grYE137 = new TGraph();
  grYE137->SetLineColor(kGray+2);
  grYE137->SetLineStyle(1);
  grYE137->SetLineWidth(1);

  if (infile.is_open()) {
    infile.close();
  }
  infile.open(Form("%sy%.0f_%.1f_miniboone_n_lim.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    mv *= 1e3;
    grYMBN->SetPoint(grYMBN->GetN(),mchi,scaling);

    //if (mv < mvValues.back() && mchi < mchiValues.back()) {
    if (mv <= 130) {
      double sensEval = ySens2DGr->Interpolate(mv/1e3,mchi/1e3);
      double clEval = yCL2DGr->Interpolate(mv/1e3,mchi/1e3);
      double newEpSens = std::pow(sensEval/newAlphaD,1.0/4.0);
      double newEpCL = std::pow(clEval/newAlphaD,1.0/4.0);
      double sens = newEpSens*newEpSens*newAlphaD*std::pow(mchi/mv,4.0);
      double cl = newEpCL*newEpCL*newAlphaD*std::pow(mchi/mv,4.0);
      grY44ns->SetPoint(grY44ns->GetN(),mchi,sens);
      grY->SetPoint(grY->GetN(),mchi,cl);

      sensEval = ySens2DGr200->Interpolate(mv/1e3,mchi/1e3);
      newEpSens = std::pow(sensEval/newAlphaD,1.0/4.0);
      sens = newEpSens*newEpSens*newAlphaD*std::pow(mchi/mv,4.0);
      grYCCM200->SetPoint(grYCCM200->GetN(),mchi,sens);

      sensEval = ySens2DGr200Low->Interpolate(mv/1e3,mchi/1e3);
      newEpSens = std::pow(sensEval/newAlphaD,1.0/4.0);
      sens = newEpSens*newEpSens*newAlphaD*std::pow(mchi/mv,4.0);
      grYCCM200Low->SetPoint(grYCCM200Low->GetN(),mchi,sens);
    } else {
      grY44ns->SetPoint(grY44ns->GetN(),mchi,9e9);
      grY->SetPoint(grY->GetN(),mchi,9e9);
      grYCCM200->SetPoint(grYCCM200->GetN(),mchi,9e9);
      grYCCM200Low->SetPoint(grYCCM200Low->GetN(),mchi,9e9);
    }
    //}
  }
  infile.close();

  //std::shared_ptr<TGraphSmooth> grSmooth = std::make_shared<TGraphSmooth>("normal");
  //grY = grSmooth->Approx(grY,"min",grY->GetN()*1000);

  infile.open(Form("%sy%.0f_%.1f_ccm_1_event.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grYCCM200_1->SetPoint(grYCCM200_1->GetN(),mchi,scaling);
  }
  infile.close();

  infile.open(Form("%sy%.0f_%.1f_ccm_10_event.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grYCCM200_10->SetPoint(grYCCM200_10->GetN(),mchi,scaling);
  }
  infile.close();

  infile.open(Form("%sy%.0f_%.1f_ccm_1000_event.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grYCCM200_1000->SetPoint(grYCCM200_1000->GetN(),mchi,scaling);
  }
  infile.close();

  infile.open(Form("%sy%.0f_%.1f_miniboone_e_lim.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grYMBe->SetPoint(grYMBe->GetN(),mchi,scaling);
  }
  infile.close();

  infile.open(Form("%sy%.0f_%.1f_lsndlim.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grYLSND->SetPoint(grYLSND->GetN(),mchi,scaling);
  }
  infile.close();

  infile.open(Form("%sy%.0f_%.1f_e137lim_withBump.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grYE137->SetPoint(grYE137->GetN(),mchi,scaling);
  }
  infile.close();

  infile.open(Form("%sy%.0f_%.1f_babar.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_rare_decay.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_babar2017.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_zprime.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_monojet.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_miniboone_n_lim.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_lsndlim.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_kpipinvisk.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_invispion.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_sensei_e.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_miniboone_e_lim.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_e137lim_withBump.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_direct_det_e.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_direct_det.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_scdms_e.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();
  infile.open(Form("%sy%.0f_%.1f_NA64.dat",otherLimitsDir.c_str(),mvmChiRatio,newAlphaD));
  AddOld(infile,grOld);
  infile.close();

  TGraph * grRelic = new TGraph();
  grRelic->SetLineColor(gkBlueGreen);
  grRelic->SetLineWidth(5);
  infile.open(Form("%sY%.0f_Relic_Density.dat",otherLimitsDir.c_str(),mvmChiRatio));
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    grRelic->SetPoint(grRelic->GetN(),mchi,scaling);
  }
  infile.close();


  for (int bin=1; bin < grBest->GetNbinsX()+1; ++bin) {
    double content = grBest->GetBinContent(bin);
    double point = (1e-6 + content)/2.0;
    double error= 1e-6 - point;
    grBest->SetBinContent(bin,point);
    grBest->SetBinError(bin,error);
  }

  //avgKeep->Scale(1.0/static_cast<double>(count));
  //avgKeep->SetLineColor(kRed);
  //avgKeep->SetLineWidth(2);
  //leg->AddEntry(avgKeep,"Average","l");
  //hs->Add(avgKeep);

  canvas->Clear();
  //TH2D * rangeHist = new TH2D("rangeHist",";m_{#chi} (GeV/c^{2});Y = #varepsilon^{2}#alpha_{D}(m_{#chi}/m_{V})^{4}",100000,1e-4,1.022,100,1e-6,2e-12);
  //rangeHist->SetBinContent(150,1);
  //rangeHist->SetMarkerSize(0);
  canvas->SetLogy(true);
  canvas->SetLogx(true);
  canvas->SetMargin(0.15,0.05,0.1,0.05);
  gStyle->SetOptStat(0);
  //rangeHist->Draw("");
  auto grBestCopy = new TH1D("grBestCopy","",binsY2.size()-1,&binsY2.front());
  grBestCopy->GetXaxis()->SetRangeUser(1.001,100);
  grBestCopy->GetYaxis()->SetTitleOffset(1.8);
  grBestCopy->GetXaxis()->SetTitleOffset(1.3);
  grBestCopy->GetXaxis()->SetTitle("m_{#chi} (MeV/c^{2})"); 
  grBestCopy->GetYaxis()->SetTitle(Form("Y = #varepsilon^{2}#alpha_{D}(m_{#chi}/m_{V})^{4} (m_{V} = %.0fm_{#chi}, #alpha_{D} = %.1f)",mvmChiRatio,newAlphaD));
  grBestCopy->SetMinimum(5e-13);
  grBestCopy->SetMaximum(2e-7);
  grBestCopy->Draw("axis");
  //grBest->Draw("same E3");
  for (auto & gr : grOld) {
    gr->DrawClone("same 03");
  }
  //for (auto & gr : grOld) {
  //  gr->SetLineColor(kGray+2);
  //  gr->Draw("same lX0");
  //}
  grRelic->Draw("same l");
  grYLSND->Draw("same l");
  grYE137->Draw("same l");
  grYMBe->Draw("same l");
  grYMBN->Draw("same l");
  grY44ns->Draw("same l");
  grY->Draw("same l");
  grYCCM200->Draw("same l");
  grYCCM200Low->Draw("same l");
  //grYCCM200_1->Draw("same l");
  //grYCCM200_10->Draw("same l");
  //grYCCM200_1000->Draw("same l");
  grBestCopy->Draw("same axis");

  TLegend * legSens = new TLegend(0.65,0.15,0.94,0.3);
  legSens->SetBorderSize(0);
  legSens->SetFillStyle(0);
  legSens->SetTextFont(103);
  legSens->SetTextSize(12);
  //legSens->SetHeader("90% Limts/Sensitivities");
  //legSens->AddEntry(grBest,"Current Best Limits","f");
  //legSens->AddEntry(grRelic,"Relic Density","l");
  //legSens->AddEntry(grYLSND,"LSND Electron Limit","l");
  //legSens->AddEntry(grYMBe,"MiniBooNE Electron Limit","l");
  //legSens->AddEntry(grYMBN,"MiniBooNE Nucleon Limit","l");
  //legSens->AddEntry(grY,"CCM120 90% Limit","l");
  //legSens->AddEntry(grY44ns,"CCM120 Sens.","l");
  legSens->SetHeader("CCM200 Signal Reach");
  legSens->AddEntry(grYCCM200_1,"> 1 Event","l");
  legSens->AddEntry(grYCCM200_10,"> 10 Event","l");
  legSens->AddEntry(grYCCM200_1000,"> 1000 Event","l");
  //legSens->Draw();
  canvas->Update();
  //TLatex latex;
  latex.SetTextFont(103);
  latex.SetTextSize(12);

  latex.SetTextColor(gkBlack);
  latex.DrawLatex(4,1.2e-8,"CCM120");
  latex.SetTextColor(gkBlack);
  latex.DrawLatex(4,5e-9,"CCM120 Sens.");
  latex.SetTextColor(gkVermilion);
  latex.DrawLatex(1.8,6e-11,"CCM200 50keV Sens.");
  latex.SetTextColor(gkRedPurple);
  latex.DrawLatex(2,9e-12,"CCM200 10keV Sens.");
  //latex.SetTextColor(gkBlack);
  //latex.DrawLatex(1.5,7e-11,"#splitline{CCM200}{> 1000 Events}");
  //latex.SetTextColor(gkBlack);
  //latex.DrawLatex(1.5,4e-12,"CCM200 > 10 Events");
  //latex.SetTextColor(gkBlack);
  //latex.DrawLatex(10,3.25e-12,"CCM200 > 1 Events");
  //latex.SetTextColor(gkBlack);
  //latex.DrawLatex(8,2.5e-9,"CCM120 Sens.");
  latex.SetTextColor(kGray+2);
  latex.DrawLatex(3,6e-10,"MB FullN");
  latex.SetTextColor(kGray+2);
  latex.DrawLatex(2,2e-12,"LSND");
  latex.SetTextColor(kGray+2);
  latex.DrawLatex(42,7e-10,"E137");
  latex.SetTextColor(kGray+2);
  latex.DrawLatex(60,3.5e-9,"MBe");
  latex.SetTextColor(gkBlueGreen);
  latex.DrawLatex(50,1e-10,"#splitline{Relic}{Density}");
  gPad->SetGridx(true);
  gPad->SetGridy(true);
  canvas->Print(Form("%sdm_energy_pred_y_energyEff_All_ccm200%s_zoom.pdf",outputDir.c_str(),postFix.c_str()),"pdf");
  canvas->SetGrayscale(true);
  canvas->Print(Form("%sdm_energy_pred_y_energyEff_All_ccm200%s_zoom_grayscale.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

  /*
     canvas->Clear();
     gPad->SetMargin(0.12,0.18,0.1,0.05);
     gPad->SetLogz(true);
     ySens2D->GetZaxis()->SetTitleOffset(1.6);
     ySens2D->SetMinimum(minSens2D*0.5);
     ySens2D->SetMaximum(maxSens2D*1.1);
     ySens2D->Draw("colz");
     canvas->Print(Form("%sdm_ep4ad_sens_2d.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

     canvas->Clear();
     gPad->SetMargin(0.12,0.18,0.1,0.05);
     gPad->SetLogz(true);
     yCL2D->GetZaxis()->SetTitleOffset(1.6);
     yCL2D->SetMinimum(minCL2D*0.5);
     yCL2D->SetMaximum(maxCL2D*1.1);
     yCL2D->Draw("colz");
     canvas->Print(Form("%sdm_ep4ad_cl_2d.pdf",outputDir.c_str(),postFix.c_str()),"pdf");
     */

  return EXIT_SUCCESS;

  /*
     std::ofstream outfile("dm_avg_energy_dist_50keV.txt");
     outfile << "Elow\tEhigh\tContent" << std::endl;
     for (int bin = 1; bin < avgKeep->GetNbinsX(); ++bin) {
     outfile << avgKeep->GetXaxis()->GetBinLowEdge(bin) << '\t' 
     << avgKeep->GetXaxis()->GetBinLowEdge(bin+1) << '\t'
     << avgKeep->GetBinContent(bin) << std::endl;
     }
     outfile.close();
     */

}

//-------------------------------------------------------------------------------------------------
std::shared_ptr<TH1D> GetTotalHist(float mv, float mchi, int axis, std::string histName, 
    std::string masterDir, std::string cvDir, std::vector<std::string> sysDir,
    std::vector<std::shared_ptr<TH1D>> * sys)
{
  THnSparseD * totalSparseHist = 0;

  std::string name = Form("%s/%s/Histograms/darkMatter_dp%.1fdm%.1f",masterDir.c_str(),cvDir.c_str(),mv,mchi);

  std::shared_ptr<TFile> file = std::make_shared<TFile>(std::string(name+"_coherent.root").c_str(),"READ");
  if (!file || file == nullptr) {
    return nullptr;
  }
  if (!file->IsOpen()) {
    return nullptr;
  }
  file->GetObject(histName.c_str(),totalSparseHist);
  if (!totalSparseHist) {
    std::cerr << "total space hist could not be found from " << histName << std::endl;
    return nullptr;
  }

  std::shared_ptr<TH1D>  hist = std::shared_ptr<TH1D>(totalSparseHist->Projection(axis));
  hist->SetDirectory(0);
  delete totalSparseHist;

  if (sys) {
    if (!sys->empty()) {
      sys->clear();
    }
    for (auto & dir : sysDir) {
      name = Form("%s/%s/Histograms/darkMatter_dp%.1fdm%.1f",masterDir.c_str(),dir.c_str(),mv,mchi);

      file = std::make_shared<TFile>(std::string(name+"_coherent.root").c_str(),"READ");
      if (!file || file == nullptr) {
        return nullptr;
      }
      if (!file->IsOpen()) {
        return nullptr;
      }
      file->GetObject(histName.c_str(),totalSparseHist);
      if (!totalSparseHist) {
        std::cerr << "total space hist could not be found from " << histName << std::endl;
        return nullptr;
      }

      sys->emplace_back(std::shared_ptr<TH1D>(totalSparseHist->Projection(axis)));
      sys->back()->SetDirectory(0);
    } // end for sysDir
  } // end if sys

  return hist;

}

//-------------------------------------------------------------------------------------------------
void CalculateSignalVariance(std::shared_ptr<TH1D> hist, const double omFrac, const double potFrac, const double pi0Frac,
    const double qfFrac) 
{
  for (int binX = 1; binX <= gNumRecoBins; ++binX) {
    double events = hist->GetBinContent(binX);
    double error = events + std::pow(events*omFrac,2.0) + 
      std::pow(events*potFrac,2.0) + std::pow(pi0Frac*events,2.0) + std::pow(events*qfFrac,2.0);
    hist->SetBinError(binX,std::sqrt(error));
  }
}

//-------------------------------------------------------------------------------------------------
void FindBest(std::ifstream & infile, TH1D *  best)
{
  double mchi = 0;
  double mv = 0;
  double scaling = 0;
  while (infile >> mv >> mchi >> scaling) {
    //if (mchi < 1e-3) {
    //  continue;
    //}
    int bin = best->FindBin(mchi);
    double content = best->GetBinContent(bin);
    if (content == 0) {
      best->SetBinContent(bin,scaling);
      continue;
    }

    best->SetBinContent(bin,std::min(scaling,content));
  }
}

//-------------------------------------------------------------------------------------------------
void AddOld(std::ifstream & infile, std::vector<TGraphAsymmErrors*> & grOld)
{
  grOld.push_back(new TGraphAsymmErrors());
  grOld.back()->SetFillColor(gkLightGray);
  grOld.back()->SetLineColor(gkLightGray);
  double mchi = 0;
  double mv = 0;
  double scaling = 0;
  while (infile >> mv >> mchi >> scaling) {
    mchi *= 1e3;
    //if (mchi < 1e-3) {
    //  continue;
    //}

    grOld.back()->SetPoint(grOld.back()->GetN(),mchi,scaling);
    grOld.back()->SetPointError(grOld.back()->GetN()-1,0,0,0,9999-scaling);
  }
}

//-------------------------------------------------------------------------------------------------
void ScaleDMSys(const TMatrixD & cv, const std::vector<TMatrixD> & fracError, TMatrixD & sysError, bool print) 
{
  const int kNRows = fracError.front().GetNrows();
  for (int x = 0; x < kNRows; ++x) {
    for (int y = 0; y < kNRows; ++y) {
      sysError(x,y) = 0;
      for (auto & matrix : fracError) {
        sysError(x,y) += matrix(x,y)*cv(x,0)*cv(y,0);
      }
    }
  }
  if (print) {
    double cvTotal = cv.Sum();
    MsgInfo(MsgLog::Form("sysSum %f, frac %f",std::sqrt(sysError.Sum()),std::sqrt(sysError.Sum())/cvTotal));
  }
}

//-------------------------------------------------------------------------------------------------
double CalculateChi2(const TMatrixD & data, const TMatrixD & pred, const TMatrixD & sys)
{
  TMatrixD sysInv(sys);
  double det = sysInv.Determinant();
  if (det < 0.5) {
    MsgError("determinate is less than 0.5. Printing covariance matrix and then exiting");
    sysInv.Print();
    abort();
  }
  sysInv.InvertFast();
  TMatrixD diff(data);
  diff -= pred;

  return TMatrixD(TMatrixD(diff,TMatrixD::kTransposeMult,sysInv),TMatrixD::kMult,diff)(0,0);
}

void GetToyMC(std::vector<TMatrixD> & fakeDataVec, const TMatrixD & cv, const TMatrixD & sys, std::mt19937 & mt)
{
  std::normal_distribution<double> normal(0,1);
  TDecompChol chol(sys);
  chol.Decompose();
  auto U = chol.GetU();
  TMatrixD ran(cv);
  for (auto & fakeData : fakeDataVec) {
    fakeData.ResizeTo(cv);
    for (int x = 0; x < ran.GetNrows(); ++x) {
      ran(x,0) = normal(mt);
    }
    for (int j = 0; j < cv.GetNrows(); ++j) {
      double row = 0;
      for (int k = 0; k < cv.GetNrows(); ++k) {
        row += U(k,j)*ran(k,0);
      }
      double element = cv(j,0) + row;
      if (element > 0) {
        fakeData(j,0) = element;
      } else {
        fakeData(j,0) = 0;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CalculateSensCL(const TMatrixD dm, const TMatrixD & data, const TMatrixD & bkg, 
    const std::vector<TMatrixD> & sysFracError, const TMatrixD & bkgError, const int & index)
{

  std::vector<double> nullThrows;
  double sens = 0;
  double cl = 0;

  std::vector<double> scaleVector;
  for (double exp = -12; exp < 7; ++exp) {
    for (double digit = 1; digit < 10; ++digit) {
      scaleVector.emplace_back(digit*std::pow(10.0,exp));
    }
  }

  std::map<double,std::pair<std::shared_ptr<TMatrixD>,std::shared_ptr<TMatrixD>>> mapBkgPlusSig;

  std::shared_ptr<TMatrixD> throwHist = std::make_shared<TMatrixD>(dm.GetNrows(),dm.GetNcols());
  std::shared_ptr<TMatrixD> throwHist2 = std::make_shared<TMatrixD>(dm.GetNrows(),dm.GetNcols());
  std::shared_ptr<TMatrixD> sysError = std::make_shared<TMatrixD>(sysFracError.front().GetNrows(),sysFracError.front().GetNcols());
  std::shared_ptr<TMatrixD> sysError2 = std::make_shared<TMatrixD>(sysFracError.front().GetNrows(),sysFracError.front().GetNcols());
  //std::shared_ptr<MatrixD> zeroHist(dm.GetNrows(),dm.GetNcols());

  //throwHist->operator=(dm);
  //ScaleDMSys(*throwHist,sysFracError,*sysError,true);

  const double ep4aD = std::pow(gDefaultEpsilon,4.0)*gDefaultAlphaD;
  std::map<double,double> dataDeltaChi2;
  double minDataChi2 = 9e9;
  //double minDataScale = 0;
  for (const auto & testScale : scaleVector) {
    throwHist->operator=(dm);

    throwHist->operator*=(testScale);
    ScaleDMSys(*throwHist,sysFracError,*sysError);
    throwHist->operator+=(bkg);
    sysError->operator+=(bkgError);
    mapBkgPlusSig.emplace(testScale,std::make_pair((std::make_shared<TMatrixD>(*throwHist)),(std::make_shared<TMatrixD>(*sysError))));
    double chi2 = CalculateChi2(data,*throwHist,*sysError);

    dataDeltaChi2.emplace(testScale*ep4aD,chi2);
    if (chi2 < minDataChi2) {
      minDataChi2 = chi2;
      //minDataScale = testScale;
    }
  } // end for scaleVector looping to calculate deltaChi2 with data

  for (auto & p : dataDeltaChi2) {
    p.second -= minDataChi2;
  }

  // null throws
  std::vector<TMatrixD> fakeDataVec(gNumNullFakeData);
  std::mt19937 mt(123456789);
  while (nullThrows.size() < gNumNullFakeData) {
    GetToyMC(fakeDataVec,bkg,bkgError,mt);
    for (auto & fakeData : fakeDataVec) {
      double nullChi2 = CalculateChi2(fakeData,bkg,bkgError);

      double minChi2 = 9e9;
      //double minScale = 0;
      for (const auto & testScale : scaleVector) {
        auto itMap = mapBkgPlusSig.find(testScale);
        if (itMap == mapBkgPlusSig.end()) {
          throwHist->operator=(dm);
          throwHist->operator*=(testScale);
          ScaleDMSys(*throwHist,sysFracError,*sysError);
          throwHist->operator+=(bkg);
          sysError->operator+=(bkgError);
          mapBkgPlusSig.emplace(testScale,std::make_pair((std::make_shared<TMatrixD>(*throwHist)),(std::make_shared<TMatrixD>(*sysError))));
        } else {
          throwHist = itMap->second.first;
          sysError = itMap->second.second;
        }

        double chi2 = CalculateChi2(fakeData,*throwHist,*sysError);

        if (chi2 < minChi2) {
          minChi2 = chi2;
          //minScale = testScale;
        }
      } // end for scaleVector
      if (nullChi2 - minChi2 < 0) {
        //  --index;
        continue;
      }
      nullThrows.emplace_back(nullChi2 - minChi2);
    } // end for fakeData
  } // end while nullThrows

  std::sort(nullThrows.begin(),nullThrows.end());
  size_t nullIndex90 = nullThrows.size()*0.90;
  double null90Chi2 = nullThrows.at(nullIndex90);
  
  std::shared_ptr<TGraph> sigNullGraph = std::make_shared<TGraph>();
  std::map<double,std::vector<double>> sigNullThrows;
  std::map<double,std::vector<double>> sigSigThrows;

  for (const auto & testScale : scaleVector) {
    auto itMap = mapBkgPlusSig.find(testScale);
    if (itMap == mapBkgPlusSig.end()) {
      throwHist2->operator=(dm);
      throwHist2->operator*=(testScale);
      ScaleDMSys(*throwHist2,sysFracError,*sysError2);
      throwHist2->operator+=(bkg);
      sysError2->operator+=(bkgError);
      //for (int row = 0; row < sysError->GetNrows(); ++row) {
      //  throwHist2->operator()(row,0) += bkg(row,0);
      //  for (int col = 0; col < sysError->GetNcols(); ++col) {
      //    sysError2->operator()(row,col) += bkgError(row,col);
      //  }
      //}
      mapBkgPlusSig.emplace(testScale,std::make_pair((std::make_shared<TMatrixD>(*throwHist2)),(std::make_shared<TMatrixD>(*sysError2))));
    } else {
      throwHist2 = itMap->second.first;
      sysError2 = itMap->second.second;
    }

    auto itSigSim = sigSigThrows.find(testScale*ep4aD);
    if (itSigSim == sigSigThrows.end()) {
      sigSigThrows.emplace(testScale*ep4aD,std::vector<double>());
      itSigSim = sigSigThrows.find(testScale*ep4aD);
    }

    auto itNullSim = sigNullThrows.find(testScale*ep4aD);
    if (itNullSim == sigNullThrows.end()) {
      sigNullThrows.emplace(testScale*ep4aD,std::vector<double>());
      itNullSim = sigNullThrows.find(testScale*ep4aD);
    }

    fakeDataVec.clear();
    fakeDataVec.resize(gNumSigFakeData);
    std::random_device rd;
    mt.seed(rd());
    while (itSigSim->second.size() < gNumSigFakeData) {
      GetToyMC(fakeDataVec,*throwHist2,*sysError2,mt);
      for (auto & fakeData : fakeDataVec) {
        double nullChi2 = CalculateChi2(fakeData,bkg,bkgError);
        double sigChi2 = CalculateChi2(fakeData,*throwHist2,*sysError2);

        double minChi2 = 9e9;
        //double minScale = 0;
        for (const auto & testScale2 : scaleVector) {
          itMap = mapBkgPlusSig.find(testScale2);
          if (itMap == mapBkgPlusSig.end()) {
            throwHist->operator=(dm);
            throwHist->operator*=(testScale2);
            ScaleDMSys(*throwHist,sysFracError,*sysError);
            throwHist->operator+=(bkg);
            sysError->operator+=(bkgError);
            mapBkgPlusSig.emplace(testScale2,std::make_pair((std::make_shared<TMatrixD>(*throwHist)),(std::make_shared<TMatrixD>(*sysError))));
          } else {
            throwHist = itMap->second.first;
            sysError = itMap->second.second;
          }
          double chi2 = CalculateChi2(fakeData,*throwHist,*sysError);
          if (chi2 < minChi2) {
            minChi2 = chi2;
            //minScale = testScale2;
          }
          if (chi2 > minChi2 * 9) {
            break;
          }// end if chi2 < minChi2
        } // end for testScale2
        if (nullChi2 - minChi2 < 0) {
          //--index;
          continue;
        }
        itNullSim->second.emplace_back(nullChi2 - minChi2);
        itSigSim->second.emplace_back(sigChi2 - minChi2);
      } // end for fakeData
    } // end while itSigSim.size() < 1000
    std::sort(itNullSim->second.begin(),itNullSim->second.end());
    double mean = std::accumulate(itNullSim->second.begin(),itNullSim->second.end(),0.0)/static_cast<double>(itNullSim->second.size());
    double low = itNullSim->second.front();
    sigNullGraph->SetPoint(sigNullGraph->GetN(),mean,testScale);
    if (low > null90Chi2*2.0) {
      sens = sigNullGraph->Eval(null90Chi2)*ep4aD;

      std::shared_ptr<TGraph> dataGraph = std::make_shared<TGraph>();
      for (auto & p : sigSigThrows) {
        std::sort(p.second.begin(),p.second.end());
        int indexP = p.second.size()*0.9;
        double clchi290 = p.second.at(indexP);
        auto itData = dataDeltaChi2.find(p.first);
        if (itData == dataDeltaChi2.end()) {
          continue;
        }
        dataGraph->SetPoint(dataGraph->GetN(),clchi290 - itData->second,p.first);
      }

      cl = dataGraph->Eval(0);

      break;
    }
  } // end for testScale

  gMTX.lock();
  MsgInfo(MsgLog::Form("mV %.1f mChi %.1f null90Chi2 %f",std::get<0>(gResults.at(index)),std::get<1>(gResults.at(index)),null90Chi2));
  std::get<2>(gResults.at(index)) = sens;
  std::get<3>(gResults.at(index)) = cl;
  MsgInfo(MsgLog::Form("mV %.1f mChi %.1f Sens %g CL %g",
        std::get<0>(gResults.at(index)),std::get<1>(gResults.at(index)),std::get<2>(gResults.at(index)),std::get<3>(gResults.at(index))));
  gMTX.unlock();

  return;
}


//continue;

/*
  double value = 0;
   TH1D * nullHist = new TH1D("nullHist","Null Fake Data;#Delta#chi^{2}_{null};Probability Distribution Function",
   nullThrows.size(),nullThrows.front(),nullThrows.back());
   for (auto & t : nullThrows) {
   nullHist->Fill(t);
   }
   nullHist->Scale(1.0/static_cast<double>(nullThrows.size()));

   std::vector<double> yAxis;
   for (double exp = -12; exp < 7; ++exp) {
   for (double digit = 1; digit < 10; ++digit) {
   yAxis.emplace_back(digit*std::pow(10.0,exp)*ep4aD);
   }
   }
   yAxis.emplace_back(1*std::pow(10.0,7)*ep4aD);

   std::vector<double> xAxis;
   for (double exp = -12; exp < 8; ++exp) {
   for (double digit = 1; digit < 10; digit+=0.25) {
   xAxis.emplace_back(digit*std::pow(10.0,exp));
   }
   }

   TH2D * nullSigHist = new TH2D("nullSigHist","Signal Fake Data;#Delta#chi^{2}_{null};#varepsilon^{4}#alpha_{D};Row Normalized",
   xAxis.size()-1,&xAxis.front(),
   yAxis.size()-1,&yAxis.front());
   TGraph * nullSigHistMean = new TGraph();

   for (auto & p : sigNullThrows) {
   double y = p.first;
   double total = p.second.size();
   double mean = std::accumulate(p.second.begin(),p.second.end(),0.0)/total;
   nullSigHistMean->SetPoint(nullSigHistMean->GetN(),mean,y);
   for (auto & value : p.second) {
   nullSigHist->Fill(value,y,1.0/total);
   }
   }
   nullSigHist->Print();

   TH1D * sigSigHist = new TH1D("sigSigHist",";#varepsilon^{4}#alpha_{D};#Delta#chi^{2}_{sig}",yAxis.size()-1,&yAxis.front());
   TH1D * dataSigHist = new TH1D("dataSigHist",";#varepsilon^{4}#alpha_{D};#Delta#chi^{2}_{sig}",yAxis.size()-1,&yAxis.front());


   for (auto & p : dataDeltaChi2) {
   std::cout << p.first << '\t' << p.second << std::endl;
   dataSigHist->Fill(p.first,p.second);
   }
//dataSigHist->Print("all");

double lastY = 0;
for (auto & p : sigSigThrows) {
double y = p.first;
lastY = y;
int indexP = p.second.size()*0.9;
double clchi290 = p.second.at(indexP);
sigSigHist->Fill(y,clchi290);
}
//sigSigHist->Print("all");

TCanvas * canvas = new TCanvas("canvas","canvas",500,500);
gStyle->SetOptStat(0);
canvas->SetMargin(0.12,0.05,0.1,0.1);

nullHist->Draw("hist");
TLine * line = new TLine(null90Chi2,0,null90Chi2,nullHist->GetMaximum());
line->SetLineColor(kRed);
line->SetLineWidth(2);
nullHist->SetLineWidth(2);
nullHist->SetLineColor(kBlack);
nullHist->GetYaxis()->SetTitleOffset(1.6);
nullHist->GetXaxis()->SetTitleOffset(1.4);
canvas->SetLogy(true);
canvas->SetLogx(true);
line->Draw();
canvas->Update();
canvas->Print(Form("%ssample_null_throws.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

canvas->Clear();
canvas->SetMargin(0.12,0.15,0.1,0.1);
nullSigHist->Draw("colz");
nullSigHist->GetYaxis()->SetRangeUser(yAxis.front(),lastY);
nullSigHist->GetYaxis()->SetTitleOffset(1.6);
nullSigHist->GetZaxis()->SetTitleOffset(1.6);
nullSigHist->GetXaxis()->SetTitleOffset(1.4);
nullSigHistMean->SetMarkerColor(kBlack);
nullSigHistMean->SetMarkerSize(1);
nullSigHistMean->SetMarkerStyle(20);
nullSigHistMean->Draw("same p");
canvas->SetLogy(true);
canvas->SetLogz(true);
canvas->SetLogx(true);
line->DrawLine(null90Chi2,yAxis.front(),null90Chi2,yAxis.back());
canvas->Update();
canvas->Print(Form("%ssample_sig_null_throws.pdf",outputDir.c_str(),postFix.c_str()),"pdf");

THStack * hs = new THStack("hs",";#varepsilon^{4}#alpha_{D};#Delta#chi^{2}_{sig}");
hs->Add(sigSigHist,"hist");
hs->Add(dataSigHist, "hist p");

sigSigHist->SetLineColor(kRed);
dataSigHist->SetMarkerColor(kBlack);
dataSigHist->SetMarkerStyle(20);

canvas->Clear();
canvas->SetMargin(0.12,0.05,0.1,0.05);
canvas->SetLogx(true);
hs->Draw("nostack");
hs->GetHistogram()->GetXaxis()->SetRangeUser(yAxis.front(),lastY);
hs->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
hs->GetHistogram()->GetXaxis()->SetTitleOffset(1.4);

TLegend * leg = new TLegend(0.3,0.2,0.5,0.4);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(sigSigHist,"Signal Fake Data","l");
leg->AddEntry(dataSigHist,"Data","p");

leg->Draw();

canvas->Print(Form("%ssample_data_cl.pdf",outputDir.c_str(),postFix.c_str()),"pdf");


return EXIT_SUCCESS;
*/

