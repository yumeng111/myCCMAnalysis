/*!**********************************************
 * \file CCMTimeShift.cxx
 * \author R.T. Thornton (LANL)
 * \date December 3, 2020
 *
 * Main code to model the expected timing distribution
 * in the simulation.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMTimeShift.h"
#include "CCMModuleTable.h"

#include "Pulses.h"
#include "MsgLog.h"
#include "Utility.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"

//See CCMModuleTable for info
MODULE_DECL(CCMTimeShift);

//_______________________________________________________________________________________
CCMTimeShift::CCMTimeShift(const char* version) 
  : CCMModule("CCMTimeShift"),
    fPulses(nullptr),
    fRD(),
    fMT(),
    fTimeDist(kCCMBeamPi0),
    fStartNS(Utility::fgkWindowStartTime),
    fEndNS(Utility::fgkWindowEndTime),
    fBeamWidthNS(270.0),
    fTotalEvents(0),
    fOutOfWindowBeamPi0(0),
    fOutOfWindowBeamPipm(0),
    fOutOfWindowBeamMupm(0)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMTimeShift::CCMTimeShift(const CCMTimeShift& clufdr) 
: CCMModule(clufdr),
  fPulses(clufdr.fPulses),
  fRD(),
  fMT(),
  fTimeDist(clufdr.fTimeDist),
  fStartNS(clufdr.fStartNS),
  fEndNS(clufdr.fEndNS),
  fBeamWidthNS(clufdr.fBeamWidthNS),
  fTotalEvents(clufdr.fTotalEvents),
  fOutOfWindowBeamPi0(clufdr.fOutOfWindowBeamPi0),
  fOutOfWindowBeamPipm(clufdr.fOutOfWindowBeamPipm),
  fOutOfWindowBeamMupm(clufdr.fOutOfWindowBeamMupm)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMTimeShift::~CCMTimeShift()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMTimeShift::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting Time Shift for MC \"Event\"");
  }

  if (fPulses->GetNumPulses() == 0) {
    return kCCMSuccess;
  }

  double currentTimeShift = fPulses->GetTriggerTime();

  double newStartTime = 0.0;
  switch (fTimeDist) {
    case kCCMUniform: newStartTime = GetTimeFromUniform(fStartNS,fEndNS); break;
    case kCCMBeamPi0: newStartTime = GetTimeFromBeamPi0(fStartNS,fEndNS,fBeamWidthNS); break;
    case kCCMBeamPipm: newStartTime = GetTimeFromBeamPipm(fStartNS,fEndNS,fBeamWidthNS); break;
    case kCCMBeamMupm: newStartTime = GetTimeFromBeamMupm(fStartNS,fEndNS,fBeamWidthNS); break;
    default: newStartTime = 0.0;
  }

  double newStartTimeBin = (newStartTime - Utility::fgkWindowStartTime) / Utility::fgkBinWidth;
  MsgInfo(MsgLog::Form("newST = %.2f newSTB = %.2f",newStartTime,newStartTimeBin));

  fPulses->ShiftTimeOffset(newStartTimeBin - currentTimeShift);

  ++fTotalEvents;

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMTimeShift::Configure(const CCMConfig& c ) 
{

  //Initialize any parameters here
  //by reading them from the CCMConfig object.
  int seed = 0;
  c("RANSEED").Get(seed);
  if (seed == 0) {
    fMT.seed(fRD());
  } else {
    fMT.seed(seed);
  }

  std::string tempString = "";
  c("TimingDist").Get(tempString);
  for (auto & c : tempString) {
    c = toupper(c);
  }
  if (tempString.find("DEFAULT") != std::string::npos ||
      tempString.find("BEAMPI0") != std::string::npos) {
    fTimeDist = kCCMBeamPi0;
    tempString = "BEAMPI0";
  } else if (tempString.find("BEAMPIPM") != std::string::npos) {
    fTimeDist = kCCMBeamPipm;
  } else if (tempString.find("BEAMMUPM") != std::string::npos) {
    fTimeDist = kCCMBeamMupm;
  } else if (tempString.find("UNIFORM") != std::string::npos) {
    fTimeDist = kCCMUniform;
  } else {
    MsgFatal(MsgLog::Form("Must apply a valid TimingDist option, passed %s",tempString.c_str()));
  }

  c("StartNS").Get(fStartNS);
  if (fTimeDist == kCCMUniform) {
    c("EndNS").Get(fEndNS);
  } else {
    fEndNS = Utility::fgkWindowEndTime;
    c("BeamWidthNS").Get(fBeamWidthNS);
  }

  MsgInfo(MsgLog::Form("\t-Seed values %d fMT output %u",seed,fMT()));
  MsgInfo(MsgLog::Form("\t-Time Dist: %s",tempString.c_str()));
  MsgInfo(MsgLog::Form("\t-StartNS: %0.2f",fStartNS));
  MsgInfo(MsgLog::Form("\t-EndNS: %0.2f",fEndNS));
  MsgInfo(MsgLog::Form("\t-BeamWidthNS: %0.2f",fBeamWidthNS));

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMTimeShift::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Number of events: %.0f",fTotalEvents));
  MsgInfo(MsgLog::Form("Had to rerun %.0f beamPi0 numbers",fOutOfWindowBeamPi0));
  MsgInfo(MsgLog::Form("Had to rerun %.0f beamPipm numbers",fOutOfWindowBeamPipm));
  MsgInfo(MsgLog::Form("Had to rerun %.0f beamMupm numbers",fOutOfWindowBeamMupm));

  return kCCMSuccess;
}

//_______________________________________________________________________________________
double CCMTimeShift::GetTimeFromUniform(double start, double end)
{
  std::uniform_real_distribution<double> uniform(start,end);
  return uniform(fMT);
}

//_______________________________________________________________________________________
double CCMTimeShift::GetTimeFromBeamPi0(double start, double end, double beamWidth)
{
  std::array<double, 3> i{start, start+beamWidth/2.0, start+beamWidth};
  std::array<double, 3> w{0, 1, 0};
  auto beamDist = std::piecewise_linear_distribution<double>{i.begin(), i.end(), w.begin()};
  double time = beamDist(fMT);
  while (time > end) {
    time = beamDist(fMT);
    ++fOutOfWindowBeamPi0;
  }

  return time;
}

//_______________________________________________________________________________________
double CCMTimeShift::GetTimeFromBeamPipm(double start, double end, double beamWidth)
{
  std::exponential_distribution<double> pipmDist(1.0/26.033); // in ns
  double time = GetTimeFromBeamPi0(start,end,beamWidth)+pipmDist(fMT);
  while (time > end) {
    time = GetTimeFromBeamPi0(start,end,beamWidth)+pipmDist(fMT);
    ++fOutOfWindowBeamPipm;
  }
  return time;
}

//_______________________________________________________________________________________
double CCMTimeShift::GetTimeFromBeamMupm(double start, double end, double beamWidth)
{
  std::exponential_distribution<double> mupmDist(1.0/2.1969811e3); // in ns
  double time = GetTimeFromBeamPipm(start,end,beamWidth)+mupmDist(fMT);
  while (time > end) {
    time = GetTimeFromBeamPipm(start,end,beamWidth)+mupmDist(fMT);
    ++fOutOfWindowBeamMupm;
  }
  return time;
}

