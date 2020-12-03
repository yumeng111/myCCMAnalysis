/*!**********************************************
 * \file CCMPMTResponse.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to modle the response of the PMTs
 * in the simulation.
 ***********************************************/

#include "CCMConfig.h"
#include "CCMConfigParam.h"
#include "CCMPMTResponse.h"
#include "CCMModuleTable.h"

#include "PMTInfoMap.h"
#include "PMTInformation.h"
#include "MCTruth.h"
#include "SinglePulse.h"
#include "Pulses.h"
#include "MsgLog.h"
#include "Utility.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"

#include <array>
#include <cctype>

//See CCMModuleTable for info
MODULE_DECL(CCMPMTResponse);

//_______________________________________________________________________________________
CCMPMTResponse::CCMPMTResponse(const char* version) 
  : CCMModule("CCMPMTResponse"),
    fMCTruth(nullptr),
    fPulses(nullptr),
    fRD(),
    fMT(),
    fUniform(0,1),
    fSPEWeights(),
    fWaveforms(),
    fWaveformsHist(),
    fFile(nullptr),
    fHighEnergy(4.0),
    fLowEnergy(2.0),
    fQE(0.2),
    fSquareWF(false),
    fTriangleWF(true),
    fTotalHits(0),
    fAfterQE(0),
    fAfterADC(0),
    fTotalPulses(0),
    fAvgPulseLength(0),
    fAvgPulseIntegral(0)
{
  //Default constructor
  this->SetCfgVersion(version);
}

//_______________________________________________________________________________________
CCMPMTResponse::CCMPMTResponse(const CCMPMTResponse& clufdr) 
: CCMModule(clufdr),
  fMCTruth(clufdr.fMCTruth),
  fPulses(clufdr.fPulses),
  fRD(),
  fMT(),
  fUniform(0,1),
  fSPEWeights(clufdr.fSPEWeights),
  fWaveforms(clufdr.fWaveforms),
  fWaveformsHist(clufdr.fWaveformsHist),
  fFile(clufdr.fFile),
  fHighEnergy(clufdr.fHighEnergy),
  fLowEnergy(clufdr.fLowEnergy),
  fQE(clufdr.fQE),
  fSquareWF(clufdr.fSquareWF),
  fTriangleWF(clufdr.fTriangleWF),
  fTotalHits(clufdr.fTotalHits),
  fAfterQE(clufdr.fAfterQE),
  fAfterADC(clufdr.fAfterADC),
  fTotalPulses(clufdr.fTotalPulses),
  fAvgPulseLength(clufdr.fAvgPulseLength),
  fAvgPulseIntegral(clufdr.fAvgPulseIntegral)
{
  // copy constructor
}

//_______________________________________________________________________________________
CCMPMTResponse::~CCMPMTResponse()
{ 
  // destructor
}

//_______________________________________________________________________________________
CCMResult_t CCMPMTResponse::ProcessTrigger()
{
  if (MsgLog::GetGlobalDebugLevel() >= 1) {
    MsgDebug(1,"Starting PMTResponse for MC \"Event\"");
  }

  // first get the number of hits observed by the MC event
  const size_t kNumHits = fMCTruth->GetHitNumber();

  // if the number of hits is equal to 0 then reset the Pulses object
  // and return as no pulses will be found
  if (kNumHits == 0) {
    fPulses->Reset();
    return kCCMSuccess;
  }

  // if the SPEWeights object has not been filled fill it
  // if it has been filled, reset the waveforms
  if (fSPEWeights.empty()) {
    FillSPEWeights();
    fFile = TFile::Open("check_pmtresponse_waveforms.root","RECREATE");
  } else {
    for (auto & p : fWaveforms) {
      std::fill(p.second.begin(),p.second.end(),0.0);
    }
    fFile->cd();
    for (auto & p : fWaveformsHist) {
      p.second->Write();
      p.second->Reset("ICESM");
    }
  }

  // variables used in the for loop
  int row = 0;
  int col = 0;
  int key = 0;
  double energy = 0.0;
  double time = 0.0;
  double angle = 0.0;
  double adc = 0.0;
  double length = 0.0;
  bool passedQF = true;

  fTotalHits += kNumHits;

  // loop over the hits and build the waveform
  for (size_t hit = 0; hit < kNumHits; ++hit) {
    row = fMCTruth->GetHitRow(hit);
    col = fMCTruth->GetHitCol(hit);

    // if the key value is less than 0 then the PMT is currently
    // not in the master map
    // A mapping shoule be created for the geometry
    // that is being simulated
    key = PMTInfoMap::ConvertRowColToKey(row,col);
    if (!PMTInfoMap::IsActive(key) || key < 0) {
      continue;
    }

    // NOTE: will be added in a later version
    //passedQF = fMCTruth->GetPassedQF(hit);
    if (!passedQF) {
      continue;
    }

    time = fMCTruth->GetHitTime(hit);
    energy = fMCTruth->GetHitEnergy(hit);
    angle = fMCTruth->GetHitAngle(hit);

    if (!PMTQE(energy,angle)) {
      continue;
    }

    ++fAfterQE;

    GetADCValueAndLength(key,adc,length);

    if (adc <= 0 || length <= 0) {
      continue;
    }

    ++fAfterADC;

    AddToWaveform(key,time,adc,length);

  } // end for size_t hit


  // use the waveform that was built and
  // find the pulses
  fPulses->Reset();
  for (auto & p : fWaveforms) {
    fPulses->SetKey(p.first);
    fPulses->MCFilter(p.second,0.0);
  }

  fTotalPulses += fPulses->GetNumPulses();
  for (size_t pulse = 0; pulse < fPulses->GetNumPulses(); ++pulse) {
    fAvgPulseLength += fPulses->GetPulseLength(pulse);
    fAvgPulseIntegral += fPulses->GetPulseIntegral(pulse);
  }

  return kCCMSuccess;
}

//_______________________________________________________________________________________
void CCMPMTResponse::Configure(const CCMConfig& c ) 
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

  c("QE").Get(fQE);
  c("HighEnergy").Get(fHighEnergy);
  c("LowEnergy").Get(fLowEnergy);

  std::string tempString = "";
  c("WaveformRepresentation").Get(tempString);

  for (auto & c : tempString) {
    c = toupper(c);
  }
  if (tempString.find("DEFAULT") != std::string::npos ||
      tempString.find("TRIANGLE") != std::string::npos) {
    fSquareWF = false;
    fTriangleWF = true;
  } else if (tempString.find("SQUARE") != std::string::npos) {
    fSquareWF = true;
    fTriangleWF = false;
  }

  MsgInfo(MsgLog::Form("\t-Seed values %d fMT output %u",seed,fMT()));
  MsgInfo(MsgLog::Form("\t-QE: %0.2f",fQE));
  MsgInfo(MsgLog::Form("\t-HighEnergy: %0.2f",fHighEnergy));
  MsgInfo(MsgLog::Form("\t-LowEnegy: %0.2f",fLowEnergy));
  MsgInfo(MsgLog::Form("\t-Waveform Rep. Triangle %s",(fTriangleWF) ? "true" : "false"));
  MsgInfo(MsgLog::Form("\t-Waveform Rep. Square %s",(fSquareWF) ? "true" : "false"));

  fIsInit = true;
}

//_______________________________________________________________________________________
CCMResult_t CCMPMTResponse::EndOfJob() 
{ 
  MsgInfo(MsgLog::Form("Total Hits %.0f",fTotalHits));
  MsgInfo(MsgLog::Form("After QE %.0f",fAfterQE));
  MsgInfo(MsgLog::Form("After ADC %.0f",fAfterADC));
  MsgInfo(MsgLog::Form("Total Pulses %.0f",fTotalPulses));
  MsgInfo(MsgLog::Form("Avg Pulse Length %.3f",fAvgPulseLength*Utility::fgkBinWidth/fTotalPulses));
  MsgInfo(MsgLog::Form("Avg Pulse Integral %.3f",fAvgPulseIntegral/fTotalPulses));
  MsgInfo(MsgLog::Form("Average Number Pulses Per Active PMT %.3f",fTotalPulses/static_cast<double>(fSPEWeights.size())));

  if (fFile) {
    delete fFile;
  }

  return kCCMSuccess;
}

//_______________________________________________________________________________________
bool CCMPMTResponse::PMTQE(double energy, double angle)
{
  if (energy > fHighEnergy || energy < fLowEnergy) {
    return false;
  }

  //the following will sample [0,1)
  //if you want [0,1] please let me know
  double testval = fUniform(fMT); 

  if (angle < 0) {
    angle = 0.0;
  }

  //0.2 is the QE
  double prob = fQE*std::sqrt(angle);

  if (testval < prob) {
    return true;
  }

  return false;
}

//_______________________________________________________________________________________
void CCMPMTResponse::FillSPEWeights()
{
  const size_t kMinKey = PMTInfoMap::GetMinKey();
  const size_t kMaxKey = PMTInfoMap::GetMaxKey();

  for (size_t key = kMinKey; key < kMaxKey; ++key) {
    if (!PMTInfoMap::IsActive(key)) {
      continue;
    }

    auto pmt = PMTInfoMap::GetPMTInfo(key);
    if (!pmt) {
      continue;
    }

    if (pmt->IsVeto()) {
      continue;
    }

    double adcToPE = pmt->GetADCToPE();
    double error = pmt->GetADCToPERMS();

    fSPEWeights.emplace(key,std::make_pair(adcToPE,error));
    fWaveforms.emplace(key,std::vector<double>(Utility::fgkNumBins,0.0));
    fWaveformsHist.emplace(key,std::make_shared<TH1D>(Form("hist%zu",key),"",Utility::fgkNumBins,
          0,Utility::fgkNumBins*Utility::fgkBinWidth));
  } // end for size_t key

  return;
}

//_______________________________________________________________________________________
void CCMPMTResponse::GetADCValueAndLength(size_t key, double & adc, double & length)
{
  auto itSPEWeight = fSPEWeights.find(key);
  if (itSPEWeight == fSPEWeights.end()) {
    adc = 0.0;
    length = 0.0;
    return;
  }

  std::array<double,3> energyLengthPars = {1.5609,5.77514,0.0997306};
  std::array<std::array<double,3>,3> cholDecomp;
  cholDecomp[0][0] = 0.03092;
  cholDecomp[0][1] = 0.02219;
  cholDecomp[0][2] = 0.002114;
  cholDecomp[1][0] = 0;
  cholDecomp[1][1] = 0.04763;
  cholDecomp[1][2] = 0.0007355;
  cholDecomp[2][0] = 0;
  cholDecomp[2][1] = 0;
  cholDecomp[2][2] = 0.000387;

  double adcToPE = itSPEWeight->second.first;
  double error = itSPEWeight->second.second;
  std::normal_distribution<double> gaus(adcToPE,error);

  adc = gaus(fMT);
  length = adc;

  std::normal_distribution<double> normal(0,1);
  std::array<double,3> random;
  for (int r = 0; r < 3; ++r) {
    random[r] = normal(fMT);
  }
  for (int col=0; col < 3; ++col) {
    double value = 0;
    for (int row=0; row < 3; ++row) {
      value += cholDecomp[row][col]*random[row];
    }
    energyLengthPars[col] += value;
  }

  length *= energyLengthPars[0] +
    energyLengthPars[1]*std::exp(-energyLengthPars[2]*adc);

}

//_______________________________________________________________________________________
void CCMPMTResponse::AddToWaveform(int key, double time, double adc, double length)
{
  auto itWaveform = fWaveforms.find(key);
  if (itWaveform == fWaveforms.end()) {
    return;
  }

  auto itWaveformHist = fWaveformsHist.find(key);

  // time and length must be in ns
  double startBin = std::floor(time / Utility::fgkBinWidth);
  double numBins = std::ceil(length / Utility::fgkBinWidth);
  double lastBin = std::min(startBin+numBins,static_cast<double>(Utility::fgkNumBins));

  if (fSquareWF && !fTriangleWF) {
    double adcPerBin = adc/static_cast<double>(numBins);

    for (int bin = startBin; bin < lastBin; ++bin) {
      itWaveform->second.at(bin) += adcPerBin;
    }

  } else if (fTriangleWF && !fSquareWF) {
    double middle = startBin+numBins/2.0;
    double integral = 0.0;

    for (int bin = startBin; bin < lastBin; ++bin) {
      integral = 0.0;
      if (bin < middle) {
        integral = 2.0*adc/numBins*(static_cast<double>(bin) - startBin)/(middle - startBin);
      } else if (bin != middle) {
        integral = 2.0*adc/numBins*(startBin+numBins - static_cast<double>(bin))/(startBin + numBins - middle);
      } else {
        integral = 2.0*adc/numBins;
      }
      itWaveform->second.at(bin) += integral;
      itWaveformHist->second->AddBinContent(bin+1,integral);
    } // end for bin = startBin 
  } // end if which type of waveform

  return;

}

