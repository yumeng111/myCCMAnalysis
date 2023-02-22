/*!**********************************************
 * \file CCMPMTResponse.cxx
 * \author R.T. Thornton (LANL), E. Dunton (Columbia)
 * \date November 19, 2020
 *
 * Main code to modle the response of the PMTs
 * in the simulation.
 ***********************************************/

#include <array>
#include <cctype>

#include "CCMAnalysis/CCMSimulationUtils/CCMPMTResponse.h"

#include "CCMAnalysis/CCMDataStructures/Pulses.h"
#include "CCMAnalysis/CCMDataStructures/MCTruth.h"
#include "CCMAnalysis/CCMDataStructures/SinglePulse.h"
#include "CCMAnalysis/CCMFramework/CCMConfig.h"
#include "CCMAnalysis/CCMFramework/CCMConfigParam.h"
#include "CCMAnalysis/CCMFramework/CCMModuleTable.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMUtils/PMTInfoMap.h"
#include "CCMAnalysis/CCMUtils/PMTInformation.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"

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
    fTRandom(std::make_shared<TRandom3>(0)),
    fSPEWeights(),
    fSPEHists(),
    fWaveforms(),
    fPMTRelEff(),
    //fWaveformsHist(),
    fPMTSPEFile(nullptr),
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
  fTRandom(clufdr.fTRandom),
  fSPEWeights(clufdr.fSPEWeights),
  fSPEHists(clufdr.fSPEHists),
  fWaveforms(clufdr.fWaveforms),
  fPMTRelEff(clufdr.fPMTRelEff),
  //fWaveformsHist(clufdr.fWaveformsHist),
  fPMTSPEFile(clufdr.fPMTSPEFile),
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
    return kCCMFailure;
  }

  // if the SPEWeights object has not been filled fill it
  // if it has been filled, reset the waveforms
  if (fSPEWeights.empty() && fSPEHists.empty()) {
    FillSPEWeights();
  } else {
    for (auto & p : fWaveforms) {
      std::fill(p.second.begin(),p.second.end(),0.0);
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

  fTotalHits += kNumHits;

  // loop over the hits and build the waveform
  for (size_t hit = 0; hit < kNumHits; ++hit) {
    if (!fMCTruth->GetPassedQF(hit)) {
      continue;
    }

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

    time = fMCTruth->GetHitTime(hit);
    energy = fMCTruth->GetHitEnergy(hit);
    angle = fMCTruth->GetHitAngle(hit);

    if (!PMTQE(key,energy,angle)) {
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

  if (fPulses->GetNumPulses() == 0) {
    //MsgWarning("No Pulses Found");
    return kCCMFailure;
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
    fTRandom->SetSeed(0);
  } else {
    if (seed > 0) {
      fMT.seed(seed);
      fTRandom->SetSeed(seed);
    } else {
      seed = std::abs(seed);
      std::string env = std::getenv("CCMPROJECT");
      std::string fileName = env + "/calibrationFiles/2019/listOfRandomSeeds.txt";
      std::ifstream infile(fileName.c_str());
      std::vector<unsigned int> seeds((std::istream_iterator<unsigned int>(infile)),
                                       std::istream_iterator<unsigned int>());
      infile.close();
      if (static_cast<unsigned int>(seed) < seeds.size()) {
        fMT.seed(seeds.at(seed));
        fTRandom->SetSeed(seeds.at(seed));
        seed = seeds.at(seed);
      } else {
        MsgFatal(MsgLog::Form("Size of seeds list %zu and requested seed is %d",seeds.size(),seed));
      }
    }
  }
  gRandom = fTRandom.get();

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

  c("PMTSPEFile").Get(tempString);
  if (!tempString.empty()) {
    fPMTSPEFile = TFile::Open(tempString.c_str(),"READ");
  }

  c("PMTRelEff").Get(tempString);
  if (!tempString.empty()) {
    MsgInfo(MsgLog::Form("\t-Getting PMTRelEff from %s",tempString.c_str()));
    std::ifstream file(tempString.c_str());
    std::copy(std::istream_iterator<std::pair<int,double>>(file),
              std::istream_iterator<std::pair<int,double>>(),
              std::insert_iterator<std::map<int,double>>(fPMTRelEff,fPMTRelEff.begin()));
    file.close();
    MsgInfo(MsgLog::Form("\t\t-fPMTRelEff size = %zu",fPMTRelEff.size()));
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
  MsgInfo(MsgLog::Form("Average Number Pulses Per Active PMT %.3f",
        fTotalPulses/static_cast<double>(std::max(fSPEWeights.size(),fSPEHists.size()))));

  //if (fFile) {
  //  delete fFile;
  //}

  return kCCMSuccess;
}

//_______________________________________________________________________________________
bool CCMPMTResponse::PMTQE(size_t key, double energy, double angle)
{
  if (energy > fHighEnergy || energy < fLowEnergy) {
    return false;
  }

  //the following will sample [0,1)
  //if you want [0,1] please let me know

  if (!fPMTRelEff.empty()) {
    auto relEff = fPMTRelEff.find(key);
    if (relEff == fPMTRelEff.end()) {
      return false;
    }

    // since the maximum relEff is equal to 1
    // we want to throw away PE if a random number
    // is greater than the relEff. So the PMT
    // that is seen to have the maximum gain
    // will have no hits thrown away since it
    // has a relEff of 1
    double testval = fUniform(fMT); 
    if (testval > relEff->second) {
      return false;
    }
  }

  double testval = fUniform(fMT);
  if (angle < 0) {
    angle = 0.0;
  }

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

    if (!fPMTSPEFile) {
      double adcToPE = pmt->GetADCToPE();
      double error = pmt->GetADCToPERMS();
      double thresh = pmt->GetADCThreshold();

      fSPEWeights.emplace(key,std::make_tuple(adcToPE,error,thresh));
      fWaveforms.emplace(key,std::vector<double>(Utility::fgkNumBins,0.0));
    } else {
      TH2D * tempHist = nullptr;
      fPMTSPEFile->GetObject(Form("subtract%zu",key),tempHist);
      if (!tempHist) {
        continue;
      }
      fSPEHists.emplace(key,std::make_shared<TH2D>(*tempHist));
      fWaveforms.emplace(key,std::vector<double>(Utility::fgkNumBins,0.0));
      delete tempHist;
    }
  } // end for size_t key

  if (fPMTSPEFile) {
    delete fPMTSPEFile;
  }

  return;
}

//_______________________________________________________________________________________
void CCMPMTResponse::GetADCValueAndLength(size_t key, double & adc, double & length)
{
  if (fSPEHists.empty()) {
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

    double adcToPE = std::get<0>(itSPEWeight->second);
    double error = std::get<1>(itSPEWeight->second);
    double thresh = std::get<2>(itSPEWeight->second);
    std::normal_distribution<double> gaus(adcToPE,error);

    adc = gaus(fMT);
    while (adc < thresh) {
      adc = gaus(fMT);
    }
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

    return;
  }

  auto itSPEHist = fSPEHists.find(key);
  if (itSPEHist == fSPEHists.end()) {
    adc = 0.0;
    length = 0.0;
    return;
  }
  itSPEHist->second->GetRandom2(adc,length);
}

//_______________________________________________________________________________________
void CCMPMTResponse::AddToWaveform(int key, double time, double adc, double length)
{
  auto itWaveform = fWaveforms.find(key);
  if (itWaveform == fWaveforms.end()) {
    return;
  }

  //auto itWaveformHist = fWaveformsHist.find(key);

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
      //itWaveformHist->second->AddBinContent(bin+1,integral);
    } // end for bin = startBin 
  } // end if which type of waveform

  return;

}

