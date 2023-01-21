/*!**********************************************
 * \file AccumWaveform.cxx
 * \brief Source code for the #AccumWaveform class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/

#include "CCMAnalysis/CCMDataStructures/AccumWaveform.h"

#ifndef __CINT__
#include <icetray/I3Logging.h>
#include <icetray/serialization.h>
#endif // __CINT__

#include <map>
#include <cmath>
#include <limits>
#include <utility>
#include <iterator>
#include <algorithm>

ClassImp(AccumWaveform)

/*!**********************************************
 * \fn AccumWaveform::AccumWaveform()
 * \brief Default constructor
 ***********************************************/
AccumWaveform::AccumWaveform() : TObject()
{
  fEventNumber = 0;
  fComputerSecIntoEpoch = 0;
  fComputerNSIntoSec = 0;
  fBeamTime = 0;
  fBeamIntegral= 0;
  fBeamLength= 0;

  Reset();
}

/*!**********************************************
 * \fn AccumWaveform::AccumWaveform(const AccumWaveform & p)
 * \brief Copy constructor calles #operator=
 * \param[in] p The object to copy
 ***********************************************/
AccumWaveform::AccumWaveform(const AccumWaveform & p) : TObject(p)
{
  this->operator=(p);
}

/*!**********************************************
 * \fn AccumWaveform::~AccumWaveform()
 * \brief Destructor calles #Reset
 ***********************************************/
AccumWaveform::~AccumWaveform()
{
  Reset();
}

/*!**********************************************
 * \fn void AccumWaveform::Reset()
 * \brief Resets all the variables, vectors and arrays. Calls #ClearAccumWaveform
 ***********************************************/
void AccumWaveform::Reset()
{
  fTriggerTime = 0.0;
  ClearAccumWaveform();
}

/*!**********************************************
 * \fn void AccumWaveform::ClearAccumWaveform()
 * \brief Resets vectors
 ***********************************************/
void AccumWaveform::ClearAccumWaveform()
{
  for (auto & p : fIntegralTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fIntegralDer) {
    p.second.fill(0.f);
  }
  for (auto & p : fPulsesTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoBottomTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoTopTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoCLeftTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoCRightTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoCBackTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoCFrontTime) {
    p.second.fill(0.f);
  }
  for (auto & p : fVetoTotalTime) {
    p.second.fill(0.f);
  }

  for (auto & p: fPMTWaveform) {
    for (auto & a : p.second) {
      a.fill(0.f);
    }
  }
  for (auto & p: fPMTWaveform) {
    for (auto & a : p.second) {
      a.fill(0.f);
    }
  }
}


/*!**********************************************
 * \fn AccumWaveform & AccumWaveform::operator=(const AccumWaveform & rhs) 
 * \brief Set the current #AccumWaveform object to the one passed
 * \param[in] rhs The object to copy
 ***********************************************/
AccumWaveform & AccumWaveform::operator=(const AccumWaveform & rhs) 
{
  this->fEventNumber = rhs.fEventNumber;
  this->fComputerSecIntoEpoch= rhs.fComputerSecIntoEpoch;
  this->fComputerNSIntoSec= rhs.fComputerNSIntoSec;
  this->fBeamTime= rhs.fBeamTime;
  this->fBeamIntegral= rhs.fBeamIntegral;
  this->fBeamLength= rhs.fBeamLength;
  this->fTriggerTime = rhs.fTriggerTime;
  this->fPulsesTime = rhs.fPulsesTime;
  this->fIntegralTime = rhs.fIntegralTime;
  this->fIntegralDer = rhs.fIntegralDer;

  this->fVetoBottomTime = rhs.fVetoBottomTime;
  this->fVetoTopTime = rhs.fVetoTopTime;
  this->fVetoCRightTime = rhs.fVetoCRightTime;
  this->fVetoCLeftTime = rhs.fVetoCLeftTime;
  this->fVetoCFrontTime = rhs.fVetoCFrontTime;
  this->fVetoCBackTime = rhs.fVetoCBackTime;
  this->fVetoTotalTime = rhs.fVetoTotalTime;

  this->fPMTWaveform = rhs.fPMTWaveform;
  this->fPMTWaveformCount = rhs.fPMTWaveformCount;

  return *this;
}

/*!**********************************************
 * \fn AccumWaveform & AccumWaveform::operator+=(const AccumWaveform & rhs) 
 * \brief Add the passed #AccumWaveform object to the current one
 * \param[in] rhs The object to add from
 ***********************************************/
AccumWaveform & AccumWaveform::operator+=(const AccumWaveform & rhs) 
{
  AddMaps(this->fPulsesTime,rhs.fPulsesTime);
  AddMaps(this->fIntegralTime,rhs.fIntegralTime);
  AddMaps(this->fIntegralDer,rhs.fIntegralDer);
  AddMaps(this->fVetoBottomTime,rhs.fVetoBottomTime);
  AddMaps(this->fVetoTopTime,rhs.fVetoTopTime);
  AddMaps(this->fVetoCRightTime,rhs.fVetoCRightTime);
  AddMaps(this->fVetoCLeftTime,rhs.fVetoCLeftTime);
  AddMaps(this->fVetoCFrontTime,rhs.fVetoCFrontTime);
  AddMaps(this->fVetoCBackTime,rhs.fVetoCBackTime);
  AddMaps(this->fVetoTotalTime,rhs.fVetoTotalTime);
  AddMaps(this->fPMTWaveform,rhs.fPMTWaveform);
  AddMaps(this->fPMTWaveformCount,rhs.fPMTWaveformCount);

  return *this;

}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::AddMaps(MapDAQWF1D & lhs, const MapDAQWF1D & rhs)
{
  for (auto & p : lhs) {
    auto rhsIt = rhs.find(p.first);
    if (rhsIt == rhs.end()) {
      continue;
    }
    std::transform(p.second.begin(),p.second.end(),
        rhsIt->second.begin(),p.second.begin(),std::plus<float>());
  }
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::AddMaps(MapDAQWF2D & lhs, const MapDAQWF2D & rhs)
{
  for (auto & p : lhs) {
    auto rhsIt = rhs.find(p.first);
    if (rhsIt == rhs.end()) {
      continue;
    }
    const size_t kNumPMTs = p.second.size();
    const size_t kNumRhs = rhsIt->second.size();
    
    size_t limit = kNumPMTs;
    
    if (kNumPMTs > kNumRhs) {
      limit = kNumRhs;
    }
    
    //std::cout << kNumPMTs << '\t' << rhsIt->second.size() << std::endl;
    for (size_t pmt = 0; pmt < limit; ++pmt) {
      std::transform(p.second.at(pmt).begin(),p.second.at(pmt).end(),
          rhsIt->second.at(pmt).begin(),p.second.at(pmt).begin(),std::plus<float>());
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::AddMethod(CCMAccumWaveformMethod_t method)
{
  int methodInt = static_cast<int>(method);
  if (fIntegralTime.find(methodInt) != fIntegralTime.end()) {
    return;
  }

  fIntegralTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fIntegralDer.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fPulsesTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoBottomTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoTopTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoCLeftTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoCRightTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoCBackTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoCFrontTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fVetoTotalTime.emplace(methodInt,std::array<float,Utility::fgkNumBins>());
  fPMTWaveform.emplace(methodInt,std::vector<std::array<float,Utility::fgkNumBins>>());
  fPMTWaveformCount.emplace(methodInt,std::vector<std::array<float,Utility::fgkNumBins>>());

  fIntegralTime[methodInt].fill(0.f);
  fIntegralDer[methodInt].fill(0.f);
  fPulsesTime[methodInt].fill(0.f);
  fVetoBottomTime[methodInt].fill(0.f);
  fVetoTopTime[methodInt].fill(0.f);
  fVetoCLeftTime[methodInt].fill(0.f);
  fVetoCRightTime[methodInt].fill(0.f);
  fVetoCBackTime[methodInt].fill(0.f);
  fVetoCFrontTime[methodInt].fill(0.f);
  fVetoTotalTime[methodInt].fill(0.f);

  for (size_t pmt = 0; pmt < Utility::fgkNumPMTs; ++pmt) {
    fPMTWaveform[methodInt].push_back(std::array<float,Utility::fgkNumBins>());
    fPMTWaveformCount[methodInt].push_back(std::array<float,Utility::fgkNumBins>());

    fPMTWaveform[methodInt].back().fill(0.f);
    fPMTWaveformCount[methodInt].back().fill(0.f);
  }
}

//-------------------------------------------------------------------------------------------------
std::array<float,Utility::fgkNumBins> * AccumWaveform::Get(CCMAccumWaveformMethod_t method, 
    CCMAccWaveform_t waveformType, int pmtID)
{
  if (method == kCCMAccumWaveformTotalID) {
    MsgFatal("Cannot get kCCMAccumWaveformTotalID waveform as it does not exists");
  }

  int methodInt = static_cast<int>(method);
  auto it = fPulsesTime.find(methodInt);
  if (it == fPulsesTime.end()) {
    AddMethod(method);
  }

  std::array<float,Utility::fgkNumBins> * tempA = nullptr;
  switch (waveformType) {
    case kCCMPulsesTimeID:       tempA = &(fPulsesTime.find(methodInt)->second); break;
    case kCCMIntegralTimeID:     tempA = &(fIntegralTime.find(methodInt)->second); break;
    case kCCMIntegralDerID:      tempA = &(fIntegralDer.find(methodInt)->second); break;
    case kCCMVetoBottomTimeID:   tempA = &(fVetoBottomTime.find(methodInt)->second); break;
    case kCCMVetoTopTimeID:      tempA = &(fVetoTopTime.find(methodInt)->second); break;
    case kCCMVetoCRightTimeID:   tempA = &(fVetoCRightTime.find(methodInt)->second); break;
    case kCCMVetoCLeftTimeID:    tempA = &(fVetoCLeftTime.find(methodInt)->second); break;
    case kCCMVetoCFrontTimeID:   tempA = &(fVetoCFrontTime.find(methodInt)->second); break;
    case kCCMVetoCBackTimeID:    tempA = &(fVetoCBackTime.find(methodInt)->second); break;
    case kCCMVetoTotalTimeID:    tempA = &(fVetoTotalTime.find(methodInt)->second); break;
    case kCCMPMTWaveformID:      tempA = &(fPMTWaveform.find(methodInt)->second.at(pmtID)); break;
    case kCCMPMTWaveformCountID: tempA = &(fPMTWaveformCount.find(methodInt)->second.at(pmtID)); break;
    default: break;
  }

  if (tempA == nullptr) {
    MsgFatal("CCMAccumWaveform_t could not be found");
  }

  return tempA;
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::SetIndex(size_t index, float weight, CCMAccumWaveformMethod_t method, 
    CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  waveform->at(index) = weight;

  return;
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::FillIndex(size_t index, float weight, CCMAccumWaveformMethod_t method, 
    CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  waveform->at(index) += weight;

  return;
}

//-------------------------------------------------------------------------------------------------
float AccumWaveform::GetIndex(size_t index, CCMAccumWaveformMethod_t method, 
    CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  return waveform->at(index);
}

//-------------------------------------------------------------------------------------------------
float AccumWaveform::Integrate(size_t start, size_t end, CCMAccumWaveformMethod_t method, 
    CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  return std::accumulate(waveform->begin()+start,waveform->begin()+end,0.f);
}

//-------------------------------------------------------------------------------------------------
int AccumWaveform::FindFirstNoneEmptyBin(size_t start, size_t end, 
    CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  return Utility::FindFirstNoneEmptyBin<float>(waveform->begin(),waveform->begin()+start,waveform->begin()+end);
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::DumpInfo()
{
  MsgInfo(MsgLog::Form("Beam Time %d Integral %.2f Length %.2f",fBeamTime,fBeamIntegral,fBeamLength));
  for (auto & p : fIntegralTime) {
    MsgInfo(MsgLog::Form("method %s integralTime integral %f",
          Utility::ConvertCCMAccumWaveformMethodToString(static_cast<CCMAccumWaveformMethod_t>(p.first)).c_str(),
          std::accumulate(p.second.begin(),p.second.end(),0.f)));
  }
  for (auto & p : fPMTWaveform) {
    MsgInfo(MsgLog::Form("method %s PMTWaveform size %zu",
          Utility::ConvertCCMAccumWaveformMethodToString(static_cast<CCMAccumWaveformMethod_t>(p.first)).c_str(),
          p.second.size()));
  }
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::Max(size_t & loc, float & value, size_t start, size_t end, 
    CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  if (start >= waveform->size() || end >= waveform->size()) {
    MsgFatal(MsgLog::Form("Start (%zu) or end (%zu) is out of range (%zu)",start,end,waveform->size()));
  }
  auto it = std::max_element(waveform->begin()+start, waveform->begin()+end);
  loc = std::distance(waveform->begin(),it);
  value = *it;
  return;
}

//-------------------------------------------------------------------------------------------------
void AccumWaveform::Min(size_t & loc, float & value, size_t start, size_t end,
    CCMAccumWaveformMethod_t method, CCMAccWaveform_t waveformType, int pmtID)
{
  auto waveform = Get(method,waveformType,pmtID);
  auto it = std::min_element(waveform->begin()+start, waveform->begin()+end);
  loc = std::distance(waveform->begin(),it);
  value = *it;
  return;
}

#ifndef __CINT__
template <class Archive>
void
AccumWaveform::serialize(Archive& ar, unsigned version) {
    if (version > legacy_accum_waveform_version_)
        log_fatal("Attempting to read version %u from file but running version %u of AccumWaveform class.", version, legacy_accum_waveform_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("fEventNumber", fEventNumber);
    ar & make_nvp("fComputerSecIntoEpoch", fComputerSecIntoEpoch);
    ar & make_nvp("fComputerNSIntoSec", fComputerNSIntoSec);

    ar & make_nvp("fBeamTime", fBeamTime);
    ar & make_nvp("fBeamIntegral", fBeamIntegral);
    ar & make_nvp("fBeamLength", fBeamLength);

    ar & make_nvp("fTriggerTime", fTriggerTime);

    ar & make_nvp("fPulsesTime", fPulsesTime);
    ar & make_nvp("fIntegralTime", fIntegralTime);
    ar & make_nvp("fIntegralDer", fIntegralDer);

    ar & make_nvp("fVetoBottomTime", fVetoBottomTime);
    ar & make_nvp("fVetoTopTime", fVetoTopTime);
    ar & make_nvp("fVetoCRightTime", fVetoCRightTime);
    ar & make_nvp("fVetoCLeftTime", fVetoCLeftTime);
    ar & make_nvp("fVetoCFrontTime", fVetoCFrontTime);
    ar & make_nvp("fVetoCBackTime", fVetoCBackTime);
    ar & make_nvp("fVetoTotalTime", fVetoTotalTime);

    ar & make_nvp("fPMTWaveform", fPMTWaveform);
    ar & make_nvp("fPMTWaveformCount", fPMTWaveformCount);
}

I3_SERIALIZABLE(AccumWaveform, LegacyAccumWaveform);
#endif // __CINT__
