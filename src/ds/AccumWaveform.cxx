/*!**********************************************
 * \file AccumWaveform.cxx
 * \brief Source code for the #AccumWaveform class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "AccumWaveform.h"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <utility>

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

//-------------------------------------------------------------------------------------------------
void AccumWaveform::AddMethod(CCMAccumWaveformMethod_t method)
{
  if (fIntegralTime.find(method) != fIntegralTime.end()) {
    return;
  }

  fIntegralTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fIntegralDer.emplace(method,std::array<float,Utility::fgkNumBins>());
  fPulsesTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoBottomTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoTopTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoCLeftTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoCRightTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoCBackTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoCFrontTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fVetoTotalTime.emplace(method,std::array<float,Utility::fgkNumBins>());
  fPMTWaveform.emplace(method,std::vector<std::array<float,Utility::fgkNumBins>>());
  fPMTWaveformCount.emplace(method,std::vector<std::array<float,Utility::fgkNumBins>>());

  fIntegralTime[method].fill(0.f);
  fIntegralDer[method].fill(0.f);
  fPulsesTime[method].fill(0.f);
  fVetoBottomTime[method].fill(0.f);
  fVetoTopTime[method].fill(0.f);
  fVetoCLeftTime[method].fill(0.f);
  fVetoCRightTime[method].fill(0.f);
  fVetoCBackTime[method].fill(0.f);
  fVetoCFrontTime[method].fill(0.f);
  fVetoTotalTime[method].fill(0.f);

  for (size_t pmt = 0; pmt < Utility::fgkNumPMTs; ++pmt) {
    fPMTWaveform[method].push_back(std::array<float,Utility::fgkNumBins>());
    fPMTWaveformCount[method].push_back(std::array<float,Utility::fgkNumBins>());

    fPMTWaveform[method].back().fill(0.f);
    fPMTWaveformCount[method].back().fill(0.f);
  }
}

//-------------------------------------------------------------------------------------------------
std::array<float,Utility::fgkNumBins> * AccumWaveform::Get(CCMAccumWaveformMethod_t method, 
    CCMAccWaveform_t waveformType, int pmtID)
{
  if (method == kCCMAccumWaveformTotalID) {
    MsgFatal("Cannot get kCCMAccumWaveformTotalID waveform as it does not exists");
  }

  auto it = fPulsesTime.find(method);
  if (it == fPulsesTime.end()) {
    AddMethod(method);
  }

  std::array<float,Utility::fgkNumBins> * tempA = nullptr;
  switch (waveformType) {
    case kCCMPulsesTimeID:       tempA = &(fPulsesTime.find(method)->second); break;
    case kCCMIntegralTimeID:     tempA = &(fIntegralTime.find(method)->second); break;
    case kCCMIntegralDerID:      tempA = &(fIntegralDer.find(method)->second); break;
    case kCCMVetoBottomTimeID:   tempA = &(fVetoBottomTime.find(method)->second); break;
    case kCCMVetoTopTimeID:      tempA = &(fVetoTopTime.find(method)->second); break;
    case kCCMVetoCRightTimeID:   tempA = &(fVetoCRightTime.find(method)->second); break;
    case kCCMVetoCLeftTimeID:    tempA = &(fVetoCLeftTime.find(method)->second); break;
    case kCCMVetoCFrontTimeID:   tempA = &(fVetoCFrontTime.find(method)->second); break;
    case kCCMVetoCBackTimeID:    tempA = &(fVetoCBackTime.find(method)->second); break;
    case kCCMVetoTotalTimeID:    tempA = &(fVetoTotalTime.find(method)->second); break;
    case kCCMPMTWaveformID:      tempA = &(fPMTWaveform.find(method)->second.at(pmtID)); break;
    case kCCMPMTWaveformCountID: tempA = &(fPMTWaveformCount.find(method)->second.at(pmtID)); break;
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
          Utility::ConvertCCMAccumWaveformMethodToString(p.first).c_str(),
          std::accumulate(p.second.begin(),p.second.end(),0.f)));
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

