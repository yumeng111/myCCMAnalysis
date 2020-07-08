/*!**********************************************
 * \file SimplifiedEvent.cxx
 * \brief Source code for the #SimplifiedEvent class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "SimplifiedEvent.h"
#include "MsgLog.h"

#include <iostream>
#include <cmath>
#include <iterator>

ClassImp(SimplifiedEvent)

//-------------------------------------------------------------------------------------------------
SimplifiedEvent::SimplifiedEvent() : TObject()
{
  Reset();
}

//-------------------------------------------------------------------------------------------------
SimplifiedEvent::SimplifiedEvent(const SimplifiedEvent & rhs) : TObject(rhs)
{
  fEventFinderMethod       = rhs.fEventFinderMethod;
  fAccumWaveformMethod     = rhs.fAccumWaveformMethod;
  fThreshold               = rhs.fThreshold;
  fStartTime               = rhs.fStartTime;
  fLength                  = rhs.fLength;
  fLargestPMTFraction      = rhs.fLargestPMTFraction;
  fCoatedStartTime         = rhs.fCoatedStartTime;
  fUncoatedStartTime       = rhs.fUncoatedStartTime;
  fVetoBottomStartTime     = rhs.fVetoBottomStartTime;
  fVetoTopStartTime        = rhs.fVetoTopStartTime;
  fVetoBackStartTime       = rhs.fVetoBackStartTime;
  fVetoRightStartTime      = rhs.fVetoRightStartTime;
  fVetoLeftStartTime       = rhs.fVetoLeftStartTime;
  fVetoFrontStartTime      = rhs.fVetoFrontStartTime;
  fHitXPosition            = rhs.fHitXPosition;
  fHitYPosition            = rhs.fHitYPosition;
  fHitZPosition            = rhs.fHitZPosition;
  for (int i = 0; i <2; ++i) {
    fCoatedInt[i]          = rhs.fCoatedInt[i];
    fUncoatedInt[i]        = rhs.fUncoatedInt[i];

    fNCoated[i]            = rhs.fNCoated[i];
    fNUncoated[i]          = rhs.fNUncoated[i];
    fNVetoBottom[i]        = rhs.fNVetoBottom[i];
    fNVetoTop[i]           = rhs.fNVetoTop[i];
    fNVetoBack[i]          = rhs.fNVetoBack[i];
    fNVetoLeft[i]          = rhs.fNVetoLeft[i];
    fNVetoRight[i]         = rhs.fNVetoRight[i];
    fNVetoFront[i]         = rhs.fNVetoFront[i];
  }

  fPMTHits             = rhs.fPMTHits;
  fPMTHitsVec         = rhs.fPMTHitsVec;

  fAccWaveformCount    = rhs.fAccWaveformCount;
  fAccWaveformInt      = rhs.fAccWaveformInt;
  fPMTWaveform         = rhs.fPMTWaveform;
  fPMTWaveformItr     = rhs.fPMTWaveformItr;
}

//-------------------------------------------------------------------------------------------------
SimplifiedEvent::~SimplifiedEvent()
{
  // destructor
}

/*!**********************************************
 * \fn void SimplifiedEvent::Reset()
 * \brief Rests all the variables back to default
 ***********************************************/
void SimplifiedEvent::Reset()
{
  fEventFinderMethod       = static_cast<int>(kCCMDynamicLengthEventID);
  fAccumWaveformMethod     = static_cast<int>(kCCMAccumWaveformTotalID);

  fThreshold          = 0.f;
  fLargestPMTFraction = 1.f;


  // Time variables
  // If the variable is -1 then the corresponding multiplicity count should be equal to 0
  fStartTime               = -1;
  fLength                  = -1;
  fCoatedStartTime         = -1;
  fUncoatedStartTime       = -1;
  fVetoBottomStartTime     = -1;
  fVetoTopStartTime        = -1;
  fVetoBackStartTime       = -1;
  fVetoRightStartTime      = -1;
  fVetoLeftStartTime       = -1;
  fVetoFrontStartTime      = -1;

  // Integral variables
  // The integral is calculated for the event only if the PMT was above threshold
  // independent of the value of the integral or when the PMT crossed threshold
  fCoatedInt[0]         = 0.f;
  fUncoatedInt[0]       = 0.f;

  fCoatedInt[1]         = 0.f;
  fUncoatedInt[1]       = 0.f;

  // Position variables
  fHitXPosition = 0.f;
  fHitYPosition = 0.f;
  fHitZPosition = 0.f;

  // Multiplicity variables
  fNCoated[0]     = 0;
  fNUncoated[0]   = 0;
  fNVetoBottom[0] = 0;
  fNVetoTop[0]    = 0;
  fNVetoBack[0]   = 0;
  fNVetoLeft[0]   = 0;
  fNVetoRight[0]  = 0;
  fNVetoFront[0]  = 0;
  fNCoated[1]     = 0;
  fNUncoated[1]   = 0;
  fNVetoBottom[1] = 0;
  fNVetoTop[1]    = 0;
  fNVetoBack[1]   = 0;
  fNVetoLeft[1]   = 0;
  fNVetoRight[1]  = 0;
  fNVetoFront[1]  = 0;

  fPMTHits = 0;
  fPMTHitsVec.clear();

  if (!fAccWaveformCount.empty()) {
    fAccWaveformCount.clear();
  }
  if (!fAccWaveformInt.empty()) {
    fAccWaveformInt.clear();
  }

  if (!fPMTWaveform.empty()) {
    fPMTWaveform.empty();
  }
}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetNumTank(bool prompt) 
{
  return fNCoated[prompt] + fNUncoated[prompt];
}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetNumVeto(bool prompt) 
{
  return fNVetoBottom[prompt] + fNVetoTop[prompt] + GetNumVetoSide(prompt);
}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetNumVetoSide(bool prompt) 
{
  return fNVetoBack[prompt] + fNVetoLeft[prompt] + fNVetoFront[prompt] + fNVetoRight[prompt];
}

//-------------------------------------------------------------------------------------------------
double SimplifiedEvent::GetStartTimeTank() 
{
  if (fCoatedStartTime == -1) {
    return fUncoatedStartTime;
  }

  if (fUncoatedStartTime == -1) {
    return fCoatedStartTime;
  }

  return std::min(fCoatedStartTime,fUncoatedStartTime);
}

//-------------------------------------------------------------------------------------------------
double SimplifiedEvent::GetStartTimeVeto() 
{
  bool bottom = fNVetoBottom[0] || fNVetoBottom[1];
  bool top    = fNVetoTop[0] || fNVetoTop[1];
  bool back   = fNVetoBack[0] || fNVetoBack[1];
  bool right  = fNVetoRight[0] || fNVetoRight[1];
  bool left   = fNVetoLeft[0] || fNVetoLeft[1];
  bool front  = fNVetoFront[0] || fNVetoFront[1];

  if (!bottom && !top && !back && !right && !left && !front) {
    return -1;
  }

  double startTime = 0;

  if (bottom) {
    startTime = std::min(startTime,fVetoBottomStartTime);
  }

  if (top) {
    startTime = std::min(startTime,fVetoTopStartTime);
  }

  if (back) {
    startTime = std::min(startTime,fVetoBackStartTime);
  }

  if (right) {
    startTime = std::min(startTime,fVetoRightStartTime);
  }

  if (left) {
    startTime = std::min(startTime,fVetoLeftStartTime);
  }

  if (front) {
    startTime = std::min(startTime,fVetoFrontStartTime);
  }

  return startTime;

}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetIntegralTank(bool prompt) 
{
  return fCoatedInt[prompt] + fUncoatedInt[prompt];
}

/*!**********************************************
 * \fn void SimplifiedEvent::AddWaveforms(const std::vector<float> & count, const std::vector<float> & integral)
 * \brief Set the accumulated waveforms to what was passed
 * \param[in] count The PMT multiplicity accumulated waveform
 * \param[in] integral The integral charge accumulated waveform
 ***********************************************/
void SimplifiedEvent::AddWaveforms(const std::vector<float> & count, const std::vector<float> & integral)
{
  fAccWaveformCount = count;
  fAccWaveformInt = integral;
}

void SimplifiedEvent::AddPMTWaveforms(const int key, const std::vector<int> & count, const std::vector<float> & integral)
{
  auto iter = fPMTWaveform.find(key);
  if (iter == fPMTWaveform.end()) {
    fPMTWaveform.emplace(key,std::make_pair(integral,count));
    return;
  }
  iter->second.second = count;
  iter->second.first = integral;
}


void SimplifiedEvent::ResetPMTWaveformItr()
{
  fPMTWaveformItr = fPMTWaveform.cbegin();
}

void SimplifiedEvent::ResetPMTWaveformItr() const
{
  fPMTWaveformItr = fPMTWaveform.begin();
}

bool SimplifiedEvent::NextPMTWaveform()
{
  if (fPMTWaveformItr == fPMTWaveform.cend()) {
    return false;
  }

  std::advance(fPMTWaveformItr,1);

  if (fPMTWaveformItr == fPMTWaveform.end()) {
    return false;
  }

  return true;
}

bool SimplifiedEvent::NextPMTWaveform() const
{
  if (fPMTWaveformItr == fPMTWaveform.end()) {
    return false;
  }

  std::advance(fPMTWaveformItr,1);
  
  if (fPMTWaveformItr == fPMTWaveform.end()) {
    return false;
  }

  return true;
}

void SimplifiedEvent::GetPMTWaveform(int & key, std::vector<float> & vecInt, std::vector<int> & vecCount)
{
  key = fPMTWaveformItr->first;
  vecInt = fPMTWaveformItr->second.first;
  vecCount = fPMTWaveformItr->second.second;
}

void SimplifiedEvent::GetPMTWaveform(int & key, std::vector<float> & vecInt, std::vector<int> & vecCount) const
{
  key = fPMTWaveformItr->first;
  vecInt = fPMTWaveformItr->second.first;
  vecCount = fPMTWaveformItr->second.second;
}

