/*!**********************************************
 * \file SinglePulse.cxx
 * \brief Source code for the #SinglePulse class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/

#include <cmath>

#include "CCMAnalysis/CCMDataStructures/SinglePulse.h"

#include "CCMAnalysis/CCMUtils/MsgLog.h"

ClassImp(SinglePulse)

//------------------------------------------------------
SinglePulse::SinglePulse(int key) : TObject()
{
  fKey = key;

  Reset();

}

//------------------------------------------------------
SinglePulse::SinglePulse(const SinglePulse & p) : TObject(p)
{
  this->operator=(p);
}

//------------------------------------------------------
SinglePulse::~SinglePulse()
{
  Reset();
}

//------------------------------------------------------
void SinglePulse::Reset()
{
  fADCToPE = 1.f;
  fTrigOffset = 0.f;
  fPMTOffset = 0.f;
  fAmplitude = 0.f;
  fBaseline = 0.f;
  fMaxDerValue= 0.f;
  fIntegral = 0.f;
  fTime = -9999.f;
  fLength = 0.f;

  fWaveformStart = 0.f;
  fWaveformEnd = 0.f;
  if (!fWaveform.empty()) {
    fWaveform.clear();
  }

}

/*!**********************************************
 * \fn void SinglePulse::Append(const SinglePulse & rhs) 
 * \brief Add a pulse to the current pulse
 * \param[in] rhs The pulse to add to the current pulse
 ***********************************************/
void SinglePulse::Append(const SinglePulse & rhs) 
{
  fAmplitude = std::max(fAmplitude,rhs.fAmplitude);
  fMaxDerValue = std::max(fMaxDerValue,rhs.fMaxDerValue);
  fIntegral += rhs.fIntegral;
  fLength += rhs.fLength;

  if (fWaveformEnd < rhs.fWaveformEnd) {
    int start = fWaveformEnd - rhs.fWaveformStart + 1;
    int size = rhs.fWaveform.size();
    for (int loc = start; loc < size; ++loc) {
      fWaveform.push_back(rhs.fWaveform[loc]);
    }
    fWaveformEnd = rhs.fWaveformEnd;
  }
}

//------------------------------------------------------
SinglePulse & SinglePulse::operator=(const SinglePulse & rhs)
{
  this->fKey = rhs.fKey;
  this->fADCToPE = rhs.fADCToPE;
  this->fTrigOffset = rhs.fTrigOffset;
  this->fPMTOffset = rhs.fPMTOffset;
  this->fLength = rhs.fLength;
  this->fTime = rhs.fTime;
  this->fMaxDerValue= rhs.fMaxDerValue;
  this->fBaseline = rhs.fBaseline;
  this->fAmplitude = rhs.fAmplitude;
  this->fIntegral = rhs.fIntegral;
  this->fWaveformStart = rhs.fWaveformStart;
  this->fWaveformEnd = rhs.fWaveformEnd;
  this->fWaveform = rhs.fWaveform;

  return *this;
}

