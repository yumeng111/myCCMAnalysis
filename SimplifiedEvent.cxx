/*!**********************************************
 * \file SimplifiedEvent.cxx
 * \brief Source code for the #SimplifiedEvent class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "SimplifiedEvent.h"
#include "MsgLog.h"
//#include "BoardInfo.h"
//#include "ChannelInfo.h"

#include <iostream>
#include <cmath>

ClassImp(SimplifiedEvent)

float SimplifiedEvent::fgkStartDAQWindow = -9.92;
float SimplifiedEvent::fgkSampleWidth = 2e-3;

//-------------------------------------------------------------------------------------------------
SimplifiedEvent::SimplifiedEvent() : TObject()
{
  Reset();
}

//-------------------------------------------------------------------------------------------------
SimplifiedEvent::SimplifiedEvent(const SimplifiedEvent & rhs) : TObject(rhs)
{
  fComputerSecIntoEpoch    = rhs.fComputerSecIntoEpoch;
  fComputerNSIntoSec       = rhs.fComputerNSIntoSec;
  fThreshold               = rhs.fThreshold;
  fTrigger                 = rhs.fTrigger;
  fStartTime               = rhs.fStartTime;
  fLength                  = rhs.fLength;
  fLargestPMTFraction      = rhs.fLargestPMTFraction;
  fLargeSecLargePMTFraction      = rhs.fLargeSecLargePMTFraction;
  fChargeSigma             = rhs.fChargeSigma;
  fTimeSigma               = rhs.fTimeSigma;
  fCoatedStartTime         = rhs.fCoatedStartTime;
  fUncoatedStartTime       = rhs.fUncoatedStartTime;
  fCoatedFittedStartTime   = rhs.fCoatedFittedStartTime;
  fUncoatedFittedStartTime = rhs.fUncoatedFittedStartTime;
  fVetoBottomStartTime     = rhs.fVetoBottomStartTime;
  fVetoTopStartTime        = rhs.fVetoTopStartTime;
  fVetoBackStartTime       = rhs.fVetoBackStartTime;
  fVetoRightStartTime      = rhs.fVetoRightStartTime;
  fVetoLeftStartTime       = rhs.fVetoLeftStartTime;
  fVetoFrontStartTime      = rhs.fVetoFrontStartTime;
  fCoatedAmp               = rhs.fCoatedAmp;
  fUncoatedAmp             = rhs.fUncoatedAmp;
  fVetoTopAmp              = rhs.fVetoTopAmp;
  fVetoBottomAmp           = rhs.fVetoBottomAmp;
  fVetoLeftAmp             = rhs.fVetoLeftAmp;
  fVetoRightAmp            = rhs.fVetoRightAmp;
  fVetoBackAmp             = rhs.fVetoBackAmp;
  fVetoFrontAmp            = rhs.fVetoFrontAmp;
  fHitXPosition            = rhs.fHitXPosition;
  fHitYPosition            = rhs.fHitYPosition;
  fHitZPosition            = rhs.fHitZPosition;
  fHitXPosTime             = rhs.fHitXPosTime;
  fHitYPosTime             = rhs.fHitYPosTime;
  fHitZPosTime             = rhs.fHitZPosTime;
  fHitTPosTime             = rhs.fHitTPosTime;
  for (int i = 0; i <2; ++i) {
    fCoatedInt[i]          = rhs.fCoatedInt[i];
    fUncoatedInt[i]        = rhs.fUncoatedInt[i];
    fCoatedFittedInt[i]    = rhs.fCoatedFittedInt[i];
    fUncoatedFittedInt[i]  = rhs.fUncoatedFittedInt[i];
    fVetoTopInt[i]         = rhs.fVetoTopInt[i];
    fVetoBottomInt[i]      = rhs.fVetoBottomInt[i];
    fVetoLeftInt[i]        = rhs.fVetoLeftInt[i];
    fVetoRightInt[i]       = rhs.fVetoRightInt[i];
    fVetoBackInt[i]        = rhs.fVetoBackInt[i];
    fVetoFrontInt[i]       = rhs.fVetoFrontInt[i];

    fNCoated[i]            = rhs.fNCoated[i];
    fNUncoated[i]          = rhs.fNUncoated[i];
    fNVetoBottom[i]        = rhs.fNVetoBottom[i];
    fNVetoTop[i]           = rhs.fNVetoTop[i];
    fNVetoBack[i]          = rhs.fNVetoBack[i];
    fNVetoLeft[i]          = rhs.fNVetoLeft[i];
    fNVetoRight[i]         = rhs.fNVetoRight[i];
    fNVetoFront[i]         = rhs.fNVetoFront[i];
  }

  fOtherEvents             = rhs.fOtherEvents;
  fOtherEventsTime         = rhs.fOtherEventsTime;

  fPulses                  = rhs.fPulses;
  fTimeToPositionCoated    = rhs.fTimeToPositionCoated;
  fTimeToPositionUncoated  = rhs.fTimeToPositionUncoated;

  fPassedThresholdCount    = rhs.fPassedThresholdCount;
  fPassedThresholdInt      = rhs.fPassedThresholdInt;
  fPassedThresholdIntDer   = rhs.fPassedThresholdIntDer;
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
  fComputerSecIntoEpoch = 0;
  fComputerNSIntoSec = 0;

  fThreshold          = 0.f;
  fTrigger            = 0;
  fLargestPMTFraction = 1.f;
  fLargeSecLargePMTFraction = 1.f;
  fChargeSigma = 1.f;
  fTimeSigma = 1.f;


  // Time variables
  // If the variable is -1 then the corresponding multiplicity count should be equal to 0
  fStartTime               = -1;
  fLength                  = -1;
  fCoatedStartTime         = -1;
  fUncoatedStartTime       = -1;
  fCoatedFittedStartTime   = -1;
  fUncoatedFittedStartTime = -1;
  fVetoBottomStartTime     = -1;
  fVetoTopStartTime        = -1;
  fVetoBackStartTime       = -1;
  fVetoRightStartTime      = -1;
  fVetoLeftStartTime       = -1;
  fVetoFrontStartTime      = -1;

  // Amplitude variables
  // The amplitude is calculated for the event only if the PMT was above threshold
  fCoatedAmp     = 0.f;
  fUncoatedAmp   = 0.f;
  fVetoTopAmp    = 0.f;
  fVetoBottomAmp = 0.f;
  fVetoLeftAmp   = 0.f;
  fVetoRightAmp  = 0.f;
  fVetoBackAmp   = 0.f;
  fVetoFrontAmp  = 0.f;

  // Integral variables
  // The integral is calculated for the event only if the PMT was above threshold
  // independent of the value of the integral or when the PMT crossed threshold
  fCoatedInt[0]         = 0.f;
  fUncoatedInt[0]       = 0.f;
  fCoatedFittedInt[0]   = 0.f;
  fUncoatedFittedInt[0] = 0.f;
  fVetoTopInt[0]        = 0.f;
  fVetoBottomInt[0]     = 0.f;
  fVetoLeftInt[0]       = 0.f;
  fVetoRightInt[0]      = 0.f;
  fVetoBackInt[0]       = 0.f;
  fVetoFrontInt[0]      = 0.f;

  fCoatedInt[1]         = 0.f;
  fUncoatedInt[1]       = 0.f;
  fCoatedFittedInt[1]   = 0.f;
  fUncoatedFittedInt[1] = 0.f;
  fVetoTopInt[1]        = 0.f;
  fVetoBottomInt[1]     = 0.f;
  fVetoLeftInt[1]       = 0.f;
  fVetoRightInt[1]      = 0.f;
  fVetoBackInt[1]       = 0.f;
  fVetoFrontInt[1]      = 0.f;

  // Position variables
  fHitXPosition = 0.f;
  fHitYPosition = 0.f;
  fHitZPosition = 0.f;
  fHitXPosTime  = 0.f;
  fHitYPosTime  = 0.f;
  fHitZPosTime  = 0.f;
  fHitTPosTime  = 0.f;

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

  fOtherEvents = 0;
  fOtherEventsTime.clear();

  if (!fTimeToPositionCoated.empty()) {
    fTimeToPositionCoated.clear();
  }
  if (!fTimeToPositionUncoated.empty()) {
    fTimeToPositionUncoated.clear();
  }
  if (!fPulses.empty()) {
    fPulses.clear();
  }

  if (!fPassedThresholdCount.empty()) {
    fPassedThresholdCount.clear();
  }
  if (!fPassedThresholdInt.empty()) {
    fPassedThresholdInt.clear();
  }

  for (int i=0; i < 161; ++i) {
    fPassedThresholdIntDer[i].fill(0);
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
float SimplifiedEvent::GetAmplitudeTank() 
{
  return fCoatedAmp + fUncoatedAmp;
}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetAmplitudeVeto() 
{
  return fVetoBottomAmp + fVetoTopAmp + fVetoBackAmp + fVetoLeftAmp + fVetoFrontAmp + fVetoRightAmp;
}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetIntegralTank(bool prompt) 
{
  return fCoatedInt[prompt] + fUncoatedInt[prompt];
}

//-------------------------------------------------------------------------------------------------
float SimplifiedEvent::GetIntegralVeto(bool prompt) 
{
  return fVetoBottomInt[prompt] + fVetoTopInt[prompt] + fVetoBackInt[prompt] + fVetoLeftInt[prompt] + fVetoFrontInt[prompt] + fVetoRightInt[prompt];
}

/*!**********************************************
 * \fn void SimplifiedEvent::AddWaveforms(const std::vector<float> & count, const std::vector<float> & integral)
 * \brief Set the accumulated waveforms to what was passed
 * \param[in] count The PMT multiplicity accumulated waveform
 * \param[in] integral The integral charge accumulated waveform
 ***********************************************/
void SimplifiedEvent::AddWaveforms(const std::vector<float> & count, const std::vector<float> & integral)
{
  fPassedThresholdCount = count;
  fPassedThresholdInt = integral;
}
