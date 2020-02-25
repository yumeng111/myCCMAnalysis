/*!**********************************************
 * \file SimplifiedEvent.h
 * \brief Header file for the #SimplifiedEvent class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef SimplifiedEvent_h
#define SimplifiedEvent_h

#include "SinglePulse.h"
#include "TObject.h"

#include <iostream>
#include <vector>
#include <array>

//class BoardInfo;
//class ChannelInfo;

/*!**********************************************
 * \class SimplifiedEvent
 * \brief All the information needed after finding an event
 *
 * Not all parameters are filled with every event definition
 * but this is a good list of all the various parameters
 * one might want to calculate
 ***********************************************/
class SimplifiedEvent : public TObject
{
  public :
    SimplifiedEvent();
    SimplifiedEvent(const SimplifiedEvent & rhs);
    ~SimplifiedEvent();

    void Reset();

    void SetComputerSecIntoEpoch(unsigned int value) { fComputerSecIntoEpoch = value; }
    void SetComputerNSIntoSec(unsigned int value) { fComputerNSIntoSec = value; }

    unsigned int GetComputerSecIntoEpoch() { return fComputerSecIntoEpoch; }
    unsigned int GetComputerNSIntoSec() { return fComputerNSIntoSec; }

    ////////////////////////////////////
    // Set Functions
    ////////////////////////////////////
    void SetThreshold          ( float value  ) { fThreshold = value; }
    void SetTriggerNumber      ( int value    ) { fTrigger = value; }
    void SetLargestPMTFraction ( double value ) { fLargestPMTFraction = value; }
    void SetLargeSecLargePMTFraction ( double value ) { fLargeSecLargePMTFraction = value; }
    void SetChargeSigma ( double value ) { fChargeSigma = value; }
    void SetTimeSigma ( double value ) { fTimeSigma = value; }

    void SetStartTime               (double value) { fStartTime = value; }
    void SetLength                  (double value) { fLength = value; }
    void SetStartTimeCoated         (double value) { fCoatedStartTime = value; }
    void SetStartTimeUncoated       (double value) { fUncoatedStartTime = value; }
    void SetFittedStartTimeCoated   (double value) { fCoatedFittedStartTime = value; }
    void SetFittedStartTimeUncoated (double value) { fUncoatedFittedStartTime = value; }
    void SetStartTimeVetoBottom     (double value) { fVetoBottomStartTime = value; }
    void SetStartTimeVetoTop        (double value) { fVetoTopStartTime = value; }
    void SetStartTimeVetoBack       (double value) { fVetoBackStartTime = value; }
    void SetStartTimeVetoRight      (double value) { fVetoRightStartTime = value; }
    void SetStartTimeVetoLeft       (double value) { fVetoLeftStartTime = value; }
    void SetStartTimeVetoFront      (double value) { fVetoFrontStartTime = value; }

    // Amplitude variables
    void SetAmplitudeCoated     (float value) { fCoatedAmp = value; }
    void SetAmplitudeUncoated   (float value) { fUncoatedAmp = value; }
    void SetAmplitudeVetoTop    (float value) { fVetoTopAmp = value; }
    void SetAmplitudeVetoBottom (float value) { fVetoBottomAmp = value; }
    void SetAmplitudeVetoLeft   (float value) { fVetoLeftAmp = value; }
    void SetAmplitudeVetoRight  (float value) { fVetoRightAmp = value; }
    void SetAmplitudeVetoBack   (float value) { fVetoBackAmp = value; }
    void SetAmplitudeVetoFront  (float value) { fVetoFrontAmp = value; }

    // Integral variables
    void SetIntegralCoated         (float value, bool prompt = true) { fCoatedInt[prompt] = value; }
    void SetIntegralUncoated       (float value, bool prompt = true) { fUncoatedInt[prompt] = value; }
    void SetFittedIntegralCoated   (float value, bool prompt = true) { fCoatedFittedInt[prompt] = value; }
    void SetFittedIntegralUncoated (float value, bool prompt = true) { fUncoatedFittedInt[prompt] = value; }
    void SetIntegralVetoTop        (float value, bool prompt = true) { fVetoTopInt[prompt] = value; }
    void SetIntegralVetoBottom     (float value, bool prompt = true) { fVetoBottomInt[prompt] = value; }
    void SetIntegralVetoLeft       (float value, bool prompt = true) { fVetoLeftInt[prompt] = value; }
    void SetIntegralVetoRight      (float value, bool prompt = true) { fVetoRightInt[prompt] = value; }
    void SetIntegralVetoBack       (float value, bool prompt = true) { fVetoBackInt[prompt] = value; }
    void SetIntegralVetoFront      (float value, bool prompt = true) { fVetoFrontInt[prompt] = value; }

    void SetXPosition (float value) { fHitXPosition = value; }
    void SetYPosition (float value) { fHitYPosition = value; }
    void SetZPosition (float value) { fHitZPosition = value; }
    void SetXPosTime  (float value) { fHitXPosTime = value; }
    void SetYPosTime  (float value) { fHitYPosTime = value; }
    void SetZPosTime  (float value) { fHitZPosTime = value; }
    void SetTPosTime  (float value) { fHitTPosTime = value; }

    // Multiplicity variables
    void SetNumCoated     (float value, bool prompt = true) { fNCoated[prompt] = value; }
    void SetNumUncoated   (float value, bool prompt = true) { fNUncoated[prompt] = value; }
    void SetNumVetoBottom (float value, bool prompt = true) { fNVetoBottom[prompt] = value; }
    void SetNumVetoTop    (float value, bool prompt = true) { fNVetoTop[prompt] = value; }
    void SetNumVetoBack   (float value, bool prompt = true) { fNVetoBack[prompt] = value; }
    void SetNumVetoLeft   (float value, bool prompt = true) { fNVetoLeft[prompt] = value; }
    void SetNumVetoRight  (float value, bool prompt = true) { fNVetoRight[prompt] = value; }
    void SetNumVetoFront  (float value, bool prompt = true) { fNVetoFront[prompt] = value; }

    ////////////////////////////////////
    // Add Functions
    ////////////////////////////////////
    void AddPulseToPMT             (const SinglePulse & pulse)     { fPulses.push_back(pulse); }
    void AddTimeToPositionCoated   (float value) { fTimeToPositionCoated.push_back(value); }
    void AddTimeToPositionUncoated (float value) { fTimeToPositionUncoated.push_back(value); }
    void AddWaveforms(const std::vector<float> & count, const std::vector<float> & integral);

    const std::vector<float> & GetWaveformCount() { return fPassedThresholdCount; }
    const std::vector<float> & GetWaveformInt() { return fPassedThresholdInt; }

    // Amplitude variables
    void AddAmplitudeCoated     (float value) { fCoatedAmp += value; }
    void AddAmplitudeUncoated   (float value) { fUncoatedAmp += value; }
    void AddAmplitudeVetoTop    (float value) { fVetoTopAmp += value; }
    void AddAmplitudeVetoBottom (float value) { fVetoBottomAmp += value; }
    void AddAmplitudeVetoLeft   (float value) { fVetoLeftAmp += value; }
    void AddAmplitudeVetoRight  (float value) { fVetoRightAmp += value; }
    void AddAmplitudeVetoBack   (float value) { fVetoBackAmp += value; }
    void AddAmplitudeVetoFront  (float value) { fVetoFrontAmp += value; }

    // Integral variables
    void AddIntegralCoated           (float value, bool prompt = true) { fCoatedInt[prompt] += value; }
    void AddIntegralUncoated         (float value, bool prompt = true) { fUncoatedInt[prompt] += value; }
    void AddFittedIntegralCoated     (float value, bool prompt = true) { fCoatedFittedInt[prompt] += value; }
    void AddFittedIntegralUncoated   (float value, bool prompt = true) { fUncoatedFittedInt[prompt] += value; }
    void AddIntegralVetoTop          (float value, bool prompt = true) { fVetoTopInt[prompt] += value; }
    void AddIntegralVetoBottom       (float value, bool prompt = true) { fVetoBottomInt[prompt] += value; }
    void AddIntegralVetoLeft         (float value, bool prompt = true) { fVetoLeftInt[prompt] += value; }
    void AddIntegralVetoRight        (float value, bool prompt = true) { fVetoRightInt[prompt] += value; }
    void AddIntegralVetoBack         (float value, bool prompt = true) { fVetoBackInt[prompt] += value; }
    void AddIntegralVetoFront        (float value, bool prompt = true) { fVetoFrontInt[prompt] += value; }

    // Multiplicity variables
    void AddNumCoated     (int value = 1, bool prompt = true) { fNCoated[prompt] += value; }
    void AddNumUncoated   (int value = 1, bool prompt = true) { fNUncoated[prompt] += value; }
    void AddNumVetoBottom (int value = 1, bool prompt = true) { fNVetoBottom[prompt] += value; }
    void AddNumVetoTop    (int value = 1, bool prompt = true) { fNVetoTop[prompt] += value; }
    void AddNumVetoBack   (int value = 1, bool prompt = true) { fNVetoBack[prompt] += value; }
    void AddNumVetoLeft   (int value = 1, bool prompt = true) { fNVetoLeft[prompt] += value; }
    void AddNumVetoRight  (int value = 1, bool prompt = true) { fNVetoRight[prompt] += value; }
    void AddNumVetoFront  (int value = 1, bool prompt = true) { fNVetoFront[prompt] += value; }

    ////////////////////////////////////
    // Get Functions
    ////////////////////////////////////
    float GetThreshold           ( ) { return fThreshold; }
    int   GetTriggerNumber       ( ) { return fTrigger; }
    double GetLargestPMTFraction ( ) { return fLargestPMTFraction; }
    double GetLargeSecLargePMTFraction ( ) { return fLargeSecLargePMTFraction; }
    double GetChargeSigma ( ) { return fChargeSigma; }
    double GetTimeSigma ( ) { return fTimeSigma; }


    // Time variables
    double GetStartTime               () { return fStartTime; }
    double GetLength                  () { return fLength; }
    double GetStartTimeCoated         () { return fCoatedStartTime; }
    double GetStartTimeUncoated       () { return fUncoatedStartTime; }
    double GetFittedStartTimeCoated   () { return fCoatedFittedStartTime; }
    double GetFittedStartTimeUncoated () { return fUncoatedFittedStartTime; }
    double GetStartTimeVetoBottom     () { return fVetoBottomStartTime; }
    double GetStartTimeVetoTop        () { return fVetoTopStartTime; }
    double GetStartTimeVetoBack       () { return fVetoBackStartTime; }
    double GetStartTimeVetoRight      () { return fVetoRightStartTime; }
    double GetStartTimeVetoLeft       () { return fVetoLeftStartTime; }
    double GetStartTimeVetoFront      () { return fVetoFrontStartTime; }
    double GetStartTimeTank();
    double GetStartTimeVeto();

    size_t GetNumberOfCoatedPulses   ()           { return fTimeToPositionCoated.size(); }
    size_t GetNumberOfUncoatedPulses ()           { return fTimeToPositionUncoated.size(); }
    float GetTimeToPositionCoated    (int loc)    { return fTimeToPositionCoated.at(loc); }
    float GetTimeToPositionUncoated  (int loc)    { return fTimeToPositionUncoated.at(loc); }

    size_t  GetNumPulse() { return fPulses.size(); }
    const SinglePulse & GetPulse(int loc) { return fPulses.at(loc); }

    // Amplitude variables
    float GetAmplitudeCoated     () { return fCoatedAmp; }
    float GetAmplitudeUncoated   () { return fUncoatedAmp; }
    float GetAmplitudeVetoTop    () { return fVetoTopAmp; }
    float GetAmplitudeVetoBottom () { return fVetoBottomAmp; }
    float GetAmplitudeVetoLeft   () { return fVetoLeftAmp; }
    float GetAmplitudeVetoRight  () { return fVetoRightAmp; }
    float GetAmplitudeVetoBack   () { return fVetoBackAmp; }
    float GetAmplitudeVetoFront  () { return fVetoFrontAmp; }
    float GetAmplitudeTank();
    float GetAmplitudeVeto();

    // Integral variables
    float GetIntegralCoated         (bool prompt = true) { return fCoatedInt[prompt]; }
    float GetIntegralUncoated       (bool prompt = true) { return fUncoatedInt[prompt]; }
    float GetFittedIntegralCoated   (bool prompt = true) { return fCoatedFittedInt[prompt]; }
    float GetFittedIntegralUncoated (bool prompt = true) { return fUncoatedFittedInt[prompt]; }
    float GetIntegralVetoTop        (bool prompt = true) { return fVetoTopInt[prompt]; }
    float GetIntegralVetoBottom     (bool prompt = true) { return fVetoBottomInt[prompt]; }
    float GetIntegralVetoLeft       (bool prompt = true) { return fVetoLeftInt[prompt]; }
    float GetIntegralVetoRight      (bool prompt = true) { return fVetoRightInt[prompt]; }
    float GetIntegralVetoBack       (bool prompt = true) { return fVetoBackInt[prompt]; }
    float GetIntegralVetoFront      (bool prompt = true) { return fVetoFrontInt[prompt]; }
    float GetIntegralTank(bool prompt = true);
    float GetIntegralVeto(bool prompt = true);

    float GetXPosition () { return fHitXPosition; }
    float GetYPosition () { return fHitYPosition; }
    float GetZPosition () { return fHitZPosition; }
    float GetXPosTime  () { return fHitXPosTime; }
    float GetYPosTime  () { return fHitYPosTime; }
    float GetZPosTime  () { return fHitZPosTime; }
    float GetTPosTime  () { return fHitTPosTime; }

    // Multiplicity variables
    float    GetNumCoated     (bool prompt = true) { return fNCoated[prompt]; }
    float    GetNumUncoated   (bool prompt = true) { return fNUncoated[prompt]; }
    float    GetNumVetoBottom (bool prompt = true) { return fNVetoBottom[prompt]; }
    float    GetNumVetoTop    (bool prompt = true) { return fNVetoTop[prompt]; }
    float    GetNumVetoBack   (bool prompt = true) { return fNVetoBack[prompt]; }
    float    GetNumVetoLeft   (bool prompt = true) { return fNVetoLeft[prompt]; }
    float    GetNumVetoRight  (bool prompt = true) { return fNVetoRight[prompt]; }
    float    GetNumVetoFront  (bool prompt = true) { return fNVetoFront[prompt]; }
    float    GetNumTank(bool prompt = true);
    float    GetNumVeto(bool prompt = true);
    float    GetNumVetoSide(bool prompt = true);

    void  SetOtherEvents (int value, std::vector<float> & time) { fOtherEvents = value; fOtherEventsTime = time; }
    void  AddOtherEvents (int value = 1, float time = 0) { fOtherEvents += value; fOtherEventsTime.push_back(time); }
    int   GetOtherEvents () { return fOtherEvents; }
    const std::vector<float> & GetOtherEventsTime () { return fOtherEventsTime; }

  protected:
    static float fgkStartDAQWindow; //!
    static float fgkSampleWidth; //!

    float fThreshold;                    //!
    int   fTrigger;                      ///< trigger number
    int fOtherEvents;                    ///< true if there is another triggerable event
    std::vector<float> fOtherEventsTime; ///< time of all the other events
    double fLargestPMTFraction;          ///< fraction of PMT with the most charge
    double fLargeSecLargePMTFraction;    ///< fraction of PMT with the second most charge
    double fChargeSigma;                 ///< sigma of charge distribution
    double fTimeSigma;                   ///< sigma of time distribution


    // Time variables
    // If the variable is -1 then the corresponding multiplicity count should be equal to 0
    double fStartTime;               ///< time of event in bin sample
    double fLength;                  ///< length of event
    double fCoatedStartTime;         ///< time of first coated tank event
    double fUncoatedStartTime;       ///< time of first uncoated event
    double fCoatedFittedStartTime;   ///< time of first coated tank event
    double fUncoatedFittedStartTime; ///< time of first uncoated event
    double fVetoBottomStartTime;     ///< time of first veto on the bottom
    double fVetoTopStartTime;        ///< time of first veto on the top
    double fVetoBackStartTime;       ///< time of first veto on the back
    double fVetoRightStartTime;      ///< time of first veto on the right
    double fVetoLeftStartTime;       ///< time of first veto on the left
    double fVetoFrontStartTime;      ///< time of first veto on the front

    // Amplitude variables
    // The amplitude is calculated for the event only if the PMT was above threshold
    float fCoatedAmp;     ///< of coated tubes
    float fUncoatedAmp;   ///< of uncoated tubes
    float fVetoTopAmp;    ///< of top veto tubes
    float fVetoBottomAmp; ///< of bottom veto tubes
    float fVetoLeftAmp;   ///< of left veto tubes
    float fVetoRightAmp;  ///< of right veto tubes
    float fVetoBackAmp;   ///< of back veto tubes
    float fVetoFrontAmp;  ///< of front veto tubes

    // Integral variables
    // The integral is calculated for the event only if the PMT was above threshold
    // independent of the value of the integral or when the PMT crossed threshold
    float fCoatedInt[2];         ///< of coated tubes
    float fUncoatedInt[2];       ///< of uncoated tubes
    float fCoatedFittedInt[2];   ///< of coated tubes
    float fUncoatedFittedInt[2]; ///< of uncoated tubes
    float fVetoTopInt[2];        ///< of top veto tubes
    float fVetoBottomInt[2];     ///< of bottom veto tubes
    float fVetoLeftInt[2];       ///< of left veto tubes
    float fVetoRightInt[2];      ///< of right veto tubes
    float fVetoBackInt[2];       ///< of back veto tubes
    float fVetoFrontInt[2];      ///< of front veto tubes

    float fHitXPosition; ///< of tank tubes (coated + uncoated)
    float fHitYPosition; ///< of tank tubes (coated + uncoated)
    float fHitZPosition; ///< of tank tubes (coated + uncoated)

    float fHitXPosTime; ///< of tank tubes (coated + uncoated)
    float fHitYPosTime; ///< of tank tubes (coated + uncoated)
    float fHitZPosTime; ///< of tank tubes (coated + uncoated)
    float fHitTPosTime; ///< of tank tubes (coated + uncoated)

    // Multiplicity variables
    float fNCoated[2];     ///< for coated tank tubes
    float fNUncoated[2];   ///< for just uncoated tubes
    float fNVetoBottom[2]; ///< for 8" bottom veto tubes
    float fNVetoTop[2];    ///< for 1" top veto tubes
    float fNVetoBack[2];   ///< for 1" cylinder veto tubes (columes 22 through 3)
    float fNVetoLeft[2];   ///< for 1" cylinder veto tubes (columes 4 through 9)
    float fNVetoRight[2];  ///< for 1" cylinder veto tubes (columes 16 through 21)
    float fNVetoFront[2];  ///< for 1" cylinder veto tubes (columes 10 through 15)

    /// Vector of the single pulses that made up the event
    std::vector<SinglePulse> fPulses; //!
    std::vector<float> fTimeToPositionCoated; //!
    std::vector<float> fTimeToPositionUncoated; //!
    std::vector<float> fPassedThresholdCount; //
    std::vector<float> fPassedThresholdInt; //
    //static std::array<std::array<float,4460>,161> fPassedThresholdCountDer; //!
    std::array<std::array<float,4460>,161> fPassedThresholdIntDer; //!

    unsigned int fComputerSecIntoEpoch;
    unsigned int fComputerNSIntoSec;

  ClassDef(SimplifiedEvent,9)  //SimplifiedEvent class
};

#endif
