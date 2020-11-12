/*!**********************************************
 * \file SinglePulse.h
 * \brief Header file for the #SinglePulse class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef SinglePulse_h
#define SinglePulse_h

#include <vector>
#include <iostream>
#include <sys/types.h>
#include <vector>

#include "TObject.h"

/*!**********************************************
 * \class SinglePulse
 * \brief The information for a given pulse found
 ***********************************************/
class SinglePulse : public TObject
{
  public:
    SinglePulse(int key = 0);
    SinglePulse(const SinglePulse & p);
    ~SinglePulse();

    void Reset();

    void SetKey(int key) { fKey= key; }

    void SetADCToPE (float value) { fADCToPE = value; }
    void SetTriggerOffset (float trigOffset) { fTrigOffset = trigOffset; }
    void SetPMTOffset ( float pmtOffset) { fPMTOffset = pmtOffset; }

    void SetLength(float value) { fLength = value; }
    void SetMaxDerValue(float value) { fMaxDerValue= value; }
    void SetTime(float value) { fTime = value; }

    void SetBaseline(float value) { fBaseline = value; }
    void SetAmplitude(float value) { fAmplitude = value; }
    void SetIntegral(float value) { fIntegral = value; }

    size_t GetKey() { return fKey; }

    float GetADCToPE() { return fADCToPE; }
    float GetTriggerOffset() { return fTrigOffset; }
    float GetPMTOffset() { return fPMTOffset; }

    float GetLength() { return fLength; }
    float GetMaxDerValue() { return fMaxDerValue; }
    float GetTime() { return fTime; }

    float GetBaseline() { return fBaseline; }
    float GetAmplitude() { return fAmplitude; }
    float GetIntegral() { return fIntegral; }

    size_t GetKey() const { return fKey; }

    float GetADCToPE() const { return fADCToPE; }
    float GetTriggerOffset() const { return fTrigOffset; }
    float GetPMTOffset() const { return fPMTOffset; }

    float GetLength() const { return fLength; }
    float GetMaxDerValue() const { return fMaxDerValue; }
    float GetTime() const { return fTime; }

    float GetBaseline() const { return fBaseline; }
    float GetAmplitude() const { return fAmplitude; }
    float GetIntegral() const { return fIntegral; }

    //float GetSample(size_t index) { return fWaveform.at(index); }
    //float GetSample(size_t index) const { return fWaveform.at(index); }
    //void SetSample(size_t index, float value) { fWaveform.at(index) = value; }

    void SetWaveformStart(size_t value) { fWaveformStart = value; }
    void SetWaveformEnd(size_t value) { fWaveformEnd = value; }

    size_t GetWaveformStart() { return fWaveformStart; }
    size_t GetWaveformEnd() { return fWaveformEnd; }
    size_t GetWaveformStart() const { return fWaveformStart; }
    size_t GetWaveformEnd() const { return fWaveformEnd; }

    const std::vector<float> & GetWaveform() { return fWaveform; }
    const std::vector<float> & GetWaveform() const { return fWaveform; }
    void AddSample(float sample) { fWaveform.push_back(sample); }
    void SetWaveform(const std::vector<float> & waveform) { fWaveform = waveform; }

    void Append(const SinglePulse & rhs);

    SinglePulse & operator=(const SinglePulse & rhs);

  private:
    size_t fKey;

    float fADCToPE; ///< the ADC to PE conversion value used to determine fIntegral value in PEs
    float fTrigOffset; ///< trigger offset for the board
    float fPMTOffset; ///< offset because of pmt type

    float fAmplitude; ///< the maximum amplitude (baseline - minADC) for the pulse
    float fBaseline; ///< the calculate pulse baseline
    float fMaxDerValue; ///< the absolute maximum derivative value in the pulse
    float fIntegral; ///< the integral if the waveform from [time, time + length) below the baseline
    float fTime; ///< the start time of the pulse in sample number (offsets already applied)
    float fLength; ///< the length of the pulse in sample number

    size_t fWaveformStart; ///< the sample number corresponding to the start of the saved waveform (offsets not applied)
    size_t fWaveformEnd; ///< the sample number corresponding to the end of the saved waveform (offsets not applied)

    /// The waveform that made up the pulse (currently not being saved)
    std::vector<float> fWaveform; //!

  ClassDef(SinglePulse,3)

};

#endif // #ifndef SinglePulse_h

