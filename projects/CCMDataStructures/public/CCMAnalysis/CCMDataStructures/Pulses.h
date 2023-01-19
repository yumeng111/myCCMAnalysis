/*!**********************************************
 * \file Pulses.h
 * \brief Header file for the #Pulses class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef Pulses_h
#define Pulses_h

#include <array>
#include <vector>
#include <iostream>
#include <sys/types.h>
#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
#include <icetray/serialization.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3DefaultName.h>
#endif // __MAKECINT__ || __ROOTCLING__

#include "CCMAnalysis/CCMDataStructures/SinglePulse.h"
#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMUtils/PMTInfoMap.h"

static const unsigned legacy_pulses_version_ = 1;

#include "TObject.h"

/*!**********************************************
 * \class Pulses
 * \brief Container of the #SinglePulse found in a given DAQ window.
 *
 * Container for a vector of #SinglePulse pulses found in the detector
 * and the methods to find such pulses. The class is saved in the
 * data file.
 ***********************************************/
class Pulses : public TObject
#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
    , public I3FrameObject
#endif // __MAKECINT__ || __ROOTCLING__
{
  public:
    Pulses(int board = 0, int channel = 0);
    Pulses(const Pulses & p);
    ~Pulses();

    void Reset();
    void ClearPulses();

    void RemovePulsesByThreshold();

    /// \fn void SetEventNumber(unsigned int value)
    /// \brief Set the event (trigger) number of the DAQ window currently looking at
    /// \param[in] value The event number
    void SetEventNumber(unsigned int value) { fEventNumber = value; }
    /// \fn void SetComputerSecIntoEpoch(unsigned int value)
    /// \brief Set the computer time of the DAQ window in s
    /// \param[in] value The computer time in s
    void SetComputerSecIntoEpoch(unsigned int value) { fComputerSecIntoEpoch = value; }
    /// \fn void SetComputerNSIntoSec(unsigned int value)
    /// \brief Set the computer time of the DAQ window in ns
    /// \param[in] value The computer time in ns
    void SetComputerNSIntoSec(unsigned int value) { fComputerNSIntoSec = value; }

    /// \fn unsigned int GetEventNumber() const
    /// \return #fEventNumber
    unsigned int GetEventNumber() const { return fEventNumber; }
    /// \fn unsigned int GetComputerSecIntoEpoch() const
    /// \return #fComputerSecIntoEpoch
    unsigned int GetComputerSecIntoEpoch() const { return fComputerSecIntoEpoch; }
    /// \fn unsigned int GetComputerNSIntoSec() const
    /// \return #fComputerNSIntoSec
    unsigned int GetComputerNSIntoSec() const { return fComputerNSIntoSec; }

    /// \fn unsigned int GetEventNumber()
    /// \return #fEventNumber
    unsigned int GetEventNumber() { return fEventNumber; }
    /// \fn unsigned int GetComputerSecIntoEpoch()
    /// \return #fComputerSecIntoEpoch
    unsigned int GetComputerSecIntoEpoch() { return fComputerSecIntoEpoch; }
    /// \fn unsigned int GetComputerNSIntoSec()
    /// \return #fComputerNSIntoSec
    unsigned int GetComputerNSIntoSec() { return fComputerNSIntoSec; }

    /// \fn void SetBoardChannel()
    /// \brief Set the board and channel number of the pulse currently being looked at
    /// \param[in] board Digitizer board number
    /// \param[in] channel Channel number on the digitizer board
    void SetBoardChannel(int board, int channel) { fBoard = board; fChannel = channel; }
    /// \fn void SetKey()
    /// \brief Set the key for the pmt currently being looked at
    /// \param[in] key The key value
    void SetKey(int key) { PMTInfoMap::DecodeKey(key,fBoard,fChannel); }

    void DerivativeFilter(const u_int16_t input[], int length, float triggerTime, 
        int points = 5, float threshold = 0, float trigOffset = 0, float pmtOffset = 0, float adcToPE = 1);
    void MCFilter(const std::vector<double> & input, float triggerTime);
    void SmoothWaveform(const std::vector<u_int16_t> & input);
    void Sort();

    double GetTimeWindowStart() { return fgkTimeWindowStart; }
    double GetSampleTimeWidth() { return fgkSampleTimeWidth; }
    double GetTimeWindowStart() const { return fgkTimeWindowStart; }
    double GetSampleTimeWidth() const { return fgkSampleTimeWidth; }

    size_t GetKey(size_t pos);

    size_t GetNumPulses();
    float GetTriggerTime();

    float GetPulseLength(size_t pos);
    float GetPulseMaxDerValue(size_t pos);
    float GetPulseTime(size_t pos);

    float GetPulseBaseline(size_t pos);
    float GetPulseAmplitude(size_t pos);
    float GetPulseIntegral(size_t pos);

    size_t GetKey(size_t pos) const;

    size_t GetNumPulses() const;
    float GetTriggerTime() const;

    float GetPulseLength(size_t pos) const;
    float GetPulseMaxDerValue(size_t pos) const;
    float GetPulseTime(size_t pos) const;

    float GetPulseBaseline(size_t pos) const;
    float GetPulseAmplitude(size_t pos) const;
    float GetPulseIntegral(size_t pos) const;

    void AddSinglePulse(const SinglePulse & pulse);
    void RemoveSinglePulse(size_t pos);
    void RemoveSinglePulse(std::vector<SinglePulse>::iterator it);
    const SinglePulse * GetSinglePulse(size_t pos);
    const SinglePulse * GetSinglePulse(size_t pos) const;

    float LongestPulse(float time = Utility::fgkNumBins) const;
    float LargestPulse(float time = Utility::fgkNumBins) const;

    std::array<int,Utility::fgkNumBins> GetOrigWaveformAt(int i) { return fOrigArray.at(i); }
    std::array<float,Utility::fgkNumBins> GetRawWaveformAt(int i) { return fRawArray.at(i); }
    std::array<float,Utility::fgkNumBins> GetWaveformAt(int i) { return fSmoothArray.at(i); }
    std::array<float,Utility::fgkNumBins> GetWaveformDerAt(int i) { return fSmoothArrayDer.at(i); }
    const DAQWF2DArrayI & GetOrigWaveform() { return fOrigArray; }
    const DAQWF2DArrayF & GetRawWaveform() { return fRawArray; }
    const DAQWF2DArrayF & GetWaveform() { return fSmoothArray; }
    const DAQWF2DArrayF & GetWaveformDer() { return fSmoothArrayDer; }

    size_t FindFirstAfter(float time, int startLoc = 0) const;
    const SinglePulse * FindPreviousPulse(int board, int channel);
    const SinglePulse * FindPreviousPulse(unsigned int key);

    Pulses & operator=(const Pulses & rhs);
    void CopyPulses(const Pulses & rhs, int startBoard, int endBoard=10, int startChannel=0, int endChannel=15);

    void ShiftTimeOffset(const double & timeOffset);

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive & ar, unsigned version);
#endif // __MAKECINT__ || __ROOTCLING__

  private:
    static bool SortCondition(const SinglePulse & a, const SinglePulse & b);

  private:
    /// The board currently being looked at (not saved)
    size_t fBoard; //!
    /// The channel currently being looked at (not saved)
    size_t fChannel; //!

    /// The time of the window with sample number == 0 (not saved)
    /// in microseconds
    constexpr static const double fgkTimeWindowStart = -9.92; //!
    /// The sample time width (not saved)
    /// in microseconds
    constexpr static const double fgkSampleTimeWidth = 2e-3; //!

    /// The event number
    unsigned int fEventNumber;
    /// The time of the window in s
    unsigned int fComputerSecIntoEpoch;
    /// The time of the window in ns
    unsigned int fComputerNSIntoSec;

    /// The number of pulses found
    size_t fNumPulses;

    /// The trigger time of the trigger
    float fTriggerTime;

    /// The vector of pulses found
    std::vector<SinglePulse> fPulses;

    /// Array of the waveform (not saved)
    DAQWF2DArrayI fOrigArray; //!
    /// Array of the waveform after exponential moving average (not saved)
    DAQWF2DArrayF fRawArray; //!
    /// Array of the waveform after 3 point average (not saved)
    DAQWF2DArrayF fSmoothArray; //!
    /// Array of the derivative of the waveform after 3 point aveage (not saved)
    DAQWF2DArrayF fSmoothArrayDer; //!

    /// Pointer to the current #SinglePulse (not saved)
    SinglePulse * fCurrPulse; //!



  ClassDef(Pulses,legacy_pulses_version_)

};

#if !(defined(__MAKECINT__) || defined(__ROOTCLING__))
I3_DEFAULT_NAME(Pulses, LegacyPulses);
I3_POINTER_TYPEDEFS(Pulses);
I3_CLASS_VERSION(Pulses, legacy_pulses_version_);
#endif // __MAKECINT__ || __ROOTCLING__

#endif // #ifndef Pulses_h
