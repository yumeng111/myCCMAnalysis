/*!**********************************************
 * \file Pulses.h
 * \brief Header file for the #Pulses class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef Pulses_h
#define Pulses_h

#include <vector>
#include <array>
#include <iostream>
#include <sys/types.h>

#include "SinglePulse.h"
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

    void DerivativeFilter(const u_int16_t input[], int length, float triggerTime, 
        int points = 5, float threshold = 0, float trigOffset = 0, float pmtOffset = 0, float adcToPE = 1);
    void SmoothWaveform(const std::vector<u_int16_t> & input);
    void Sort();

    double GetTimeWindowStart() { return fgkTimeWindowStart; }
    double GetSampleTimeWidth() { return fgkSampleTimeWidth; }
    double GetTimeWindowStart() const { return fgkTimeWindowStart; }
    double GetSampleTimeWidth() const { return fgkSampleTimeWidth; }

    size_t GetKey(size_t pos) { return fPulses[pos].GetKey(); }

    size_t GetNumPulses() { return fNumPulses; }
    float GetTriggerTime() { return fTriggerTime; }

    float GetPulseLength(size_t pos) { return fPulses[pos].GetLength(); }
    float GetPulseMaxDerValue(size_t pos) { return fPulses[pos].GetMaxDerValue(); }
    float GetPulseTime(size_t pos) { return fPulses[pos].GetTime(); }

    float GetPulseBaseline(size_t pos) { return fPulses[pos].GetBaseline(); }
    float GetPulseAmplitude(size_t pos) { return fPulses[pos].GetAmplitude(); }
    float GetPulseIntegral(size_t pos) { return fPulses[pos].GetIntegral(); }

    size_t GetKey(size_t pos) const { return fPulses[pos].GetKey(); }

    size_t GetNumPulses() const { return fNumPulses; }
    float GetTriggerTime() const { return fTriggerTime; }

    float GetPulseLength(size_t pos) const { return fPulses[pos].GetLength(); }
    float GetPulseMaxDerValue(size_t pos) const { return fPulses[pos].GetMaxDerValue(); }
    float GetPulseTime(size_t pos) const { return fPulses[pos].GetTime(); }

    float GetPulseBaseline(size_t pos) const { return fPulses[pos].GetBaseline(); }
    float GetPulseAmplitude(size_t pos) const { return fPulses[pos].GetAmplitude(); }
    float GetPulseIntegral(size_t pos) const { return fPulses[pos].GetIntegral(); }

    void AddSinglePulse(const SinglePulse & pulse) { fPulses.push_back(pulse); ++fNumPulses;}
    void RemoveSinglePulse(size_t pos) { fPulses.erase(fPulses.begin()+pos); --fNumPulses;}
    void RemoveSinglePulse(std::vector<SinglePulse>::iterator it) { fPulses.erase(it); --fNumPulses;}
    const SinglePulse * GetSinglePulse(size_t pos) { return &fPulses[pos]; }
    const SinglePulse * GetSinglePulse(size_t pos) const { return &fPulses[pos]; }

    float LongestPulse(float time = 8000) const;
    float LargestPulse(float time = 8000) const;

    std::array<int,8000> GetOrigWaveformAt(int i) { return fOrigArray.at(i); }
    std::array<float,8000> GetRawWaveformAt(int i) { return fRawArray.at(i); }
    std::array<float,8000> GetWaveformAt(int i) { return fSmoothArray.at(i); }
    std::array<float,8000> GetWaveformDerAt(int i) { return fSmoothArrayDer.at(i); }
    const std::array<std::array<int,8000>,160> & GetOrigWaveform() { return fOrigArray; }
    const std::array<std::array<float,8000>,160> & GetRawWaveform() { return fRawArray; }
    const std::array<std::array<float,8000>,160> & GetWaveform() { return fSmoothArray; }
    const std::array<std::array<float,8000>,160> & GetWaveformDer() { return fSmoothArrayDer; }

    size_t FindFirstAfter(float time, int startLoc = 0) const;
    const SinglePulse * FindPreviousPulse(int board, int channel);
    const SinglePulse * FindPreviousPulse(unsigned int key);

    Pulses & operator=(const Pulses & rhs);
    void Copy(const Pulses & rhs, int startBoard, int endBoard=10, int startChannel=0, int endChannel=15);

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
    std::array<std::array<int,8000>,160> fOrigArray; //!
    /// Array of the waveform after exponential moving average (not saved)
    std::array<std::array<float,8000>,160> fRawArray; //!
    /// Array of the waveform after 3 point average (not saved)
    std::array<std::array<float,8000>,160> fSmoothArray; //!
    /// Array of the derivative of the waveform after 3 point aveage (not saved)
    std::array<std::array<float,8000>,160> fSmoothArrayDer; //!

    /// Pointer to the current #SinglePulse (not saved)
    SinglePulse * fCurrPulse; //!



  ClassDef(Pulses,1)

};

#endif // #ifndef Pulses_h
