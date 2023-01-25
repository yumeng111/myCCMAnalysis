/*!**********************************************
 * \file RawData_v2.h
 * \brief Header file for the #RawData_v2 class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef RawData_v2_h
#define RawData_v2_h

#include <vector>
#include <iostream>
#include <stdlib.h>

static const unsigned legacy_raw_data_version_ = 2;

#include "TObject.h"

/*!**********************************************
 * \class RawData_v2
 * \brief Container for the raw data coming out for the detector
 ***********************************************/
class RawData_v2 : public TObject
{
  public:
    RawData_v2();
    RawData_v2(size_t numBoards, size_t numChannels, size_t numSamples, u_int32_t evtNum, u_int32_t secToDay, u_int32_t nsToSec);
    RawData_v2(const RawData_v2 & rhs);
    ~RawData_v2();

    void Reset(size_t numBoards, size_t numChannels, size_t numSamples, u_int32_t evtNum, u_int32_t secToDay,u_int32_t nsToSec);
    void TruncateWaveform(size_t numberBoards = 1);

    void SetWaveform(int detector, const u_int16_t input[]);
    void SetOffset(const std::vector<float> & offsets);
    void SetSize(size_t board, const u_int16_t input[]);
    void SetChannelMask(size_t board, const u_int16_t input[]);
    void SetBoardEventNum(const u_int32_t input[]);
    u_int16_t SetTemp(size_t board, const u_int16_t input[]);
    void SetClockTime(const u_int32_t input[]);

    float GetOffset(int detector) { return fOffset[detector]; }
    std::vector<u_int16_t> & GetSamples(int detector) { return fWaveforms[detector]; }
    u_int16_t GetSample(int detector, int sample) {
        return fWaveforms[detector][sample];
    }
    u_int16_t GetSize(int board, int channel) { return fSize[board][channel]; }
    u_int16_t GetChannelMask(int board, int channel) { return fChannelMask[board][channel]; }
    u_int16_t GetTemp(int board, int channel) { return fTemp[board][channel]; }
    u_int32_t GetBoardEventNum(int board) { return fBoardEventNum[board]; }
    u_int32_t GetClockTime(int board) { return fClockTime[board]; }

    u_int32_t GetEventNumber() { return fEventNumber; }
    u_int32_t GetGPSNSIntoSec() { return fGPSNSIntoSec; }
    u_int32_t GetGPSSecIntoDay() { return fGPSSecIntoDay; }


    size_t GetNumBoards() { return fSize.size(); }
    size_t GetNumChannels() { return fSize.front().size(); }
    size_t GetNumSamples() { return fWaveforms.front().size(); }
    size_t GetNumWaveforms() { return fWaveforms.size(); }

    float GetOffset(int detector) const { return fOffset[detector]; }
    const std::vector<u_int16_t> & GetSamples(int detector) const { return fWaveforms[detector]; }
    u_int16_t GetSample(int detector, int sample) const { return fWaveforms[detector][sample]; }
    u_int16_t GetSize(int board, int channel) const { return fSize[board][channel]; }
    u_int16_t GetChannelMask(int board, int channel) const { return fChannelMask[board][channel]; }
    u_int16_t GetTemp(int board, int channel) const { return fTemp[board][channel]; }
    u_int32_t GetBoardEventNum(int board) const { return fBoardEventNum[board]; }
    u_int32_t GetClockTime(int board) const { return fClockTime[board]; }

    u_int32_t GetEventNumber() const { return fEventNumber; }
    u_int32_t GetGPSNSIntoSec() const { return fGPSNSIntoSec; }
    u_int32_t GetGPSSecIntoDay() const { return fGPSSecIntoDay; }


    size_t GetNumBoards() const { return fSize.size(); }
    size_t GetNumChannels() const { return fSize.front().size(); }
    size_t GetNumSamples() const { return fWaveforms.front().size(); }

    float GetEarliestOffset();
    int FindFirstNIMSample(int channelNumber);
    bool IsTriggerPresent(std::string triggerName);
    int GetBCMTime(double * integral = 0, double * length = 0);
    int GetBoard10ChannelOffset();

    RawData_v2 & operator=(const RawData_v2 & rhs);

  private:

    /// Number of digitizer boards
    size_t fNumBoards;
    /// Number of channels on a digitizer board
    size_t fNumChannels;
    /// Number of samples for a given waveform
    size_t fNumSamples;

    /// DAQ Window unique ID
    u_int32_t fEventNumber;
    /// Times the trigger was observed in eqch board
    std::vector<float> fOffset;
    /// Waveform for each PMT
    std::vector<std::vector<u_int16_t>> fWaveforms;
    /// The size of each waveform for each PMT by board and channel
    std::vector<std::vector<u_int16_t>> fSize;
    /// The channel mask for each PMT
    std::vector<std::vector<u_int16_t>> fChannelMask;
    /// The temperature of each digitizer channel
    std::vector<std::vector<u_int16_t>> fTemp;
    /// The event number on each board
    std::vector<u_int32_t> fBoardEventNum;
    /// The clock time of each board
    std::vector<u_int32_t> fClockTime;
    /// DAQ window time in ns
    u_int32_t fGPSNSIntoSec;
    /// DAQ window time in s
    u_int32_t fGPSSecIntoDay;

};

#endif // RawData_v2_h
