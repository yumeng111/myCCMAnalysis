/*!**********************************************
 * \file RawData.cxx
 * \brief Source code for the #RawData class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "RawData.h"
#include "Pulses.h"
#include "MsgLog.h"

#include <memory>
#include <algorithm>
#include <numeric>

ClassImp(RawData)

/*!**********************************************
 * \fn RawData::RawData()
 * \brief Default Constructor
 ***********************************************/
RawData::RawData() : TObject()
{
  fNumBoards = 11; // default settings
  fNumChannels = 16; // default settings
  fNumSamples = 8000; // default settings
  fEventNumber = 0;
  fOffset.resize(0);
  fWaveforms.resize(0);
  fSize.resize(0);
  fChannelMask.resize(0);
  fBoardEventNum.resize(0);
  fClockTime.resize(0);
  fGPSNSIntoSec = 0;
  fGPSSecIntoDay = 0;
}

/*!**********************************************
 * \fn RawData::RawData(size_t numBoards, size_t numChannels, size_t numSamples, u_int32_t evtNum, u_int32_t secToDay, u_int32_t nsToSec)
 * \brief Constructor (calls #Reset)
 * \param[in] numBoards The number of digitizer boards
 * \param[in] numChannels The number of channels on a digitizer board
 * \param[in] numSamples The number of samples for a given waveform
 * \param[in] evtNum The trigger number of the DAQ window
 * \param[in] secToDay The time of the DAQ window in s
 * \param[in] nsToSec The time of the DAQ window in ns
 ***********************************************/
RawData::RawData(size_t numBoards, size_t numChannels, size_t numSamples, u_int32_t evtNum, u_int32_t secToDay, u_int32_t nsToSec)
: TObject()
{
  Reset(numBoards,numChannels,numSamples,evtNum,secToDay,nsToSec);
}

/*!**********************************************
 * \fn void RawData::Reset(size_t numBoards, size_t numChannels, size_t numSamples, u_int32_t evtNum, u_int32_t secToDay, u_int32_t nsToSec)
 * \brief Resets the values to what is passed
 * \param[in] numBoards The number of digitizer boards
 * \param[in] numChannels The number of channels on a digitizer board
 * \param[in] numSamples The number of samples for a given waveform
 * \param[in] evtNum The trigger number of the DAQ window
 * \param[in] secToDay The time of the DAQ window in s
 * \param[in] nsToSec The time of the DAQ window in ns
 ***********************************************/
void RawData::Reset(size_t numBoards, size_t numChannels, size_t numSamples, u_int32_t evtNum, u_int32_t secToDay, u_int32_t nsToSec)
{
  fNumBoards       = numBoards;
  fNumChannels     = numChannels;
  fNumSamples      = numSamples;
  fEventNumber     = evtNum;
  fGPSNSIntoSec    = nsToSec;
  fGPSSecIntoDay   = secToDay;

  fOffset.resize(numBoards);
  fWaveforms.resize(numBoards*numChannels,std::vector<u_int16_t>(numSamples));
  fSize.resize(numBoards,std::vector<u_int16_t>(numChannels));
  fChannelMask.resize(numBoards,std::vector<u_int16_t>(numChannels));
  fTemp.resize(numBoards,std::vector<u_int16_t>(numChannels));
  fBoardEventNum.resize(numBoards);
  fClockTime.resize(numBoards);

}

/*!**********************************************
 * \fn void RawData::TruncateWaveform(size_t numBoards)
 * \brief Truncates the number of waveforms saved to the numBoards passed * #fNumChannels
 * \param[in] numBoards The number of digitizer boards to save
 ***********************************************/
void RawData::TruncateWaveform(size_t numBoards)
{
  fWaveforms.resize(numBoards*fNumChannels,std::vector<u_int16_t>(fNumSamples));
}

/*!**********************************************
 * \fn RawData::RawData(const RawData & rhs)
 * \brief Copy constuctor
 * \param[in] rhs The object to copy
 ***********************************************/
RawData::RawData(const RawData & rhs) : TObject(rhs)
{
  this->operator=(rhs);
}

/*!**********************************************
 * \fn RawData::~RawData()
 * \brief Destructor
 ***********************************************/
RawData::~RawData()
{
  // empty
}

//-------------------------------------------------------------------------------------------------
void RawData::SetWaveform(int detector, const u_int16_t input[])
{
  fWaveforms[detector].assign(input,input+fNumSamples);
}

//-------------------------------------------------------------------------------------------------
void RawData::SetOffset(const std::vector<float> & offsets)
{
  fOffset = offsets;
}

//-------------------------------------------------------------------------------------------------
void RawData::SetSize(size_t board, const u_int16_t input[])
{
  fSize[board].assign(input,input+fNumChannels);
}

//-------------------------------------------------------------------------------------------------
void RawData::SetChannelMask(size_t board, const u_int16_t input[])
{
  fChannelMask[board].assign(input,input+fNumChannels);
}

//-------------------------------------------------------------------------------------------------
u_int16_t RawData::SetTemp(size_t board, const u_int16_t input[])
{
  fTemp[board].assign(input,input+fNumChannels);
  return *std::max_element(input,input+fNumChannels);
}

//-------------------------------------------------------------------------------------------------
void RawData::SetBoardEventNum(const u_int32_t input[])
{
  fBoardEventNum.assign(input,input+fNumBoards);
}

//-------------------------------------------------------------------------------------------------
void RawData::SetClockTime(const u_int32_t input[])
{
  fClockTime.assign(input,input+fNumBoards);
}

/*!**********************************************
 * \fn int RawData::FindFirstNIMSample(int channelNumber)
 * \brief Find the first sample above threshold for a digitizer channel
 * \param[in] channelNumber The channel number based on the total number of channels saved in #fWaveforms
 * \return The sample number to the turn on of the NIM signal
 *
 * Assumes the pulse that is possibly in the window is a NIM signal
 ***********************************************/
int RawData::FindFirstNIMSample(int channelNumber)
{
  int firstSample = 0;
  auto samples = GetSamples(channelNumber);
  //size_t lengthOfVector = samples.size();
  //MsgInfo(MsgLog::Form("Length of vector = %zu",lengthOfVector));
  double avgChannelNoise = std::accumulate(samples.begin(),samples.begin()+500,0);
  avgChannelNoise /= 500.0;
  int sample = 0;
  for (const auto & adc : samples) {
    double adjustedADC = avgChannelNoise- static_cast<double>(adc);
    //std::cout << "\t\t\t" << channelNumber << '\t' << sample << '\t' << adc << '\t' << adjustedADC << std::endl;
    if (adjustedADC > 400) {
      firstSample = sample;
      break;
    }
    ++sample;
  } // end for sample < 8000

  return firstSample;
}

/*!**********************************************
 * \fn bool RawData::IsTriggerPresent(std::string tirggerName)
 * \brief Return true if the trigger is present
 * \param[in] triggerName The name of the trigger
 * \return True if the trigger is present
 *
 * The \p triggerName can be any combination of
 * - BEAM
 * - STROBE
 * - LED
 * - ALL
 ***********************************************/
bool RawData::IsTriggerPresent(std::string triggerName)
{
  int boardOffset = GetBoard10ChannelOffset();
  if (boardOffset < 0) {
    return false;
  }

  if (MsgLog::GetGlobalDebugLevel() >= 5) {
    MsgDebug(5,MsgLog::Form("Board Offset = %d",boardOffset));
  }

  int firstSampleStrobe = 0;
  int firstSampleBeam = 0;
  int firstSampleLED = 0;
  
  //WT - WORKAROUND TO ALLOW MODULES TO WORK W/OUT TRIGGERS ON LAST BOARD
   if ( triggerName.find("ALL",0) != std::string::npos){
        return true;
      }
   
  if (triggerName.find("STROBE") != std::string::npos ||
      triggerName.find("ALL") != std::string::npos) {
    firstSampleStrobe = FindFirstNIMSample(boardOffset+1);
  } 
  if (triggerName.find("BEAM") != std::string::npos ||
      triggerName.find("ALL") != std::string::npos) {
    firstSampleBeam = FindFirstNIMSample(boardOffset+2);
  } 
  if (triggerName.find("LED") != std::string::npos ||
      triggerName.find("ALL") != std::string::npos) {
    firstSampleLED = FindFirstNIMSample(boardOffset+3);
  } 

  if (firstSampleStrobe == 0 && firstSampleBeam == 0 && firstSampleLED == 0) {
    if (MsgLog::GetGlobalDebugLevel() >= 1) {
      MsgDebug(1,MsgLog::Form("Trigger name is %s and could not be found. Skipping Event", triggerName.c_str()));
    }
    return false;
  }

  return true;

}

/*!**********************************************
 * \fn float RawData::GetEarliestOffset()
 * \brief Return the time of the earliest board copy of the trigger
 * \return Return the time of the earliest board copy of the trigger
 ***********************************************/
float RawData::GetEarliestOffset()
{
  return *std::min_element(fOffset.begin(),fOffset.end());
}

/*!**********************************************
 * \fn int RawData::GetBCMTime(double * integral, double * length)
 * \brief Return the sample number to the beam current monitor turn on
 * \param[out] integral The integral of the BCM signal
 * \param[out] length The length of the BCM signal
 * \return The sample number to the beam current monitor turn on
 *
 * Uses the derivative of the waveform to find the turn on.
 * Function dependent of the Pulses class
 ***********************************************/
int RawData::GetBCMTime(double * integral, double * length)
{
  int boardOffset = GetBoard10ChannelOffset();

  double offset = GetOffset(10) - GetEarliestOffset();
  std::shared_ptr<Pulses> pulsesBeam = std::make_shared<Pulses>(0,0);
  pulsesBeam->DerivativeFilter(&fWaveforms.at(boardOffset+0).front(),GetNumSamples(),0,3,2800,offset,1.0);

  size_t beamTime = GetNumSamples()+1;
  const size_t kNumPulses = pulsesBeam->GetNumPulses();
  for (size_t pulse = 0; pulse < kNumPulses; ++pulse) {
    auto time = pulsesBeam->GetPulseTime(pulse);
    if (beamTime == GetNumSamples()+1) {
      beamTime = time;
      if (integral) {
        *integral = pulsesBeam->GetPulseIntegral(pulse);
      }
      if (length) {
        *length = pulsesBeam->GetPulseLength(pulse);
      }
      break;
    }
  }
  pulsesBeam->Reset();

  return beamTime;

}

/*!**********************************************
 * \fn int RawData::GetBoard10ChannelOffset()
 * \brief Return the number of channels board 10 starts at in the #fWaveforms vector
 * \return Return the number of channels board 10 starts at in the #fWaveforms vector
 *
 * Sometimes the data is saved with all the boards and sometimes with just the
 * 10 board has its sampels saved. This function returns the number of channels
 * to offset the index in #fWaveforms by to get the right waveform
 ***********************************************/
int RawData::GetBoard10ChannelOffset()
{
  int boardOffset = 0;
  if (GetNumBoards() == 0) {
    MsgWarning(MsgLog::Form("Something went wrong in that the number of boards = 0, waveform size %zu", GetNumWaveforms()));
    return -1;
  } else {
    //MsgInfo(MsgLog::Form("NumWaveforms %zu NumChannels %zu",GetNumWaveforms(),GetNumChannels()));
    if (GetNumWaveforms() != GetNumChannels()) {
      boardOffset = (GetNumBoards()-1)*GetNumChannels();
    }
  }

  return boardOffset;
}

/*!**********************************************
 * \fn RawData & RawData::operator=(const RawData & rhs)
 * \brief Copy assignment for the #RawData class
 * \return Return reference to the current class
 ***********************************************/
RawData & RawData::operator=(const RawData & rhs)
{
  this->fNumBoards       = rhs.fNumBoards;
  this->fNumChannels     = rhs.fNumChannels;
  this->fNumSamples      = rhs.fNumSamples;
  this->fEventNumber     = rhs.fEventNumber;
  this->fOffset          = rhs.fOffset;
  this->fWaveforms       = rhs.fWaveforms;
  this->fSize            = rhs.fSize;
  this->fChannelMask     = rhs.fChannelMask;
  this->fBoardEventNum   = rhs.fBoardEventNum;
  this->fClockTime       = rhs.fClockTime;
  this->fGPSNSIntoSec    = rhs.fGPSNSIntoSec;
  this->fGPSSecIntoDay   = rhs.fGPSSecIntoDay;

  return *this;
}


