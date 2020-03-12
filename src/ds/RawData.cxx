/*!**********************************************
 * \file RawData.cxx
 * \brief Source code for the #RawData class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "RawData.h"

#include <algorithm>

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
  //fWaveforms.resize(numBoards*numChannels,std::vector<u_int16_t>(numSamples));
  fWaveforms.resize(numChannels,std::vector<u_int16_t>(numSamples)); // only save the last board
  fSize.resize(numBoards,std::vector<u_int16_t>(numChannels));
  fChannelMask.resize(numBoards,std::vector<u_int16_t>(numChannels));
  fTemp.resize(numBoards,std::vector<u_int16_t>(numChannels));
  fBoardEventNum.resize(numBoards);
  fClockTime.resize(numBoards);

}

/*!**********************************************
 * \fn RawData::RawData(const RawData & rhs)
 * \brief Copy constuctor
 * \param[in] rhs The object to copy
 ***********************************************/
RawData::RawData(const RawData & rhs) : TObject(rhs)
{
  fNumBoards       = rhs.fNumBoards;
  fNumChannels     = rhs.fNumChannels;
  fNumSamples      = rhs.fNumSamples;
  fEventNumber     = rhs.fEventNumber;
  fOffset          = rhs.fOffset;
  fWaveforms       = rhs.fWaveforms;
  fSize            = rhs.fSize;
  fChannelMask     = rhs.fChannelMask;
  fBoardEventNum   = rhs.fBoardEventNum;
  fClockTime       = rhs.fClockTime;
  fGPSNSIntoSec    = rhs.fGPSNSIntoSec;
  fGPSSecIntoDay   = rhs.fGPSSecIntoDay;
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

