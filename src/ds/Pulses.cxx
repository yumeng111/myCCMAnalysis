/*!**********************************************
 * \file Pulses.cxx
 * \brief Source code for the #Pulses class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "Pulses.h"
#include "Utility.h"
#include "MsgLog.h"
#include "PMTInfoMap.h"
#include "PMTInformation.h"

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <utility>

ClassImp(Pulses)

/*!**********************************************
 * \fn Pulses::Pulses(int board, int channel)
 * \brief Default constructor
 * \param[in] board The digitizer board currently looking at (default = 0)
 * \param[in] channel The channel on the digitizer board currently looking at (default = 0)
 ***********************************************/
Pulses::Pulses(int board, int channel) : TObject()
{
  fBoard = board;
  fChannel = channel;

  fEventNumber = 0;
  fComputerSecIntoEpoch = 0;
  fComputerNSIntoSec = 0;

  fCurrPulse = new SinglePulse();

  for (int i=0; i < 160; ++i) {
    for (int j=0; j < 8000; ++j) {
      //fOrigArray[i][j] = 0;
      fRawArray[i][j] = 0.f;
      fSmoothArray[i][j] = 0.f;
      //fSmoothArrayDer[i][j] = 0.f;
    }
  }

  Reset();
}

/*!**********************************************
 * \fn Pulses::Pulses(const Pulses & p)
 * \brief Copy constructor calles #operator=
 * \param[in] p The object to copy
 ***********************************************/
Pulses::Pulses(const Pulses & p) : TObject(p)
{
  this->operator=(p);
}

/*!**********************************************
 * \fn Pulses::~Pulses()
 * \brief Destructor calles #Reset
 ***********************************************/
Pulses::~Pulses()
{
  Reset();

  if (fCurrPulse) {
    delete fCurrPulse;
  }
}

/*!**********************************************
 * \fn void Pulses::Reset()
 * \brief Resets all the variables, vectors and arrays. Calls #ClearPulses
 ***********************************************/
void Pulses::Reset()
{
  fTriggerTime = 0.0;
  ClearPulses();
}

/*!**********************************************
 * \fn void Pulses::ClearPulses()
 * \brief Resets vectors and arrays associated with the pulses found
 ***********************************************/
void Pulses::ClearPulses()
{
  fNumPulses = 0;

  if (!fPulses.empty()) {
    fPulses.clear();
  }

  if (fCurrPulse) {
    fCurrPulse->Reset();
  }

  for (int i=0; i < 160; ++i) {
    for (int j=0; j < 8000; ++j) {
      //fOrigArray[i][j] = 0;
      fRawArray[i][j] = 0.f;
      fSmoothArray[i][j] = 0.f;
      //fSmoothArrayDer[i][j] = 0.f;
    }
  }

}

/*!**********************************************
 * \fn void Pulses::DerivativeFilter(const u_int16_t input[], int length, float triggerTime, int points, float threshold, float trigOffset, float pmtOffset, float adcToPE)
 * \brief Find pulses based on the derivative filter
 * \param[in] input The waveform for the current PMT
 * \param[in] length The number of samples in the current waveform
 * \param[in] triggerTime The minimum trigger time for the event (this gets saved for each waveform)
 * \param[in] points The number of points to use in the derivative (I believe this variable is not longer being used)
 * \param[in] threshold The threshold applied to the potential pulse that is found (0 means no threshold)
 * \param[in] trigOffset Is the time offset of the board
 * \param[in] pmtOffset Is the 1" to 8" pmt offset
 * \param[in] adcToPE The ADC to PE calibration to be applied (1.0 means no calibration)
 *
 * Find pulses by using the derivative of the waveform. Steps:
 * -# Smooth the waveform with an exponential moving filter
 *  \f[
 *  S_{t} = \left\{\begin{array}{ll}
 *   Y_1,                                           & t = 1 \\
 *   \alpha \cdot Y_t + (1 - \alpha) \cdot S_{t-1}, & t > 1
 *   \end{array}\right.
 *  \f]
 *  followed by a white noise filter (averages the current sample with the one before and after).
 *  The smoothing is crucial for getting rid of the various random noise in the system. Even though
 *  the random noise is a few ADC counts (< 5) on average, any small change messes up the derivative
 *  filter method. We are using the derivative filter because of the overshoot problem with the bases.
 *  The derivative filter does is not dependent on the baseline as long as the baseline changes slowly.
 * -# Apply the 3-point derivative to each sample in the waveform
 *  \f[
 *  S_{t} = 0.5\left(S_{t+1}-S_{t-1}\right)
 *  \f]
 * -# Find the pulse and save properties. The pulse starts when \f$S_{t}\f$ goes negative and ends 
 *  the second time \f$S_{t}\f$ goes negative. The means the derivative had to have gone positive
 *  before before going back to zero. The charge is approximated by the absolute integral of the derivative,
 *  which is approximately equal to the amplitude of the pulse, but without having to know the baseline.
 * -# An initial cut is applied to make sure not too much noise is saved. The cut requires:
 *  - The length of the pulse is at least 10 ns (5 samples).
 *  - The maximum absolute value of the derivative within the pulse has to be greater than or equal to 0.65
 * -# Apply the \p threshold passed and convert from ADCs to PE with \p adcToPE if the value does not equal 1.
 *
 * The pulses are added to the #fPulses vector with the rest of the pulses for the given DAQ window. That is why
 * it is important to call #SetBoardChannel before calling #DerivativeFilter, otherwise the pulses will be assigned
 * the wrong PMT.
 ***********************************************/
void Pulses::DerivativeFilter(const u_int16_t input[], int length, float triggerTime, 
    int points, float threshold, float trigOffset, float pmtOffset, float adcToPE)
{
  fTriggerTime = triggerTime;

  bool crossed1 = false;
  bool crossed2 = false;
  bool tempAbove = false;

  float time = -1.0;
  float pulseLength = 0.0;
  float integral = 0.0;
  float baseline = 0.0;
  float amplitude = 0.0;

  //float prevBaseline = std::numeric_limits<float>::max();
  //float prevEnd = std::numeric_limits<float>::max();

  int key = PMTInfoMap::CreateKey(fBoard,fChannel);
  fCurrPulse->SetKey(key);

  float f1 = 0.;

  float maxPos = -2;

  float currAmplitude = 0.0;

  fCurrPulse->Reset();
  fCurrPulse->SetLength(8000);
  //fCurrPulse->SetTime(2);
  fCurrPulse->SetWaveformStart(2);
  fCurrPulse->SetWaveformEnd(8000-4);
  fCurrPulse->SetBaseline(0);

  int timeOffset = 2;
  float alpha = 2e-9/10e-9;
  float smoothWidth = 0;//alpha*input[0];
  for (int loc = 0; loc < 4; ++loc) {
    //fOrigArray[key][loc] = input[loc];
    smoothWidth += alpha*(input[loc] - smoothWidth);
    //fSmoothArray.at(loc) = smoothWidth;
    //fRawArray.at(loc) = input[loc];
    fRawArray[key][loc] = smoothWidth;
  }
  //fSmoothArray.fill(0.f);
  fSmoothArray[key][1] = fRawArray[key][0] + fRawArray[key][1] + fRawArray[key][2];
  fSmoothArray[key][2] = fRawArray[key][1] + fRawArray[key][2] + fRawArray[key][3];

  float numSamples = 3.0;

  smoothWidth = 3.f;

  for (int sampleLoc=2; sampleLoc<length-1; ++sampleLoc) {
    if (sampleLoc + 2 < length - 1) {
      //fOrigArray[key][sampleLoc+2] = input[sampleLoc+2];
      fRawArray[key][sampleLoc+2] = fRawArray[key][sampleLoc+1]*(1.f  - alpha) + alpha*input[sampleLoc+2];
      fSmoothArray[key][sampleLoc+1] = fSmoothArray[key][sampleLoc] + fRawArray[key][sampleLoc+2] - fRawArray[key][sampleLoc-1];
    }

    f1 = (fSmoothArray[key][sampleLoc+1] - fSmoothArray[key][sampleLoc-1])/smoothWidth*0.5;
    //fSmoothArrayDer[key][sampleLoc] = f1;
    //MsgInfo(MsgLog::Form("sampleLoc %d input %d raw %.2f smooth %.2f der %.2f",sampleLoc,input[sampleLoc],fRawArray[key][sampleLoc],fSmoothArray[key][sampleLoc]/smoothWidth,fSmoothArrayDer[key][sampleLoc]));

    tempAbove = f1 <= 0;

    currAmplitude = baseline - fSmoothArray[key][sampleLoc]/smoothWidth;


    if (!crossed1 && !crossed2 && tempAbove) {
      fCurrPulse->Reset();
      baseline = 0.0;
      numSamples = 3.0;

      // calculate the baseline based on the three previous samples
      // save the waveform starting with five previous samples
      bool first = true;
      int sampleLoc2 = std::max(0,sampleLoc-5);
      numSamples = sampleLoc - sampleLoc2;
      for (; sampleLoc2 < sampleLoc; ++sampleLoc2) {
        if (first) {
          fCurrPulse->SetWaveformStart(sampleLoc2+timeOffset);
          first = false;
        }
        fCurrPulse->AddSample(fSmoothArray[key][sampleLoc2]/smoothWidth);
        baseline += fSmoothArray[key][sampleLoc2]/smoothWidth;
      }
      baseline /= numSamples;

      currAmplitude = baseline - fSmoothArray[key][sampleLoc]/smoothWidth;

      amplitude = currAmplitude;
      integral = -f1;

      time = sampleLoc;
      pulseLength= 1;
      maxPos = -f1;
      fCurrPulse->AddSample(fSmoothArray[key][sampleLoc]/smoothWidth);

      crossed1 = true;
    } else  if (crossed1 && !crossed2 && tempAbove) {
      // still with a negative first derivative and have not gone positive yet
      amplitude = std::max(currAmplitude,amplitude);
      integral += -f1;
      ++pulseLength;
      maxPos = std::max(-f1,maxPos);
      fCurrPulse->AddSample(fSmoothArray[key][sampleLoc]/smoothWidth);
    } else  if (crossed1 && !crossed2 && !tempAbove) {//currAmplitude > 0) {
      // first positive derivative after origionally going negative
      amplitude = std::max(currAmplitude,amplitude);
      integral += f1;
      ++pulseLength;
      crossed2 = true;
      maxPos = std::max(f1,maxPos);
      fCurrPulse->AddSample(fSmoothArray[key][sampleLoc]/smoothWidth);
    } else  if (crossed1 && crossed2 && !tempAbove) {//currAmplitude > 0) {
      // derivative has gone negative and then position again and the waveform is still
      // below the calculated baseline
      amplitude = std::max(currAmplitude,amplitude);
      maxPos = std::max(f1,maxPos);
      integral += f1;
      ++pulseLength;
      fCurrPulse->AddSample(fSmoothArray[key][sampleLoc]/smoothWidth);
    } else  if (crossed1 && tempAbove) {
      // first derivative went negative and now we are above the calculated baseline
      // after the derivative went positive 
      crossed1 = false;
      crossed2 = false;
      //prevBaseline = baseline;
      //prevEnd = static_cast<float>(fSmoothArray[key][sampleLoc-1]/smoothWidth);

      amplitude = std::max(currAmplitude,amplitude);
      maxPos = std::max(-f1,maxPos);
      integral += -f1;
      ++pulseLength;

      // save the pulse if 
      //      the integral is above threshold
      //      the pulse length is at least 8 samples (16 ns)
      //      the maximum derivative values (maxPos) is at least 0.65 
      if (integral >= threshold && pulseLength >= 5 && maxPos >= 0.65) {
        fCurrPulse->SetTriggerOffset(trigOffset);
        fCurrPulse->SetPMTOffset(pmtOffset);
        fCurrPulse->SetADCToPE(adcToPE);
        fCurrPulse->SetTime(time+trigOffset+pmtOffset);
        if (time+trigOffset+pmtOffset < 0 || fCurrPulse->GetTime() < 0) {
          MsgInfo(MsgLog::Form("Checking time of pulse time %g trigOffset %g pmtOffset %g",time,trigOffset,pmtOffset));
        }
        fCurrPulse->SetLength(pulseLength);
        fCurrPulse->SetMaxDerValue(maxPos);
        fCurrPulse->SetAmplitude(amplitude);
        fCurrPulse->SetIntegral(integral/adcToPE);
        fCurrPulse->SetBaseline(baseline);
        for (int sampleLoc2 = sampleLoc; sampleLoc2 < sampleLoc+5 && sampleLoc2 != length-4; ++sampleLoc2) {
          fCurrPulse->AddSample(fSmoothArray[key][sampleLoc2]/smoothWidth);
        }
        if (sampleLoc+5 < length-4) {
          fCurrPulse->SetWaveformEnd(sampleLoc);
        } else {
          fCurrPulse->SetWaveformEnd(length);
        }

        //if (!fPulses.empty()) {
        //  const SinglePulse * prevPulse = FindPreviousPulse(fBoard,fChannel);
        //  if (prevPulse) {
        //    if (fCurrPulse->GetTime() <= prevPulse->GetTime()+prevPulse->GetLength()+1) {
        //      fPulses.back().Append(*fCurrPulse);
        //    } else {
        //      fPulses.push_back(*fCurrPulse);
        //      ++fNumPulses;
        //    }
        //  } else {
        //    fPulses.push_back(*fCurrPulse);
        //    ++fNumPulses;
        //  }
        //} else {
          fPulses.push_back(SinglePulse(*fCurrPulse));
          ++fNumPulses;
        //}

      } // end if integral
      //--sampleLoc;
      //if (sampleLoc > 4960) {
      //  break;
      //}
    } else {
      crossed1 = false;
      crossed2 = false;
      //if (sampleLoc > 4960) {
      //  break;
      //}
    }
  }
}

/*!**********************************************
 * \fn void Pulses::MCFilter(const std::vector<double> & input, float triggerTime)
 * \brief Find pulses based on the derivative filter
 * \param[in] input The waveform for the current PMT
 * \param[in] triggerTime Offset from t = 0 if the true event time is different from 0
 *
 * Find pulses by using the following steps: (to be used until the base is simulated)
 * - The pulse starts when the waveform goes above 0 and ends when the waveform goes back to 0
 * - The integral is summing up all the bins between the start and the end
 *
 * The pulses are added to the #fPulses vector with the rest of the pulses for the given DAQ window. That is why
 * it is important to call #SetBoardChannel before calling #MCFilter, otherwise the pulses will be assigned
 * the wrong PMT.
 ***********************************************/
void Pulses::MCFilter(const std::vector<double> & input, float triggerTime)
{
  fTriggerTime = triggerTime;

  float time = -1.0;
  float pulseLength = 0.0;
  float integral = 0.0;
  float amplitude = 0.0;

  int key = PMTInfoMap::CreateKey(fBoard,fChannel);
  fCurrPulse->SetKey(key);

  float maxPos = -2;

  fCurrPulse->Reset();
  fCurrPulse->SetLength(8000);
  //fCurrPulse->SetTime(2);
  fCurrPulse->SetWaveformStart(2);
  fCurrPulse->SetWaveformEnd(8000-4);
  fCurrPulse->SetBaseline(0);

  const size_t kLength = input.size();
  for (size_t loc = 0; loc < kLength; ++loc) {
    if (input[loc] != 0 && integral == 0) {
      fCurrPulse->Reset();

      bool first = true;
      size_t sampleLoc2 = std::max(0,static_cast<int>(loc)-5);
      for (; sampleLoc2 < loc; ++sampleLoc2) {
        if (first) {
          fCurrPulse->SetWaveformStart(sampleLoc2+triggerTime/Utility::fgkBinWidth);
          first = false;
        }
        fCurrPulse->AddSample(input[loc]);
      }

      amplitude = input[loc];
      time = loc;
      pulseLength = 1;
      maxPos = input[loc]; 
      integral += input[loc];
      fCurrPulse->AddSample(input[loc]);

    } else  if (input[loc] != 0 && integral != 0) {
      amplitude = std::max(static_cast<float>(input[loc]),amplitude);
      integral += input[loc];
      ++pulseLength;
      maxPos = std::max(static_cast<float>(input[loc]),maxPos);
      fCurrPulse->AddSample(input[loc]);
    } else  if (input[loc] == 0 && integral != 0) {
      amplitude = std::max(static_cast<float>(input[loc]),amplitude);
      integral += input[loc];
      ++pulseLength;
      maxPos = std::max(static_cast<float>(input[loc]),maxPos);
      fCurrPulse->AddSample(input[loc]);

      // save the pulse if 
      //      the pulse length is at least 5 samples (10 ns)
      if (pulseLength >= 5) {
        fCurrPulse->SetTriggerOffset(triggerTime/Utility::fgkBinWidth);
        fCurrPulse->SetPMTOffset(0);
        fCurrPulse->SetADCToPE(1.0);
        fCurrPulse->SetTime(time+triggerTime/Utility::fgkBinWidth);
        if (time+triggerTime/Utility::fgkBinWidth < 0 || fCurrPulse->GetTime() < 0) {
          MsgInfo(MsgLog::Form("Checking time of pulse time %g trigTime",time,triggerTime/Utility::fgkBinWidth));
        }
        fCurrPulse->SetLength(pulseLength);
        fCurrPulse->SetMaxDerValue(maxPos);
        fCurrPulse->SetAmplitude(amplitude);
        fCurrPulse->SetIntegral(integral);
        fCurrPulse->SetBaseline(0);
        for (size_t sampleLoc2 = loc; sampleLoc2 < loc+5 && sampleLoc2 != kLength-4; ++sampleLoc2) {
          fCurrPulse->AddSample(input[loc]);
        }
        if (loc+5 < kLength-4) {
          fCurrPulse->SetWaveformEnd(loc);
        } else {
          fCurrPulse->SetWaveformEnd(kLength);
        }

        fPulses.push_back(SinglePulse(*fCurrPulse));
        ++fNumPulses;

        integral = 0.0;
        amplitude = 0.0;
        time = 0.0;
        maxPos = 0.0;
        pulseLength = 0.0;

      } else {
        integral = 0.0;
        amplitude = 0.0;
        time = 0.0;
        maxPos = 0.0;
        pulseLength = 0.0;
      } // end if integral
    }
  }
}

/*!**********************************************
 * \fn void Pulses::SmoothWaveform(const std::vector<u_int16_t> & input)
 * \brief Smooth the waveform and set results to #fSmoothArray
 * \param[in] input The waveform to smooth
 *
 * Smooths waveform by taking a 3 point average at a given sample.
 ***********************************************/
void Pulses::SmoothWaveform(const std::vector<u_int16_t> & input)
{
  int key = PMTInfoMap::CreateKey(fBoard,fChannel);

  float alpha = 2e-9/10e-9;
  float smoothWidth = 0;//alpha*input[0];
  for (int loc = 0; loc < 4; ++loc) {
    //fOrigArray[key][loc] = input[loc];
    smoothWidth += alpha*(input[loc] - smoothWidth);
    //fSmoothArray.at(loc) = smoothWidth;
    //fRawArray.at(loc) = input[loc];
    fRawArray[key][loc] = smoothWidth;
  }
  //fSmoothArray.fill(0.f);
  fSmoothArray[key][1] = fRawArray[key][0] + fRawArray[key][1] + fRawArray[key][2];
  fSmoothArray[key][2] = fRawArray[key][1] + fRawArray[key][2] + fRawArray[key][3];

  smoothWidth = 3.f;

  int length = static_cast<int>(input.size());

  for (int sampleLoc=2; sampleLoc<length-1; ++sampleLoc) {
    if (sampleLoc + 2 < length - 1) {
      //fOrigArray[key][sampleLoc+2] = input[sampleLoc+2];
      fRawArray[key][sampleLoc+2] = fRawArray[key][sampleLoc+1]*(1.f  - alpha) + alpha*input[sampleLoc+2];
      fSmoothArray[key][sampleLoc+1] = fSmoothArray[key][sampleLoc] + fRawArray[key][sampleLoc+2] - fRawArray[key][sampleLoc-1];
    }
  }
  for (auto & value : fSmoothArray[key]) {
    value /= smoothWidth;
  }
}


/*!**********************************************
 * \fn bool Pulses::SortCondition(const SinglePulse & a, const SinglePulse & b)
 * \brief Condition to check for sorting the pulses vector
 * \param[in] a First #SinglePulse
 * \param[in] b Second #SinglePulse
 * \return True of the time of a is less than the time of b
 *
 ***********************************************/
bool Pulses::SortCondition(const SinglePulse & a, const SinglePulse & b)
{
  return a.GetTime() < b.GetTime();
}

/*!**********************************************
 * \fn void Pulses::Sort()
 * \brief The sort function for the pulses vector
 *
 * By default the pulses vector is sorted by PMT id.
 * This function will sort the vector based on the
 * condition specified in #SortCondition.
 *
 ***********************************************/
void Pulses::Sort()
{
  std::sort(fPulses.begin(),fPulses.end(),Pulses::SortCondition);
}

/*!**********************************************
 * \fn size_t Pulses::FindFirstAfter(float time, int startLoc) const
 * \brief Return which pulse is the first pulse with a time after \p time.
 * \param[in] time The time condition
 * \param[in] startLoc Where in the vector to start looking
 * \return The distance from the start of the vector to the first pulse after \p time
 *
 * Loop through the pulses vector starting at \p startLoc and find the first pulse
 * that has a start time greater than or equal to time. If \p startLoc is 0 then
 * the loop will start at the beginning. Funny results will be given if the #Sort 
 * function is not called first.
 *
 ***********************************************/
size_t Pulses::FindFirstAfter(float time, int startLoc) const
{
  std::vector<SinglePulse>::const_iterator it = std::find_if(fPulses.cbegin()+startLoc,fPulses.cend(),
      [time](const SinglePulse & pulse) { 
      return pulse.GetTime() >= time; }
      );
  return std::distance(fPulses.cbegin(),it);
}

/*!**********************************************
 * \fn float Pulses::LongestPulse(float time) const
 * \brief Find the longest pulse that happened before \p time
 * \param[in] time The time condition
 * \return The length of the longest pulse
 *
 * Unlike #FindFirstAfter this function does not assume the vector
 * of pulses has been sorted by time
 ***********************************************/
float Pulses::LongestPulse(float time) const
{
  float longestPulse = -1;
  for (const auto & pulse : fPulses) {
    if (pulse.GetTime() >= time) {
      continue;
    }
    longestPulse = std::max(longestPulse,pulse.GetLength());
  }
  return longestPulse;
}

/*!**********************************************
 * \fn float Pulses::LargestPulse(float time) const
 * \brief Find the largest pulse that happened before \p time
 * \param[in] time The time condition
 * \return The integral of the largest pulse
 *
 * Unlike #FindFirstAfter this function does not assume the vector
 * of pulses has been sorted by time
 ***********************************************/
float Pulses::LargestPulse(float time) const
{
  float largestPulse = -1;
  for (const auto & pulse : fPulses) {
    if (pulse.GetTime() >= time) {
      continue;
    }
    largestPulse = std::max(largestPulse,pulse.GetIntegral());
  }
  return largestPulse;
}

/*!**********************************************
 * \fn const SinglePulse * Pulses::FindPreviousPulse(int board, int channel)
 * \brief Find the last pulse associated with a given PMT
 * \param[in] board Digitizer board
 * \param[in] channel Channel on the digitizer board
 * \return The pointer to the #SinglePulse of the last pulse associated with the given PMT
 ***********************************************/
const SinglePulse * Pulses::FindPreviousPulse(int board, int channel)
{
  int key = PMTInfoMap::CreateKey(board, channel);
  return FindPreviousPulse(key);
}

/*!**********************************************
 * \fn const SinglePulse * Pulses::FindPreviousPulse(int key)
 * \brief Find the last pulse associated with a given PMT
 * \param[in] key The PMT unique key
 * \return The pointer to the #SinglePulse of the last pulse associated with the given PMT
 ***********************************************/
const SinglePulse * Pulses::FindPreviousPulse(unsigned int key)
{
  std::vector<SinglePulse>::const_reverse_iterator itBegin = fPulses.rbegin();
  std::vector<SinglePulse>::const_reverse_iterator itEnd = fPulses.rend();
  std::vector<SinglePulse>::const_reverse_iterator it = itBegin;

  for (; it != itEnd; ++it) {
    if ((*it).GetKey() == key) {
      return &(*it);
    }
  }

  return nullptr;
}

/*!**********************************************
 * \fn Pulses & Pulses::operator=(const Pulses & rhs) 
 * \brief Set the current #Pulses object to the one passed
 * \param[in] rhs The object to copy
 ***********************************************/
Pulses & Pulses::operator=(const Pulses & rhs) 
{
  this->fBoard = rhs.fBoard;
  this->fChannel = rhs.fChannel;
  this->fEventNumber = rhs.fEventNumber;
  this->fComputerSecIntoEpoch= rhs.fComputerSecIntoEpoch;
  this->fComputerNSIntoSec= rhs.fComputerNSIntoSec;
  this->fNumPulses = rhs.fNumPulses;
  this->fTriggerTime = rhs.fTriggerTime;
  this->fPulses = rhs.fPulses;
  this->fCurrPulse = rhs.fCurrPulse;
  return *this;
}

/*!**********************************************
 * \fn void Pulses::CopyPulses(const Pulses & rhs, int startBoard, int endBoard, int startChannel, int endChannel) 
 * \brief Set the current #Pulses object to the one passed, but only copy a range of the board and channels
 * \param[in] rhs The object to copy
 * \param[in] startBoard The board to start copying
 * \param[in] endBoard The board to end copying
 * \param[in] startChannel The channel to start copying
 * \param[in] endChannel The endChannel to start copying
 ***********************************************/
void Pulses::CopyPulses(const Pulses & rhs, int startBoard, int endBoard, int startChannel, int endChannel) 
{
  this->fBoard = rhs.fBoard;
  this->fChannel = rhs.fChannel;
  this->fEventNumber = rhs.fEventNumber;
  this->fComputerSecIntoEpoch= rhs.fComputerSecIntoEpoch;
  this->fComputerNSIntoSec= rhs.fComputerNSIntoSec;
  this->fNumPulses = rhs.fNumPulses;
  this->fTriggerTime = rhs.fTriggerTime;
  this->fCurrPulse = rhs.fCurrPulse;

  int board = 0;
  int channel = 0;
  for (auto & pulse : rhs.fPulses) {
    int key = pulse.GetKey();
    PMTInfoMap::DecodeKey(key,board,channel);
    if (board < startBoard || board > endBoard || channel < startChannel || channel > endChannel) {
      continue;
    }
    fPulses.push_back(pulse);
  }
}

/*!**********************************************
 * \fn void Pulses::RemovePulsesByThreshold()
 * \brief Loop through the pulses vector and remove pulses below threshold
 * 
 * Loop through the pulses vector and remove pulses below threshold.
 * <b>This function is not used and should be modified before it is used.
 * It currently grabs a hard coded calibration file to get the thresholds.
 * Need to make it so a vector of thresholds is passed, otherwise the function 
 * works great.</b>
 ***********************************************/
void Pulses::RemovePulsesByThreshold()
{
  TFile * calibrationFile = TFile::Open("root_out_2019ledData_run118_.root","READ");
  std::map<int,std::pair<float,float>> calibrationValues;
  TTreeReader calibrationTree("spe",calibrationFile);
  TTreeReaderValue<float> speValue(calibrationTree,"speValue");
  TTreeReaderValue<float> speThreshold(calibrationTree,"endNoiseWallFitRangeStart");
  TTreeReaderValue<int> pmtID(calibrationTree,"pmtID");
  while (calibrationTree.Next()) {
    calibrationValues.insert(std::make_pair(*pmtID,std::make_pair(*speValue,*speThreshold)));
  }
  delete calibrationFile;

  fPulses.erase(
      std::remove_if(
        fPulses.begin(),
        fPulses.end(),
        [&calibrationValues] (SinglePulse const & p) {
        //const PMTInformation * pmtInfo = PMTInfoMap::GetPMTInfo(p.GetKey());
        return p.GetIntegral() < calibrationValues[p.GetKey()].second; 
        }
        )
      );

  std::for_each(fPulses.begin(),fPulses.end(),[&calibrationValues] (SinglePulse & p) {
      p.SetIntegral(p.GetIntegral()/calibrationValues[p.GetKey()].first);
      }
      );
}

/*!**********************************************
 * \fn void Pulses::ShiftTimeOffset(const double & timeOffset)
 * \brief Shift the time of all the pulses
 * \param[in] timeOffset The time offset to apply to all pulses
 *
 * Shift the time of all the pulses by timeOffset, since the stored time
 * in the pulses is currently the starting bin,
 * timeOffset needs to be in number of bins and will be rounded
 * down to the closed bin
 ***********************************************/
void Pulses::ShiftTimeOffset(const double & timeOffset)
{
  double binOffset = std::floor(timeOffset);

  fTriggerTime += binOffset;

  for (auto & p : fPulses) {
    p.SetTime(p.GetTime()+binOffset);
    p.SetTriggerOffset(fTriggerTime);
  }


  return;
}

