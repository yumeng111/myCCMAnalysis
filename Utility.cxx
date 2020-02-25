/*!**********************************************
 * \file Utility.cxx
 * \brief Source code for the #Utility class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "Utility.h"
#include "RawData.h"

/*!**********************************************
 * \fn int Utility::FindFirstSample(int channelNumber, const RawData * rawData)
 * \brief Find the first sample above threshold for a digitizer channel
 * \param[in] channelNumber The channel number based on the total number of channels saved in the #RawData class
 * \param[in] rawData Pointer to the current #RawData information to look at
 *
 * Assumes the pulse that is possibly in the window is a NIM signal
 ***********************************************/
int Utility::FindFirstSample(int channelNumber, const RawData * rawData)
{
  int firstSample = 0;
  auto samples = rawData->GetSamples(channelNumber);
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

