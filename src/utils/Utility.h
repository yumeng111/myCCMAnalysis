/*!**********************************************
 * \file Utility.h
 * \brief Header file for the #Utility class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#ifndef Utility_h
#define Utility_h

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <vector>
#include <array>
#include <map>
#include <cstring>
#include <numeric>
#include <iterator>
#include <utility>

/*!**********************************************
 * \enum HistInfo_t
 * \brief Enumeration object used in #MakeSPECalibrationHists class
 ***********************************************/
enum HistInfo_t
{			
	kRate2DHist,
	kSPE2DHist,
	kPEHists,
	kADC2DHist,
	kFlangeHists1,
	kFlangeHists2,
	kSPE1DHist,
	kRate1DHist,
	kNoisePeakHist,
	kNoiseEndHist,
	kFitEndHist,
	kFitRMSHist,
	kNoiseIntegralHist,
	kTailIntegralHist,
	kRootFile,
};

typedef enum {
  kRawDataID = 0, kPulsesID = 1, kEventsID = 2, kAccumWaveformID = 3, kMCTruthID = 4, kNEventBranch = 5
} CCMEventBranchID_t;

/*!**********************************************
 * \struct HighVoltage_t
 * \brief Contains information about the HV status of a PMT
 ***********************************************/
struct HighVoltage_t
{
  int status; ///< What is the status message from the CAEN high voltage supply
  int power; ///< What is the power from the CAEN high voltage supply

  /// \fn HighVoltage_t(ints, int p)
  /// \brief Constructor for this struct
  /// \param[in] s Value to set #status to
  /// \param[in] p Value to set #power to
  HighVoltage_t(int s, int p) : status(s), power(p) {}
};

/// \brief enum used to return the operation's result.
typedef enum {
  kCCMSuccess = 1, kCCMWarning = -1, kCCMError = -2, kCCMFailure = -3, kCCMDoNotWrite = 0, kCCMNewRun = -4
} CCMOldResult_t;
typedef int32_t CCMResult_t;

typedef enum {
  kCCMDynamicLengthEventID = 0,
  kCCMFixedLengthEventID = 1,
  kCCMNoneEventID = 1,
} CCMEventFinderID_t;

typedef enum {
    kCCMPulsesTimeID = 0,
    kCCMIntegralTimeID = 1,
    kCCMIntegralDerID = 2,
    kCCMVetoBottomTimeID = 3,
    kCCMVetoTopTimeID = 4,
    kCCMVetoCRightTimeID = 5,
    kCCMVetoCLeftTimeID = 6,
    kCCMVetoCFrontTimeID = 7,
    kCCMVetoCBackTimeID = 8,
    kCCMVetoTotalTimeID = 9,
    kCCMPMTWaveformID = 10,
    kCCMPMTWaveformCountID = 11
} CCMAccWaveform_t;

typedef enum {
  kCCMAccumWaveformTriangleID = 0,
  kCCMAccumWaveformStartID = 1,
  kCCMAccumWaveformTrianglePulseCutID = 2,
  kCCMAccumWaveformStartPulseCutID = 3,
  kCCMAccumWaveformTotalID = 4
} CCMAccumWaveformMethod_t;

namespace std {
  // I am not happy that I had to put these stream operators in std namespace.
  // I had to because otherwise std iterators cannot find them 
  // - you know this annoying C++ lookup rules...
  // I know one solution is to create new type inter-operable with this pair...
  // Just to lazy to do this - anyone knows workaround?
  extern istream& operator >> (istream& is, pair<int, double>& ps);
  extern ostream& operator << (ostream& os, const pair<const int, double>& ps);
}

/*!**********************************************
 * \namespace Utility
 * \brief A class with random algorithms that might be useful for different codes
 ***********************************************/
namespace Utility
{
  // Set the number of bins and bin width
  // This is hard coded and should be taken 
  // from either data_structures.hh or a data base
  // since they (in principle) could/will change
  
  /// The number of bins in the DAQ window
  constexpr const int    fgkNumBins = 8000;

  /// The width of each bin in the DAQ window in ns
  constexpr const double fgkBinWidth = 2.0;

  /// The number of digitizer channels that contains all the PMTs
  constexpr const double fgkNumPMTs = 160;

  /// The time in ns corresponding to the start of the DAQ window
  constexpr const double fgkWindowStartTime = -9920;

  /// The time in ns corresponding to the end of the DAQ window
  constexpr const double fgkWindowEndTime = fgkWindowStartTime + static_cast<double>(fgkNumBins)*fgkBinWidth;

  /// The time in bin number corresponding to the offset between CCM and FP3 detector
  constexpr const double fgkFP3Offset = 15;

  template<class T>
   void LinearUnweightedLS(size_t nPoints, T * x, T * y, T & par0, T & par1);

  template<class T>
   void FindQuartiles(std::vector<T> & vec, double & q1, double & q2, double & q3);

  template<class T>
   float Smooth(const T * vector, int start, int length);

  inline const char* tstamp()
  {
    //======================================================================
    // Provide a nicely formatted, current, time stamp string
    //======================================================================
    static char tbuff[32];
    time_t t;
    t = time(0);
    strcpy(tbuff, ctime(&t));
    tbuff[24] = '\0';
    return tbuff;
  }

  template<class T>
   T & Add(T & left, const T & right);

  extern void ParseStringForRunNumber(std::string name, int & run, int & subrun, 
      int * day = 0, int * month = 0, int * time = 0);

  template<typename T>
   int FindFirstNoneEmptyBin(typename std::vector<T>::iterator begin, 
        typename std::vector<T>::iterator start, typename std::vector<T>::iterator end);

  template<typename T>
   int FindFirstNoneEmptyBin(typename std::array<T,Utility::fgkNumBins>::iterator begin, 
        typename std::array<T,Utility::fgkNumBins>::iterator start, typename std::array<T,Utility::fgkNumBins>::iterator end);

  extern std::string ConvertCCMEventFinderIDToString(CCMEventFinderID_t evtFinder);
  extern CCMEventFinderID_t ConvertStringToCCMEventFinderID(std::string name);
  extern std::string ConvertCCMAccumWaveformMethodToString(CCMAccumWaveformMethod_t accumWaveform);
  extern CCMAccumWaveformMethod_t ConvertStringToCCMAccumWaveformMethod(std::string name);

  extern double ShiftTime(int start, int beamTime, bool applyFP3Offset = true);

  extern std::vector<std::string> GetListOfFiles(const char * file_regexp);
  extern void IndirectFileList(const char* file, std::vector<std::string>& infileList);

};

/*!**********************************************
 * \fn void Utility::LinearUnweightedLS(size_t nPoints, T * x, T * y, T & par0, T & par1)
 * \brief Calcualte a linear least squares
 * \param[in] nPoints The number of data points
 * \param[in] x The x value of those points
 * \param[in] y The y value of those points
 * \param[out] par0 The result for the 0th parameter
 * \param[out] par1 The result for the 1th parameter
 *
 * Assumes the fit is \f$y = par0 + x*par1\f$
 ***********************************************/
template<class T>
void Utility::LinearUnweightedLS(size_t nPoints, T * x, T * y, T & par0, T & par1)
{
  T sumX = 0.0;
  T sumX2 = 0.0;
  T sumY = 0.0;
  T sumXY = 0.0;

  for (size_t point = 0; point < nPoints; ++point) {
    sumX += x[point];
    sumX2 += std::pow(x[point],2.f);
    sumY += y[point];
    sumXY += x[point]*y[point];
  }

  par0 = (sumX2*sumY - sumXY*sumX)/(static_cast<T>(nPoints)*sumX2 - sumX*sumX);
  par1 = (static_cast<T>(nPoints)*sumXY - sumX*sumY)/(static_cast<T>(nPoints)*sumX2 - sumX*sumX);

  if (std::isnan(sumX) || std::isnan(sumX2) || std::isnan(sumY) || std::isnan(sumXY) || 
      std::isinf(sumX) || std::isinf(sumX2) || std::isinf(sumY) || std::isinf(sumXY) || 
      static_cast<T>(nPoints)*sumX2 - sumX*sumX == 0
      ) {
    sumX2 = 0.0;
    for (size_t point = 0; point < nPoints; ++point) {
      sumX2 += std::pow(x[point],2.f);
    }
  }

}

/*!**********************************************
 * \fn void Utility::FindQuartiles(std::vector<T> & vec, double & q1, double & q2, double & q3)
 * \brief Calcualte the quartiles given the data in \p vec
 * \param[in] vec A vector of the data to look at
 * \param[out] q1 First quartile 
 * \param[out] q2 Second quartile
 * \param[out] q3 Third quartile
 ***********************************************/
template<class T>
void Utility::FindQuartiles(std::vector<T> & vec, double & q1, double & q2, double & q3)
{
  size_t size = vec.size();

  if (size == 0) {
    return;
  }

  std::sort(vec.begin(), vec.end());

  int mid = size/2;
  q2 = size % 2 == 0 ? (vec[mid] + vec[mid-1])/2 : vec[mid];

  std::vector<T> first(vec.begin(),vec.begin()+mid);
  std::vector<T> third(vec.begin()+mid,vec.end());

  int side_length = 0;

  if (size % 2 == 0) 
  {
    side_length = size/2;
  }
  else {
    side_length = (size-1)/2;
  }

  q1 = (size/2) % 2 == 0 ? (first[side_length/2]/2 + first[(side_length-1)/2])/2 : first[side_length/2];
  q3 = (size/2) % 2 == 0 ? (third[side_length/2]/2 + third[(side_length-1)/2])/2 : third[side_length/2];

  return;
}

/*!**********************************************
 * \fn float Utility::Smooth(const T * vector, int start, int length)
 * \brief Smooth a point in a vector between \p start and \p start + \p length
 * \param[in] vector A vector of the data to look at
 * \param[in] start The start position in the vector
 * \param[in] length The number of indexes in the vector to look at
 * \return The smoothed value
 ***********************************************/
template<class T>
float Utility::Smooth(const T * vector, int start, int length)
{
  auto avg = std::accumulate(vector+start,vector+start+length,static_cast<T>(0));
  return static_cast<float>(avg)/static_cast<float>(length);
}

/*!**********************************************
 * \fn T & Utility::Add(T & left, const T & right)
 * \brief Add \p left and \p right and set the answer to \p left
 * \param[in,out] left The left argument of the addition equation
 * \param[in] right The right argument of the addition equation
 * \return The sum of \p left and \p right
 ***********************************************/
template<class T>
T & Utility::Add(T & left, const T & right)
{
  return left += right;
}

/*!**********************************************
 * \brief Loop the vector and find the first index that does not equal 0
 * \param[in] begin Iterator to the start of the vector
 * \param[in] start Iterator to where to start looking in the vector
 * \param[in] end Iterator to where to stop looping in the vector
 * \return The index of the first non-zero bin. Return value is -1 if all the bins equals zero
 ***********************************************/
template <class T>
int Utility::FindFirstNoneEmptyBin(typename std::vector<T>::iterator begin, 
    typename std::vector<T>::iterator start, typename std::vector<T>::iterator end)
{
  auto time = std::find_if_not(start,end,[](const T & i){return i == static_cast<T>(0);});
  if (std::distance(time,end) == 0) {
    return -1;
  }
  return std::distance(begin,time);
}

/*!**********************************************
 * \brief Loop the vector and find the first index that does not equal 0
 * \param[in] begin Iterator to the start of the vector
 * \param[in] start Iterator to where to start looking in the vector
 * \param[in] end Iterator to where to stop looping in the vector
 * \return The index of the first non-zero bin. Return value is -1 if all the bins equals zero
 ***********************************************/
template<typename T>
int Utility::FindFirstNoneEmptyBin(typename std::array<T,Utility::fgkNumBins>::iterator begin, 
    typename std::array<T,Utility::fgkNumBins>::iterator start, typename std::array<T,Utility::fgkNumBins>::iterator end)
{
  auto time = std::find_if_not(start,end,[](const T & i){return i == static_cast<T>(0);});
  if (std::distance(time,end) == 0) {
    return -1;
  }
  return std::distance(begin,time);
}

typedef std::map<int,std::vector<std::array<float,Utility::fgkNumBins>>> MapDAQWF2D;
typedef std::map<int,std::array<float,Utility::fgkNumBins>> MapDAQWF1D;

#endif // #ifndef Utility_h

