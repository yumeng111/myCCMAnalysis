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
#include <cstring>
#include <numeric>

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
  kRawDataID = 0, kPulsesID = 1, kEventsID = 2, kNEventBranch = 3
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

/*!**********************************************
 * \class Utility
 * \brief A class with random algorithms that might be useful for different codes
 ***********************************************/
class Utility
{
  public:
    template<class T>
    static void LinearUnweightedLS(size_t nPoints, T * x, T * y, T & par0, T & par1);

    template<class T>
      static void FindQuartiles(std::vector<T> & vec, double & q1, double & q2, double & q3);

    template<class T>
      static float Smooth(const T * vector, int start, int length);

    static const char* tstamp()
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
    static T & Add(T & left, const T & right);

    static void ParseStringForRunNumber(std::string name, int & run, int & subrun);

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

#endif // #ifndef Utility_h

