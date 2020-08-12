/*!************************************************************************************************
 * \file CCMBeamInfo.h
 * \brief Header file of the #CCMBeamInfo class
 * \author R. T. Thornton (LANL)
 * \date May 7, 2020
 **************************************************************************************************/
#ifndef CCMBeamInfo_h
#define CCMBeamInfo_h

#include <map>
#include <cstring>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>

typedef std::chrono::system_clock::time_point TimePoint;
typedef std::map<std::time_t,double> MapTimeDouble;


/*!************************************************************************************************
 * \class CCMBeamInfo
 * \brief Contains the beam current as a function of time
 **************************************************************************************************/
class CCMBeamInfo
{
  public:
    static void LoadTable(std::string fileName);
    static double GetCurrent(std::string time, std::string format = "%d-%b-%Y %X");
    static double GetCurrent(std::time_t time);
    static void PrintTable(std::string options);

    static std::time_t ConvertStringToTime(std::string time, std::string format = "%d-%b-%Y %X", bool print = false);
    static std::time_t ConvertToTime(int year, int mon, int day, int hour, int min, int sec, bool print = false);
    static std::time_t ConvertTimeSinceEPOCHtoTime(unsigned long int time);

    static void Reset();
    static int Next();
    static int Previous();
    static void Find(std::string time, std::string format = "%d-%b-%Y %X");
    static void Find(std::time_t time);

    static std::pair<std::time_t, double> GetBeamInfo();

  private:
    static MapTimeDouble fgCurrentInfo; ///< map to all the current as a function of time
    static std::tm fgLocalTM;

    static MapTimeDouble::iterator fgCurrentInfoIter; ///< iterator to current location in the map
    static MapTimeDouble::iterator fgCurrentInfoIterBegin; ///< iterator to the beginning of the map
    static MapTimeDouble::iterator fgCurrentInfoIterEnd; ///< iterator to the end of the map
};

#endif // CCMBeamInfo_h
