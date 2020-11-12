/*!************************************************************************************************
 * \file CCMBeamInfo.cxx
 * \brief Functions of the #CCMBeamInfo class
 * \author R. T. Thornton (LANL)
 * \date May 7, 2020
 **************************************************************************************************/
#include "CCMBeamInfo.h"
#include "MsgLog.h"

#include <fstream>
#include <iterator>
#include <sstream>
#include <iomanip>

MapTimeDouble CCMBeamInfo::fgCurrentInfo;
std::tm CCMBeamInfo::fgLocalTM = {0,0,0,1,1,119,0,0,1};

MapTimeDouble::iterator CCMBeamInfo::fgCurrentInfoIter;
MapTimeDouble::iterator CCMBeamInfo::fgCurrentInfoIterBegin;
MapTimeDouble::iterator CCMBeamInfo::fgCurrentInfoIterEnd;

/*!************************************************************************************************
 * \fn void CCMBeamInfo::LoadTable(std::string fileName)
 * \brief Loads the beam current information from  the file into #fgCurrentInfo
 * \param[in] fileName The name of the file to parse
 **************************************************************************************************/
void CCMBeamInfo::LoadTable(std::string fileName)
{
  std::string day = ""; // day of the year
  std::string time = ""; // time of day
  double amount = 0.0; // 1000*current
  double current = 0.0; // the value of the current
  std::string units = ""; // the units the corrent is in
  double last = 0; // do not know what this is used for

  std::string tempString = "";


  std::ifstream infile(fileName.c_str());
  MsgDebug(2,MsgLog::Form("Status of infile = %d",infile.good()));
  getline(infile,tempString);
  getline(infile,tempString);
  getline(infile,tempString);
  bool first = true;
  bool second = true;
  while (infile >> day >> time >> amount >> current >> units >> last) {
    auto loc = time.find(".");
    time.erase(time.begin()+loc,time.end());

    day += " "+time;

    time_t time = ConvertStringToTime(day);
    fgCurrentInfo.emplace(time,current);
    if (first || second) {
      MsgDebug(1,MsgLog::Form("Input: %s, %zu,\n sec %d min %d hour %d mday %d mon %d year %d isdst %d",
            day.c_str(),time,
            fgLocalTM.tm_sec,
            fgLocalTM.tm_min,
            fgLocalTM.tm_hour,
            fgLocalTM.tm_mday,
            fgLocalTM.tm_mon,
            fgLocalTM.tm_year,
            fgLocalTM.tm_isdst));
      if (first) {
        first = false;
      } else {
        second = false;
      }
    }
  }
  infile.close();

  Reset();
}

/*!************************************************************************************************
 * \fn double CCMBeamInfo::GetCurrent(std::string time,std::string format)
 * \brief Return the average current down to the second recorded by the accelerator
 * \param[in] time The time as a string object
 * \param[in] format The format the time is in
 * \return The current in microAmps
 **************************************************************************************************/
double CCMBeamInfo::GetCurrent(std::string time, std::string format)
{
  std::time_t localTime = ConvertStringToTime(time,format,true);
  MsgDebug(2,MsgLog::Form("Input: %s, %zu,\n sec %d min %d hour %d mday %d mon %d year %d isdst %d",
        time.c_str(),localTime,
        fgLocalTM.tm_sec,
        fgLocalTM.tm_min,
        fgLocalTM.tm_hour,
        fgLocalTM.tm_mday,
        fgLocalTM.tm_mon,
        fgLocalTM.tm_year,
        fgLocalTM.tm_isdst));
  return GetCurrent(localTime);
}

/*!************************************************************************************************
 * \fn double CCMBeamInfo::GetCurrent(std::time_t time)
 * \brief Return the average current down to the second recorded by the accelerator
 * \param[in] time The time_t object containing the time of the event
 * \return The current in microAmps
 **************************************************************************************************/
double CCMBeamInfo::GetCurrent(std::time_t time)
{
  MsgDebug(2,MsgLog::Form("Looking for time %zu",time));
  auto iter = fgCurrentInfo.find(time);
  if (iter == fgCurrentInfo.end()) {
    return -1;
  }

  return iter->second;
}


/*!************************************************************************************************
 * \fn void CCMBeamInfo::PrintTable(std::string options)
 * \brief Prints the beam info table
 * \param[in] options If options is "short" than only the first 5 lines from the table are printed
 **************************************************************************************************/
void CCMBeamInfo::PrintTable(std::string options)
{
  MsgInfo(MsgLog::Form("Number of seconds loaded = %zu",fgCurrentInfo.size()));
  size_t count = 0;
  bool stopEarly = options.find("short") != std::string::npos;
  for (auto & entry : fgCurrentInfo) {
    MsgInfo(MsgLog::Form("%zu: %g uAmps",entry.first,entry.second));

    ++count;
    if (stopEarly && count == 10) {
      break;
    }
  }
}

/*!************************************************************************************************
 * \fn std::time_t CCMBeamInfo::ConvertStringToTime(std::string time, std::string format, bool print)
 * \brief Converts the time in the format %d-%b-%Y %H:%M:%S
 * \param[in] time The string to parse
 * \param[in] format The format the time is in default is given in the details of the function
 * \param[in] print Flag to print the parsed string (default = false)
 * \return The time as a time_t object.
 *
 * Converts the time to a time_t object (default format is %d-%b-%Y %X)
 * - %d = month day (1-31)
 * - %b = Month short name (e.g. Jan, Feb, May, Apr, etc.)
 * - %Y = Year in 4-digit format (e.g. 2020)
 * - %H = hour in 2 digit format (e.g. 01, 03, 15)
 * - %M = minute in 2 digit format (e.g. 01, 03, 15)
 * - %S = second in 2 digit format (e.g. 01, 03, 15)
 **************************************************************************************************/
std::time_t CCMBeamInfo::ConvertStringToTime(std::string time, std::string format, bool print)
{
  std::stringstream ss(time);
  ss >> std::get_time(&fgLocalTM,format.c_str());
  fgLocalTM.tm_isdst = -1; // to force mktime to determine if the time of day was during DST

  return std::mktime(&fgLocalTM);

}

/*!************************************************************************************************
 * \fn std::time_t CCMBeamInfo::ConvertToTime(int year, int mon, int day, int hour, int min, int sec)
 * \brief Converts the already parsed time into time_t object
 * \param[in] year Number of years since 1900
 * \param[in] mon The month number [0-11]
 * \param[in] day The day int the month [1-31]
 * \param[in] hour The hour in the day [0-23]
 * \param[in] min The minute in the hour [0-59]
 * \param[in] sec The second in the minute [0-59]
 * \param[in] print Flag to print the parsed string (default = false)
 * \return The time as a time_t object.
 **************************************************************************************************/
std::time_t CCMBeamInfo::ConvertToTime(int year, int mon, int day, int hour, int min, int sec, bool print)
{
  if (print) {
    MsgDebug(2,MsgLog::Form("previous %02d-%02d-%04d %02d:%02d:%02d",
          fgLocalTM.tm_mday,fgLocalTM.tm_mon,fgLocalTM.tm_year,fgLocalTM.tm_hour,fgLocalTM.tm_min,fgLocalTM.tm_sec));
  }

  fgLocalTM.tm_year = year;
  fgLocalTM.tm_mday = day;
  fgLocalTM.tm_mon = mon;
  fgLocalTM.tm_hour = hour;
  fgLocalTM.tm_min = min;
  fgLocalTM.tm_sec = sec;
  fgLocalTM.tm_isdst = -1;

  if (print) {
    MsgDebug(2,MsgLog::Form("%02d-%02d-%04d %02d:%02d:%02d",day,mon,year,hour,min,sec));
    MsgDebug(2,MsgLog::Form("%02d-%02d-%04d %02d:%02d:%02d",
          fgLocalTM.tm_mday,fgLocalTM.tm_mon,fgLocalTM.tm_year,fgLocalTM.tm_hour,fgLocalTM.tm_min,fgLocalTM.tm_sec));
  }

  return std::mktime(&fgLocalTM);
}

/*!************************************************************************************************
 * \fn std::time_t CCMBeamInfo::ConvertTimeSinceEPOCHtoTime(unsigned long int time)
 * \brief Converts the number of seconds since EPOCH to time_t object
 * \param[in] time The number of seconds since EPOCH
 * \return The time as a time_t object.
 **************************************************************************************************/
std::time_t CCMBeamInfo::ConvertTimeSinceEPOCHtoTime(unsigned long int time)
{
  return static_cast<std::time_t>(time);
}

/*!************************************************************************************************
 * \fn void CCMBeamInfo::Reset()
 * \brief Reset iterators
 **************************************************************************************************/
void CCMBeamInfo::Reset()
{
  fgCurrentInfoIterBegin = fgCurrentInfo.begin();
  fgCurrentInfoIterEnd = fgCurrentInfo.end();
  fgCurrentInfoIter = fgCurrentInfoIterBegin;
}

/*!************************************************************************************************
 * \fn int CCMBeamInfo::Next()
 * \brief Advance the iterator by one
 * \return 0 if already at the end of the map, 1 otherwise
 **************************************************************************************************/
int CCMBeamInfo::Next()
{
  if (fgCurrentInfoIter == fgCurrentInfoIterEnd) {
    return 0;
  }

  std::advance(fgCurrentInfoIter,1);
  return 1;
}

/*!************************************************************************************************
 * \fn int CCMBeamInfo::Previous()
 * \brief Get the previous iterator
 * \return 0 if already at the end of the map, 1 otherwise
 **************************************************************************************************/
int CCMBeamInfo::Previous()
{
  if (fgCurrentInfoIter == fgCurrentInfoIterBegin) {
    return 0;
  }

  fgCurrentInfoIter = std::prev(fgCurrentInfoIter,1);
  return 1;
}

/*!************************************************************************************************
 * \fn void CCMBeamInfo::Find(std::string time, std::string format)
 * \brief Move the iterator to the position right before the time passed (see #CCMBeamInfo::GetCurrent)
 * \param[in] time The time to searh for
 * \param[in] format The format the time is in
 **************************************************************************************************/
void CCMBeamInfo::Find(std::string time, std::string format)
{
  Find(ConvertStringToTime(time,format));
}

/*!************************************************************************************************
 * \fn void CCMBeamInfo::Find(std::time_t time)
 * \brief Move the iterator to the position right before the time passed (see #CCMBeamInfo::GetCurrent)
 * \param[in] time The time to searh for
 **************************************************************************************************/
void CCMBeamInfo::Find(std::time_t time)
{
  fgCurrentInfoIter = fgCurrentInfo.lower_bound(time);
  Previous();
}

/*!************************************************************************************************
 * \fn std::pair<std::time_t,double> CCMBeamInfo::GetBeamInfo()
 * \brief Get the beam info for the current iterator position
 * \return The time,current pair the iterator is pointing to
 **************************************************************************************************/
std::pair<std::time_t, double> CCMBeamInfo::GetBeamInfo()
{
  return *fgCurrentInfoIter;
}

