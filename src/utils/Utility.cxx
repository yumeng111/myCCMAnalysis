/*!**********************************************
 * \file Utility.cxx
 * \brief Source code for the #Utility class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "Utility.h"
#include "MsgLog.h"

/*!**********************************************
 * \fn void Utility::ParseStringForRunNumber(std::string name, int & run int & subrun)
 * \brief Parses the passed string to get the run and subrun numbers
 * \param[in] name The string to parse
 * \param[out] run The output run number
 * \param[out] subrun The output subrun number
 *
 * Assumes the format of the name is the same as to which was used
 * during the 2019 run cycle, which was
 * PDSout_run<RUNNUMBER>-<SUBRUNNUMBER>_<DATE>.bin(or root)
 *
 * where RUNNUMBER and SUBRUNNUMBER were 6 digit numbers
 * and DATE was in the format of %m-%d-%H%M. This function
 * currently does nothing with the DATE
 ***********************************************/
void Utility::ParseStringForRunNumber(std::string name, int & run, int & subrun)
{
  MsgDebug(2,MsgLog::Form("Input string: %s",name.c_str()));
  // find where "run" occurs in the string
  size_t locOfRun = name.find("PDSout_run");
  if (locOfRun == std::string::npos) {
    MsgDebug(2,"run number cannot be found from file name");
    run = 0;
    subrun = 0;
    return;
  }
  MsgDebug(2,MsgLog::Form("Location of run = %zu", locOfRun));

  // get the run number and the subrun number in string format
  std::string runNumber = "";
  std::string subRunNumber = "";
  runNumber.assign(name,locOfRun+3+7,6);
  subRunNumber.assign(name,locOfRun+3+7+6+1,6);

  MsgDebug(2,MsgLog::Form("String form runNumber = %s, subRunNumber = %s",runNumber.c_str(),subRunNumber.c_str()));

  // convert the strings to int
  run = std::stoi(runNumber);
  subrun = std::stoi(subRunNumber);

  MsgDebug(2,MsgLog::Form("Integer form runNumber %d subRunNumber = %d",run,subrun));

  return;
}

