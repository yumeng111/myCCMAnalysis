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

/*!**********************************************
 * \brief Get the name of #CCMEventFinderID_t object in string form
 * \param[in] evtFinder The #CCMEventFinderID_t object
 * \return The corresponding name as a string
 ***********************************************/
std::string Utility::ConvertCCMEventFinderIDToString(CCMEventFinderID_t evtFinder)
{
  switch (evtFinder) {
    case kCCMDynamicLengthEventID: return "CCMDynamicLengthEvent";
    case kCCMFixedLengthEventID: return "CCMFixedLengthEvent";
    default: return "NONE";
  }

  return "NONE";
}

/*!**********************************************
 * \brief Converts a string to #CCMEventFinderID_t type
 * \param[in] name The name corresponding to the #CCMEventFinderID_t type
 * \return The corresponding #CCMEventFinderID_t type
 ***********************************************/
CCMEventFinderID_t Utility::ConvertStringToCCMEventFinderID(std::string name)
{
  if (name.compare("CCMDynamicLengthEvent") == 0) {
    return kCCMDynamicLengthEventID;
  }

  if (name.compare("CCMFixedLengthEvent") == 0) {
    return kCCMFixedLengthEventID;
  }

  MsgWarning(MsgLog::Form("%s does not match current CCMEventFinderID_t types returing kCCMNoneEvent",name.c_str()));

  return kCCMNoneEventID;
}

/*!**********************************************
 * \brief Get the name of #CCMAccumWaveformMethod_t object in string form
 * \param[in] accumWaveform The #CCMAccumWaveformMethod_t object
 * \return The corresponding name as a string
 ***********************************************/
std::string Utility::ConvertCCMAccumWaveformMethodToString(CCMAccumWaveformMethod_t accumWaveform)
{
  switch (accumWaveform) {
    case kCCMAccumWaveformTriangleID: return "CCMAccumWaveformTriangle";
    case kCCMAccumWaveformStartID: return "CCMAccumWaveformStart";
    case kCCMAccumWaveformTrianglePulseCutID: return "CCMAccumWaveformTrianglePulseCut";
    case kCCMAccumWaveformStartPulseCutID: return "CCMAccumWaveformStartPulseCut";
    case kCCMAccumWaveformTotalID: return "CCMAccumWaveformTotal";
    default: return "NONE";
  }

  return "NONE";
}

/*!**********************************************
 * \brief Converts a string to #CCMAccumWaveformMethod_t type
 * \param[in] name The name corresponding to the #CCMAccumWaveformMethod_t type
 * \return The corresponding #CCMAccumWaveformMethod_t type
 ***********************************************/
CCMAccumWaveformMethod_t Utility::ConvertStringToCCMAccumWaveformMethod(std::string name)
{
  if (name.compare("CCMAccumWaveformTriangle") == 0) {
    return kCCMAccumWaveformTriangleID;
  }

  if (name.compare("CCMAccumWaveformStart") == 0) {
    return kCCMAccumWaveformStartID;
  }

  if (name.compare("CCMAccumWaveformTrianglePulseCut")  == 0) {
    return kCCMAccumWaveformTrianglePulseCutID;
  }

  if (name.compare("CCMAccumWaveformStartPulseCut") == 0) {
    return kCCMAccumWaveformStartPulseCutID;
  }

  MsgWarning(MsgLog::Form("%s does not match current CCMAccumWaveformMethod_t types returing kCCMAccumWaveformTotalID",name.c_str()));

  return kCCMAccumWaveformTotalID;
}

/*!**********************************************
 *  \brief Shift the time of the event based on the BCM and FP3 time offsets
 *  \param[in] start The start bin of the pulse
 *  \param[in] beamOffset The time in the DAQ window the BCM was observed
 *  \return The time of the event in ns (input is bin count) 
 *
 *  Shift the time of the pulse to account for the jitter of the BCM. #Utility::fgkFP3Offset 
 *  is used to account for the time difference between the PMTs in CCM and the 
 *  EJ301 detector in FP3 
 *
 *  If beamTime == 0 then the trigger was STROBE or LED so shift the time of the event 
 *  based on the DAQ window true start time (#Utility::fgkWindowStartTime)
 ***********************************************/
double Utility::ShiftTime(int time, int beamOffset)
{
  // shift the time of the pulse to account for the jitter of the BCM
  // the shift by Utility::fgkFP3Offset is to account for the time difference
  // between the PMTs in CCM and the EJ301 detector in FP3
  double shiftTime = static_cast<double>(time - beamOffset + fgkFP3Offset)*fgkBinWidth;

  // if fBeamTime == 0 then the trigger was STROBE or LED
  // so shift the time of the event based on the DAQ window
  // true start time #fgkWindowStartTime
  if (beamOffset == 0) {
    shiftTime += fgkWindowStartTime;
  }

  return shiftTime;
}

