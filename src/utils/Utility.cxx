/*!**********************************************
 * \file Utility.cxx
 * \brief Source code for the #Utility class
 * \author R. T. Thornton (LANL)
 * \date February 25, 2020
 ***********************************************/
#include "Utility.h"
#include "MsgLog.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>

/*!**********************************************
 * \fn void Utility::ParseStringForRunNumber(std::string name, int & run int & subrun, int * month, int * day, int * time)
 * \brief Parses the passed string to get the run and subrun numbers
 * \param[in] name The string to parse
 * \param[out] run The output run number
 * \param[out] subrun The output subrun number
 * \param[out] month The output month
 * \param[out] day The output day
 * \param[out] time The output time
 *
 * Assumes the format of the name is the same as to which was used
 * during the 2019 run cycle, which was
 * PDSout_run<RUNNUMBER>-<SUBRUNNUMBER>_<DATE>.bin(or root)
 *
 * where RUNNUMBER and SUBRUNNUMBER were 6 digit numbers
 * and DATE was in the format of %m-%d-%H%M. This function
 * currently does nothing with the DATE
 ***********************************************/
void Utility::ParseStringForRunNumber(std::string name, int & run, int & subrun, int * month, int * day, int * time)
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
  std::string dayS = "";
  std::string monthS= "";
  std::string timeS= "";
  size_t startLocation = locOfRun+3+7;
  runNumber.assign(name,startLocation,6);
  startLocation += 7;
  subRunNumber.assign(name,startLocation,6);
  startLocation += 7;
  monthS.assign(name,startLocation,2);
  startLocation += 3;
  dayS.assign(name,startLocation,2);
  startLocation += 3;
  timeS.assign(name,startLocation,4);

  if (MsgLog::GetGlobalDebugLevel() >= 2) {
    MsgDebug(2,MsgLog::Form("String form runNumber = %s, subRunNumber = %s month = %s day = %s time = %s",
          runNumber.c_str(),subRunNumber.c_str(),monthS.c_str(),dayS.c_str(),timeS.c_str()));
  }

  // convert the strings to int
  run = std::stoi(runNumber);
  subrun = std::stoi(subRunNumber);

  if (MsgLog::GetGlobalDebugLevel() >= 3) {
    MsgDebug(3,MsgLog::Form("Integer form runNumber %d subRunNumber = %d",run,subrun));
  }

  if (day) {
    *day = std::stoi(dayS);
    if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("Integer form day = %d",*day));
    }
  }
  if (month) {
    *month = std::stoi(monthS);
    if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("Integer form month = %d",*month));
    }
  }
  if (time) {
    *time = std::stoi(timeS);
    if (MsgLog::GetGlobalDebugLevel() >= 3) {
      MsgDebug(3,MsgLog::Form("Integer form time = %d",*time));
    }
  }

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
 *  \param[in] applyFP3Offset If true apply the FP3Offset correction (default: true)
 *  \return The time of the event in ns (input is bin count) 
 *
 *  Shift the time of the pulse to account for the jitter of the BCM. #Utility::fgkFP3Offset 
 *  is used to account for the time difference between the PMTs in CCM and the 
 *  EJ301 detector in FP3 
 *
 *  If beamTime == 0 then the trigger was STROBE or LED so shift the time of the event 
 *  based on the DAQ window true start time (#Utility::fgkWindowStartTime)
 ***********************************************/
double Utility::ShiftTime(int time, int beamOffset, bool applyFP3Offset)
{
  // shift the time of the pulse to account for the jitter of the BCM
  // the shift by Utility::fgkFP3Offset is to account for the time difference
  // between the PMTs in CCM and the EJ301 detector in FP3
  int temp = time - beamOffset;
  if (applyFP3Offset) {
    temp += fgkFP3Offset;
  }
  double shiftTime = static_cast<double>(temp)*fgkBinWidth;

  // if fBeamTime == 0 then the trigger was STROBE or LED
  // so shift the time of the event based on the DAQ window
  // true start time #fgkWindowStartTime
  if (beamOffset == 0) {
    shiftTime += fgkWindowStartTime;
  }

  return shiftTime;
}

/*!**********************************************
 *  \brief Undo the shift time that was applied to calculate the start tiem of the event
 *  \param[in] start The start bin of the pulse
 *  \param[in] beamOffset The time in the DAQ window the BCM was observed
 *  \param[in] applyFP3Offset If true apply the FP3Offset correction (default: true)
 *  \return The time of the event in ns (input is bin count) 
 *
 *  Shift the time of the pulse to account for the jitter of the BCM. #Utility::fgkFP3Offset 
 *  is used to account for the time difference between the PMTs in CCM and the 
 *  EJ301 detector in FP3 
 *
 *  If beamTime == 0 then the trigger was STROBE or LED so shift the time of the event 
 *  based on the DAQ window true start time (#Utility::fgkWindowStartTime)
 ***********************************************/
double Utility::UndoShiftTime(double time, int beamOffset, bool applyFP3Offset)
{
  // shift the time of the pulse to account for the jitter of the BCM
  // the shift by Utility::fgkFP3Offset is to account for the time difference
  // between the PMTs in CCM and the EJ301 detector in FP3
  double temp = time;
  if (beamOffset == 0) {
    temp -= fgkWindowStartTime;
  }
  temp /= fgkBinWidth;

  if (applyFP3Offset) {
    temp -= fgkFP3Offset;
  }

  temp -= beamOffset;

  return temp;
}


/*!**********************************************
 *  \brief Gets the list of files
 *  \param[in] file_regexp The command to pass to ls, can contain wild cards
 *  \return A vector containg all the files that match file_regexp
 ***********************************************/
std::vector<std::string> Utility::GetListOfFiles(const char * file_regexp)
{
  // Cheap way to get list of files matching the expression (assumes
  // some flavor of unix
  char tmp_name[128] = {"/tmp/IOMOD.XXXXXX"};
  int unused = mkstemp(tmp_name);
  std::string temp_file_regexp(file_regexp);
  bool inLustre = temp_file_regexp.find("lustre") != std::string::npos;
  std::string cmd;
  if (!inLustre) {
    cmd  = "(ls ";
  } else {
    cmd = "(lfs find ";
  }
  cmd += file_regexp;
  cmd += " > ";
  cmd += tmp_name;
  cmd += ") >& /dev/null";
  MsgInfo(MsgLog::Form("Command Used %s",cmd.c_str()));
  unused = system(cmd.c_str());

  std::vector<std::string> localVec;
  // Open the temp file and get the list of files which matched
  static const int s = 256;
  char             buff[s];
  FILE*            fp     = fopen(tmp_name,"r");
  std::string      file;
  while ( fgets(buff,s-1,fp) != NULL ) {
    int len = strlen(buff);
    while (buff[len]==' ' || buff[len]=='\n' || buff[len]=='\0') --len;
    buff[len+1]='\0';
    file = buff;
    localVec.emplace_back(file);
  }
  fclose(fp);

  // Remove the temporary file
  std::remove(tmp_name);

  return localVec;
}

/*!**********************************************
 *  \brief Reads in a list of files and puts them in a vector of strings
 *  \param[in] file The file name that contains the list of files
 *  \param[out] infileList The vector containing the list of files
 ***********************************************/
void Utility::IndirectFileList(const char* file, std::vector<std::string>& infileList)
{
  // Open the file and pull the file names in
  std::ifstream infile;
  infile.open(file);
  if (!infile) {
    MsgError(MsgLog::Form("File %s not found.",file));
    exit(1);
  }

  std::string fname;
  std::string line;
  while (infile >> line) {
    // from the indirect file
    fname = line;
    infileList.push_back(fname);
  }
  infile.close();
}

std::istream& std::operator >> (std::istream& is, std::pair<int, double>& ps)
{
  return is >> ps.first >> ps.second;
}
std::ostream& std::operator << (std::ostream& os, const std::pair<const int, double>& ps)
{
  return os << ps.first << "==>>" << ps.second;
}
