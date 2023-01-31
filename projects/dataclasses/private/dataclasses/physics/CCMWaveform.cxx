#include <icetray/serialization.h>
#include <dataclasses/physics/CCMWaveform.h>

#include <algorithm>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>


using namespace std;
using namespace boost::lambda;


CCMWaveform::StatusCompound::~StatusCompound() {}

template <class Archive>
void CCMWaveform::StatusCompound::save(Archive& ar, unsigned version) const
{
  ar & make_nvp("interval", interval_);
  ar & make_nvp("status", status_);
}

template <class Archive>
void CCMWaveform::StatusCompound::load(Archive& ar, unsigned version)
{
  if (version>ccmwaveform_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMWaveform class",version,ccmwaveform_version_);

  ar & make_nvp("interval", interval_);
  ar & make_nvp("status", status_);
}

I3_SPLIT_SERIALIZABLE(CCMWaveform::StatusCompound);

unsigned CCMWaveform::GetStatus(const vector<StatusCompound>& waveformInfo)
{
  unsigned retVal = VIRGINAL;

  std::vector<StatusCompound>::const_iterator it = waveformInfo.begin();
  for ( ; it != waveformInfo.end(); it++)
    retVal |= it->GetStatus();

  return retVal;
}

unsigned
CCMWaveform::GetStatus() const
{
	return GetStatus(waveformInfo_);
}

CCMWaveform::~CCMWaveform() {}


template <class Archive>
void CCMWaveform::save(Archive& ar, unsigned version) const
{
  if (version>ccmwaveform_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMWaveform class.",version,ccmwaveform_version_);

  ar & make_nvp("startTime", startTime_);
  ar & make_nvp("binWidth", binWidth_);
  ar & make_nvp("waveform", waveform_);
  ar & make_nvp("waveformInformation", waveformInfo_);
  ar & make_nvp("source", source_);
}

template <class Archive>
void CCMWaveform::load(Archive& ar, unsigned version)
{
  if (version>ccmwaveform_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMWaveform class.",version,ccmwaveform_version_);

  ar & make_nvp("startTime", startTime_);
  ar & make_nvp("binWidth", binWidth_);
  ar & make_nvp("waveform", waveform_);
  ar & make_nvp("waveformInformation", waveformInfo_);

  ar & make_nvp("source", source_);
}

bool
operator==(const CCMWaveform& lhs, const CCMWaveform& rhs)
{
  return lhs.GetStartTime() == rhs.GetStartTime() 
    && lhs.GetBinWidth() == rhs.GetBinWidth() 
    && lhs.GetWaveform() == rhs.GetWaveform()
    && lhs.GetWaveformInformation() == rhs.GetWaveformInformation();
}

std::ostream& CCMWaveform::Print(std::ostream& oss) const
{
  std::string srcstr;
  if (GetSource() == CCMWaveform::V1730) srcstr.append("V1730");
  oss << "[CCMWaveform:\n"
      << "  StartTime : " << GetStartTime() << '\n'
      << "     Source : " << srcstr << "\n]";
  return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMWaveform& wf)
{
  return(wf.Print(oss));
}

std::ostream& CCMWaveform::StatusCompound::Print(std::ostream& oss) const
{
  std::string srcstr;
  if (GetStatus() == CCMWaveform::VIRGINAL) srcstr.append("VIRGINAL");
  if (GetStatus() == CCMWaveform::COMBINED) srcstr.append("COMBINED");
  if (GetStatus() == CCMWaveform::SATURATED) srcstr.append("SATURATED");
  if (GetStatus() == CCMWaveform::UNDERSHOT) srcstr.append("UNDERSHOT");

  oss << "[CCMWaveform::StatusCompound: \n"
      << "                        Range : " << GetInterval().first << "--" << GetInterval().second << '\n'
      << "                       Status : " << srcstr << "\n]";
  return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMWaveform::StatusCompound& sc)
{
  return(sc.Print(oss));
}

I3_SPLIT_SERIALIZABLE(CCMWaveform);

I3_SERIALIZABLE(CCMWaveformSeriesMap);
