#include <icetray/serialization.h>
#include <dataclasses/physics/CCMWaveform.h>

#include <algorithm>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>

CCMStatusCompound::~CCMStatusCompound() {}

std::ostream& operator<<(std::ostream& oss, CCMStatusCompound const & sc) {
  return(sc.Print(oss));
}

template <class Archive>
void CCMStatusCompound::save(Archive& ar, unsigned version) const {
  ar & make_nvp("interval", interval_);
  ar & make_nvp("status", status_);
}

template <class Archive>
void CCMStatusCompound::load(Archive& ar, unsigned version) {
  if (version>ccmwaveform_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMWaveform class",version,ccmwaveform_version_);

  ar & make_nvp("interval", interval_);
  ar & make_nvp("status", status_);
}

std::ostream& CCMStatusCompound::Print(std::ostream& oss) const {
  std::string srcstr;
  if (GetStatus() == CCMStatus::VIRGINAL) srcstr.append("VIRGINAL");
  if (GetStatus() == CCMStatus::COMBINED) srcstr.append("COMBINED");
  if (GetStatus() == CCMStatus::SATURATED) srcstr.append("SATURATED");
  if (GetStatus() == CCMStatus::UNDERSHOT) srcstr.append("UNDERSHOT");

  oss << "[CCMWaveform::StatusCompound: \n"
      << "                        Range : " << GetInterval().first << "--" << GetInterval().second << '\n'
      << "                       Status : " << srcstr << "\n]";
  return oss;
}

I3_SPLIT_SERIALIZABLE(CCMStatusCompound);
I3_SPLIT_SERIALIZABLE(CCMWaveformUInt16);
I3_SERIALIZABLE(CCMWaveformUInt16Series);
