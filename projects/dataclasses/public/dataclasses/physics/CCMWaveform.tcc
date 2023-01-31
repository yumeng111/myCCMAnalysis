#include <icetray/serialization.h>
#include <dataclasses/physics/CCMWaveform.h>

#include <algorithm>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/foreach.hpp>

template<typename T>
unsigned CCMWaveform<T>::GetStatus(const std::vector<CCMWaveform<T>::StatusCompound>& waveformInfo) {
  unsigned retVal = VIRGINAL;

  typename std::vector<CCMWaveform<T>::StatusCompound>::const_iterator it = waveformInfo.begin();
  for ( ; it != waveformInfo.end(); it++)
    retVal |= it->GetStatus();

  return retVal;
}

template<typename T>
unsigned CCMWaveform<T>::GetStatus() const {
	return GetStatus(waveformInfo_);
}

template<typename T>
CCMWaveform<T>::~CCMWaveform() {}

template <typename T>
template <class Archive>
void CCMWaveform<T>::save(Archive& ar, unsigned version) const {
  if (version>ccmwaveform_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMWaveform class.",version,ccmwaveform_version_);

  ar & make_nvp("startTime", startTime_);
  ar & make_nvp("binWidth", binWidth_);
  ar & make_nvp("waveform", waveform_);
  ar & make_nvp("waveformInformation", waveformInfo_);
  ar & make_nvp("source", source_);
}

template <typename T>
template <class Archive>
void CCMWaveform<T>::load(Archive& ar, unsigned version) {
  if (version>ccmwaveform_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMWaveform class.",version,ccmwaveform_version_);

  ar & make_nvp("startTime", startTime_);
  ar & make_nvp("binWidth", binWidth_);
  ar & make_nvp("waveform", waveform_);
  ar & make_nvp("waveformInformation", waveformInfo_);

  ar & make_nvp("source", source_);
}

template<typename T>
bool operator==(const CCMWaveform<T>& lhs, const CCMWaveform<T>& rhs) {
  return lhs.GetStartTime() == rhs.GetStartTime()
    && lhs.GetBinWidth() == rhs.GetBinWidth()
    && lhs.GetWaveform() == rhs.GetWaveform()
    && lhs.GetWaveformInformation() == rhs.GetWaveformInformation();
}

template<typename T>
std::ostream& CCMWaveform<T>::Print(std::ostream& oss) const {
  std::string srcstr;
  if (GetSource() == CCMSource::V1730) srcstr.append("V1730");
  oss << "[CCMWaveform:\n"
      << "  StartTime : " << GetStartTime() << '\n'
      << "     Source : " << srcstr << "\n]";
  return oss;
}

template<typename T>
std::ostream& operator<<(std::ostream& oss, const CCMWaveform<T>& wf) {
  return(wf.Print(oss));
}

