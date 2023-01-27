#include <icetray/serialization.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/external/CompareFloatingPoint.h>
#include <string>

using CompareFloatingPoint::Compare;

CCMRecoPulse::~CCMRecoPulse() {}

template <class Archive> 
void 
CCMRecoPulse::serialize(Archive& ar, unsigned version)
{
	if (version>ccmrecopulse_version_)
		log_fatal("Attempting to read version %u from file but running version %u of CCMRecoPulse class.",version,ccmrecopulse_version_);

	ar & make_nvp("Time", time_);
	ar & make_nvp("PulseCharge", charge_);
	ar & make_nvp("Width", width_);
}


bool 
CCMRecoPulse::operator==(const CCMRecoPulse& rhs) const
{
  return (time_ == rhs.time_ || (std::isnan(time_) && std::isnan(rhs.time_)))
      && (charge_ == rhs.charge_ || (std::isnan(charge_) && std::isnan(rhs.charge_)))
      && (width_ == rhs.width_ || (std::isnan(width_) && std::isnan(rhs.width_)));
}

std::ostream& CCMRecoPulse::Print(std::ostream& oss) const{
  oss << "[CCMRecoPulse:\n"
      << "             Time : " << GetTime() << std::endl
      << "           Charge : " << GetCharge() << std::endl
      << "            Width : " << GetWidth()  << "\n]";
  return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMRecoPulse& p){
  return(p.Print(oss));
}


I3_SERIALIZABLE(CCMRecoPulse);

I3_SERIALIZABLE(CCMRecoPulseSeriesMap);
I3_SERIALIZABLE(CCMRecoPulseMap);
