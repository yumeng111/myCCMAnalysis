#include <map>
#include <icetray/serialization.h>
#include <dataclasses/calibration/CCMCalibration.h>
#include "dataclasses/TankKey.h"

CCMCalibration::CCMCalibration(){}

CCMCalibration::~CCMCalibration(){}

template <class Archive>
void 
CCMCalibration::save(Archive& ar, unsigned version) const
{
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("domcal",domCal);
  ar & make_nvp("StartTime",startTime);
  ar & make_nvp("EndTime",endTime);
}

template <class Archive>
void 
CCMCalibration::load(Archive& ar, unsigned version) {
  if (version>ccmcalibration_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMCalibration class.",
              version,ccmcalibration_version_);
  
  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("domcal",domCal);
  ar & make_nvp("StartTime",startTime);
  ar & make_nvp("EndTime",endTime);
}

I3_SPLIT_SERIALIZABLE(CCMCalibration);

