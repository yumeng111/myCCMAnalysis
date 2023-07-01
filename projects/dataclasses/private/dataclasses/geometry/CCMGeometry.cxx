#include <icetray/serialization.h>
#include <dataclasses/geometry/CCMGeometry.h>

CCMGeometry::~CCMGeometry() {}

template <class Archive>
void
CCMGeometry::serialize(Archive& ar, unsigned version) {
  if (version > ccmgeometry_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMGeometry class.", version, ccmgeometry_version_);

  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("PMTGeo", pmt_geo);
  ar & make_nvp("PMTChannelFormat", pmt_channel_map);
  ar & make_nvp("TriggerChannelMap", trigger_channel_map);
  ar & make_nvp("TriggerCopyMap", trigger_copy_map);
  ar & make_nvp("TriggerToTriggerCopyMap", trigger_to_trigger_copy_map);
  ar & make_nvp("StartTime", startTime);
  ar & make_nvp("EndTime", endTime);
}

const CCMGeometry& CCMGeometry::operator=(const CCMGeometry& geometry) {
  pmt_geo = geometry.pmt_geo;
  startTime = geometry.startTime;
  endTime = geometry.endTime;

  return *this;
}

I3_SERIALIZABLE(CCMGeometry);
