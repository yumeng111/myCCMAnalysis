#include <icetray/serialization.h>
#include <dataclasses/geometry/CCMGeometry.h>

CCMGeometry::~CCMGeometry() {}

template <class Archive>
void
CCMGeometry::serialize(Archive& ar, unsigned version) {
  if (version > ccmgeometry_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMGeometry class.", version, ccmgeometry_version_);

  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("OMGeo", omgeo);
  ar & make_nvp("StartTime", startTime);
  ar & make_nvp("EndTime", endTime);
}

const CCMGeometry& CCMGeometry::operator=(const CCMGeometry& geometry) {
  omgeo = geometry.omgeo;
  startTime = geometry.startTime;
  endTime = geometry.endTime;

  return *this;
}

I3_SERIALIZABLE(CCMGeometry);
