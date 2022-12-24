
#ifndef CCMGEOMETRY_H_INCLUDED
#define CCMGEOMETRY_H_INCLUDED

#include <vector>
#include "dataclasses/Utility.h"
#include "dataclasses/geometry/CCMOMGeo.h"
#include "dataclasses/I3Time.h"
#include <icetray/I3DefaultName.h>
#include "icetray/I3FrameObject.h"
#include "dataclasses/I3Map.h"

static const unsigned ccmgeometry_version_ = 0;

class CCMGeometry : public I3FrameObject {
public:
  CCMGeometry(){};
  ~CCMGeometry();

  //Map of all OMs based on their OMKey
  I3Map<OMKey, CCMOMGeo> om_map;

  I3Time startTime;
  I3Time endTime;

  const CCMGeometry& operator=(const CCMGeometry& geometry);

  bool operator==(const CCMGeometry& rhs)
  {
    return (omgeo == rhs.omgeo &&
            startTime == rhs.startTime &&
            endTime == rhs.endTime);
  }
  bool operator!=(const CCMGeometry& rhs)
  {
    return !operator==(rhs);
  }

  friend class icecube::serialization::access;

  template <class Archive> void serialize(Archive & ar, unsigned version);
};

I3_DEFAULT_NAME(CCMGeometry);
I3_POINTER_TYPEDEFS(CCMGeometry);
I3_CLASS_VERSION(CCMGeometry, ccmgeometry_version_);

#endif // CCMGEOMETRY_H_INCLUDED
