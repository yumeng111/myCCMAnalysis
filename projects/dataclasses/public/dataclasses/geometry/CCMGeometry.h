
#ifndef CCMGEOMETRY_H_INCLUDED
#define CCMGEOMETRY_H_INCLUDED

#include <vector>
#include "dataclasses/Utility.h"
#include "dataclasses/geometry/CCMOMGeo.h"
#include "dataclasses/I3Time.h"
#include <icetray/I3DefaultName.h>
#include "icetray/I3FrameObject.h"
#include "icetray/CCMTriggerKey.h"
#include "dataclasses/I3Map.h"

static const unsigned ccmgeometry_version_ = 0;

class CCMGeometry : public I3FrameObject {
public:
  CCMGeometry(){};
  ~CCMGeometry();

  //Map of all OMs based on their OMKey
  I3Map<CCMPMTKey, CCMOMGeo> pmt_geo;
  I3Map<CCMPMTKey, uint32_t> pmt_channel_map;
  I3Map<CCMTriggerKey, uint32_t> trigger_channel_map;
  I3Map<CCMPMTKey, CCMTriggerKey> trigger_copy_map;

  I3Time startTime;
  I3Time endTime;

  const CCMGeometry& operator=(const CCMGeometry& geometry);

  bool operator==(const CCMGeometry& rhs)
  {
    return (pmt_geo == rhs.pmt_geo &&
            pmt_channel_map == rhs.pmt_channel_map &&
            trigger_channel_map == rhs.trigger_channel_map &&
            trigger_copy_map == rhs.trigger_copy_map &&
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
