
#ifndef CCMOMGEO_H_INCLUDED
#define CCMOMGEO_H_INCLUDED

#include "dataclasses/I3Position.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/I3Orientation.h"
#include "dataclasses/I3Map.h"
#include <string>
#include <iostream>
#include <sstream>
#include "icetray/OMKey.h"
#include "dataclasses/Utility.h"


static const unsigned ccmomgeo_version_ = 0;


/**
 * List the names of enumeration members defined in this file
 * here. These can be used for e.g. pybindings, which require
 * the names of the enumeration members to be known. This list
 * should be updated whenver members or new enums are added to
 * the class.
 */
#define CCMOMGEO_H_CCMOMGeo_OMType        \
  (UnknownType)(CCM8inCoated)(CCM8inUncoated)(CCM1in)

//Simple struct to contain all pertinent OM info.
//See CCMGeometry.h for more info

class CCMOMGeo
{
public:
    enum OMType {UnknownType = 0, CCM8inCoated = 10, CCM8inUncoated = 20, CCM1in = 30,};

    CCMOMGeo():omtype(UnknownType){}

    ~CCMOMGeo();

    /**
     * the OM's (or PMT's) x,y,z position
     */
    I3Position position;

    /**
     * Orientation of the OM (or PMT)
     */
    I3Orientation orientation;

    /**
     * InIce? IceTop? AMANDA OM? (see enum above)
     */
    OMType omtype;

    /**
     * Effective collection area (use I3Units)
     */
    double area;

    /**
     * Gets the I3Direction from the I3Orientation
     */
    inline I3Direction GetDirection() const {return orientation.GetDir();}

    bool operator==(const CCMOMGeo& rhs) const
    {
      return (position == rhs.position &&
              orientation == rhs.orientation &&
              omtype == rhs.omtype);
    }
    bool operator!=(const CCMOMGeo& rhs) const
    {
      return !operator==(rhs);
    }

    std::ostream& Print(std::ostream&) const;

private:
    friend class icecube::serialization::access;
    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

std::ostream& operator<<(std::ostream&, const CCMOMGeo&);

I3_POINTER_TYPEDEFS(CCMOMGeo);
I3_CLASS_VERSION(CCMOMGeo, ccmomgeo_version_);

typedef I3Map<OMKey, CCMOMGeo> CCMOMGeoMap;
I3_POINTER_TYPEDEFS(CCMOMGeoMap);

#endif //CCMOMGEO_H_INCLUDED
