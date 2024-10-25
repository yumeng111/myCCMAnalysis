#ifndef WLSLocation_H_INCLUDED
#define WLSLocation_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include "dataclasses/I3Vector.h"

static const unsigned wlslocation_version_ = 0;

class WLSLocation : public I3FrameObject {
    public:
    enum class WLSLoc : int8_t {
        Unknown = 0,
        PMT = 1,
        FoilTop = 2,
        FoilBottom = 3,
        FoilSides = 4
    };
    
    static const std::unordered_map<WLSLocation::WLSLoc, std::string> wlsLocationToProcessName;

    WLSLoc wls_loc; 

    SET_LOGGER("WLSLocation");

    bool operator==(const WLSLocation& rhs) const {
        return wls_loc == rhs.wls_loc;
    }


    WLSLocation(WLSLoc wls_loc_) : wls_loc(wls_loc_) {}
    WLSLocation() : wls_loc(WLSLocation::WLSLoc::Unknown) {}

    std::ostream& Print(std::ostream&) const;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

};

I3_CLASS_VERSION(WLSLocation,wlslocation_version_);

typedef I3Vector<WLSLocation> WLSLocationSeries;

std::ostream& operator<<(std::ostream&, const WLSLocation&);
I3_POINTER_TYPEDEFS(WLSLocation);
I3_POINTER_TYPEDEFS(WLSLocationSeries);
I3_DEFAULT_NAME(WLSLocation);
#endif
