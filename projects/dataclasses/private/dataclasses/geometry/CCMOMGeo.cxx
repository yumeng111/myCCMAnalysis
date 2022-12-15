#include <icetray/serialization.h>
#include <dataclasses/geometry/CCMOMGeo.h>

CCMOMGeo::~CCMOMGeo() { }

template <class Archive>
void
CCMOMGeo::serialize(Archive& ar, unsigned version)
{
    if (version>ccmomgeo_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMOMGeo class.",version,ccmomgeo_version_);

    ar & make_nvp("Position", position);
    ar & make_nvp("Orientation", orientation);
    ar & make_nvp("OMType", omtype);
    ar & make_nvp("Area", area);
}

I3_SERIALIZABLE(CCMOMGeo);
I3_SERIALIZABLE(CCMOMGeoMap);

std::ostream& CCMOMGeo::Print(std::ostream& os) const{
    os << "[CCMOMGeo Position: " << position << '\n'
       << "      Orientation: " << orientation << '\n'
       << "           OMType: " << omtype << '\n'
       << "             Area: " << area << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const CCMOMGeo& g){
    return(g.Print(os));
}
