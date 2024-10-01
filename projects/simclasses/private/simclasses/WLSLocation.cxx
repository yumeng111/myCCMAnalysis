#include "simclasses/WLSLocation.h"
#include <icetray/serialization.h>
#include <ostream>

const std::unordered_map<WLSLocation::WLSLoc, std::string> WLSLocation::wlsLocationToProcessName = {{WLSLocation::WLSLoc::Unknown, "Unknown"},
                                                                                                    {WLSLocation::WLSLoc::PMT, "PMT"},
                                                                                                    {WLSLocation::WLSLoc::FoilTop, "FoilTop"},
                                                                                                    {WLSLocation::WLSLoc::FoilBottom, "FoilBottom"},
                                                                                                    {WLSLocation::WLSLoc::FoilSides, "FoilSides"}};



template <class Archive>
void WLSLocation::save(Archive& ar, unsigned version) const {
    if (version > wlslocation_version_)
        log_fatal("Attempting to save version %u from file but running version %u of WLSLocation class.", version, wlslocation_version_);
    
    ar & make_nvp("wls_loc",wls_loc);
}

template <class Archive>
void WLSLocation::load(Archive& ar, unsigned version) {
    if (version > wlslocation_version_)
        log_fatal("Attempting to read version %u from file but running version %u of WLSLocation class.", version, wlslocation_version_);

    ar & make_nvp("wls_loc",wls_loc);
}


std::ostream& WLSLocation::Print(std::ostream& os) const{
    os << "[ WLSLocation::"
        << "\n  WLS Location :" << wlsLocationToProcessName.at(wls_loc) 
        << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& oss, WLSLocation const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, WLSLocation & bcm) {
    return(bcm.Print(oss));
}

I3_SPLIT_SERIALIZABLE(WLSLocation);
I3_SERIALIZABLE(WLSLocationSeries);

