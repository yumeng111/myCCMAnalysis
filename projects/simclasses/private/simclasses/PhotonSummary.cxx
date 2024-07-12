#include "simclasses/PhotonSummary.h"
#include <icetray/serialization.h>
#include <ostream>

std::ostream& PhotonSummary::Print(std::ostream& os) const{
    os << "[ PhotonSummary::"
        << "\n  Distance UV  :" << distance_uv 
        << "\n  Distance Visible :" << distance_visible 
        << "\n  Number WLS :" << n_wls 
        << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& oss, PhotonSummary const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, PhotonSummary & bcm) {
    return(bcm.Print(oss));
}

I3_SERIALIZABLE(PhotonSummary);
I3_SERIALIZABLE(PhotonSummarySeries);

