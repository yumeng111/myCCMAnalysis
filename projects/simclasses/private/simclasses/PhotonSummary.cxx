#include "simclasses/PhotonSummary.h"
#include <icetray/serialization.h>
#include <ostream>

std::ostream& PhotonSummary::Print(std::ostream& os) const{
    os << "[ PhotonSummary::"
        << "\n  G4 Distance UV  :" << g4_distance_uv 
        << "\n  G4 Distance Visible :" << g4_distance_visible 
        << "\n  Calculated Distance UV  :" << calculated_distance_uv 
        << "\n  Calculated Distance Visible :" << calculated_distance_visible 
        << "\n  G4 Global Time :" << g4_time
        << "\n  Calculated Time :" << calculated_time
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

