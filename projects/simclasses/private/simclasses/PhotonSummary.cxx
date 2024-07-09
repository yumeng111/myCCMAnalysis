#include "simclasses/PhotonSummary.h"
#include <icetray/serialization.h>
#include <ostream>

const std::unordered_map<PhotonSummary::CreationProcess, std::string> PhotonSummary::CreationProcesstoName = {{PhotonSummary::CreationProcess::Unknown, "Unknown"},
                                                                                                      {PhotonSummary::CreationProcess::Scintillation, "Scintillation"},
                                                                                                      {PhotonSummary::CreationProcess::Cerenkov, "Cerenkov"},
                                                                                                      {PhotonSummary::CreationProcess::OpWLS, "OpWLS"}};

std::ostream& PhotonSummary::Print(std::ostream& os) const{
    os << "[ PhotonSummary::"
        << "\n  ParentID :" << parent_id 
        << "\n  TrackID :" << track_id
        << "\n  Time :" << time 
        << "\n  Position :" << position
        << "\n  Creation Process :" << CreationProcesstoName.at(creation_process) 
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

