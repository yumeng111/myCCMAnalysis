#include <simclasses/CCMMCPE.h>
#include <ostream>

I3_SERIALIZABLE(CCMMCPESeriesMap);
I3_SERIALIZABLE(CCMMCPESeries);

const std::unordered_map<CCMMCPE::PhotonSource, std::string> CCMMCPE::photonSourceToProcessName = {{CCMMCPE::PhotonSource::Unknown, "Unknown"},
                                                                                                      {CCMMCPE::PhotonSource::Scintillation, "Scintillation"},
                                                                                                      {CCMMCPE::PhotonSource::Cerenkov, "Cerenkov"},
                                                                                                      {CCMMCPE::PhotonSource::OpWLS, "OpWLS"}};
std::ostream& operator<<(std::ostream& os, const CCMMCPE& pe) {
    return(pe.Print(os));
}

std::ostream& CCMMCPE::Print(std::ostream& os) const{
    os << "[ CCMMCPE::"
        << "\n  ParentID :" << parent_id 
        << "\n  TrackID :" << track_id
        << "\n  Global Time :" << global_time
        << "\n  Local Time :" << local_time
        << "\n  Wavelength :" << wavelength
        << "\n  Distance UV :" << distance_uv
        << "\n  Distance Visible :" << distance_visible
        << "\n  Position :" << position 
        << "\n  Direction :" << direction 
        << "\n PhotonSource :" << photonSourceToProcessName.at(photon_source) 
        //std::ostream & operator<<(std::ostream & os, CCMMCPE::PhotonSource const & p) {
        //    if (CCMMCPE::PhotonSourceNames.find(p) != CCMMCPE::PhotonSourceNames.end())
        //        << "\n  PhotonSource :" << CCMMCPE::PhotonSourceNames.at(p);
        << " ]";
    return os;
}


