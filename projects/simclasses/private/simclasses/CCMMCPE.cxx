#include <simclasses/CCMMCPE.h>
#include <ostream>
#include <icetray/serialization.h>

I3_SERIALIZABLE(CCMMCPESeriesMap);
I3_SERIALIZABLE(CCMMCPESeries);

const std::unordered_map<CCMMCPE::PhotonSource, std::string> CCMMCPE::photonSourceToProcessName = {{CCMMCPE::PhotonSource::Unknown, "Unknown"},
                                                                                                      {CCMMCPE::PhotonSource::Scintillation, "Scintillation"},
                                                                                                      {CCMMCPE::PhotonSource::Cerenkov, "Cerenkov"},
                                                                                                      {CCMMCPE::PhotonSource::OpWLS, "OpWLS"}};
template <class Archive>
void CCMMCPE::save(Archive& ar, unsigned version) const {
    if (version > ccmmcpe_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMMCPE class.", version, ccmmcpe_version_);

    ar & make_nvp("parent_id",parent_id);
    ar & make_nvp("track_id",track_id);
    ar & make_nvp("global_time",global_time);
    ar & make_nvp("local_time",local_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("distance_uv",distance_uv);
    ar & make_nvp("distance_visible",distance_visible);
    ar & make_nvp("position",position);
    ar & make_nvp("direction",direction);
    ar & make_nvp("photon_source",photon_source);
}

template <class Archive>
void CCMMCPE::load(Archive& ar, unsigned version) {
    if (version > ccmmcpe_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMMCPE class.", version, ccmmcpe_version_);
    
    // uncomment to load in old version of ccmmcpe
    ar & make_nvp("parent_id",parent_id);
    ar & make_nvp("track_id",track_id);
    ar & make_nvp("time",global_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("position",position);
    ar & make_nvp("direction",direction);
    ar & make_nvp("photon_source",photon_source);

    // uncomment to load in new version of ccmmcpe 
    //ar & make_nvp("parent_id",parent_id);
    //ar & make_nvp("track_id",track_id);
    //ar & make_nvp("global_time",global_time);
    //ar & make_nvp("local_time",local_time);
    //ar & make_nvp("wavelength",wavelength);
    //ar & make_nvp("distance_uv",distance_uv);
    //ar & make_nvp("distance_visible",distance_visible);
    //ar & make_nvp("position",position);
    //ar & make_nvp("direction",direction);
    //ar & make_nvp("photon_source",photon_source);


}


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
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(CCMMCPE);


