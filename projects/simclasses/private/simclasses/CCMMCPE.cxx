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
    ar & make_nvp("n_wls",n_wls);
    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("g4_time",g4_time);
    //ar & make_nvp("calculated_time",calculated_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    //ar & make_nvp("g4_distance_visible",g4_distance_visible);
    //ar & make_nvp("calculated_distance_uv",calculated_distance_uv);
    //ar & make_nvp("calculated_distance_visible",calculated_distance_visible);
    //ar & make_nvp("position",position);
    //ar & make_nvp("direction",direction);
    ar & make_nvp("photon_source",photon_source);
}

template <class Archive>
void CCMMCPE::load(Archive& ar, unsigned version) {
    if (version > ccmmcpe_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMMCPE class.", version, ccmmcpe_version_);
    
    ar & make_nvp("parent_id",parent_id);
    ar & make_nvp("track_id",track_id);
    if (ccmmcpe_version_ == 0){
        n_wls = 0;
        n_photons_per_wls = {};
        wls_loc = WLSLocationSeries();
    } else {
        ar & make_nvp("n_wls",n_wls);
        ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
        ar & make_nvp("wls_loc",wls_loc);
    }
    ar & make_nvp("g4_time",g4_time);
    //ar & make_nvp("calculated_time",calculated_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    //ar & make_nvp("g4_distance_visible",g4_distance_visible);
    //ar & make_nvp("calculated_distance_uv",calculated_distance_uv);
    //ar & make_nvp("calculated_distance_visible",calculated_distance_visible);
    //ar & make_nvp("position",position);
    //ar & make_nvp("direction",direction);
    ar & make_nvp("photon_source",photon_source);

}


std::ostream& operator<<(std::ostream& os, const CCMMCPE& pe) {
    return(pe.Print(os));
}

std::ostream& CCMMCPE::Print(std::ostream& os) const{
    //std::vector<std::string> wls_loc_string = {};
    //for (size_t w = 0; w < wls_loc.size(); w++){
    //    wls_loc_string.push_back(WLSLocation::wlsLocationToProcessName.at(wls_loc.at(w)));
    //}
    os << "[ CCMMCPE::"
        << "\n  ParentID :" << parent_id 
        << "\n  TrackID :" << track_id
        << "\n  NWLS :" << n_wls 
        << "\n  NPhotons Per WLS :" << n_photons_per_wls 
        << "\n  WLS Location :" << wls_loc 
        << "\n  G4 Time :" << g4_time
        //<< "\n  Calculated Time :" << calculated_time
        << "\n  Wavelength :" << wavelength
        << "\n  G4 Distance UV :" << g4_distance_uv
        //<< "\n  G4 Distance Visible :" << g4_distance_visible
        //<< "\n  Calculated Distance UV :" << calculated_distance_uv
        //<< "\n  Calculated Distance Visible :" << calculated_distance_visible
        //<< "\n  Position :" << position 
        //<< "\n  Direction :" << direction 
        << "\n PhotonSource :" << photonSourceToProcessName.at(photon_source) 
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(CCMMCPE);


