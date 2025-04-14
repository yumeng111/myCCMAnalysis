#include <simclasses/CCMMCPE.h>
#include <ostream>
#include <icetray/serialization.h>

I3_SERIALIZABLE(CCMMCPESeriesMap);
I3_SERIALIZABLE(CCMMCPESeries);

const std::unordered_map<CCMMCPE::PhotonSource, std::string> CCMMCPE::photonSourceToProcessName = {{CCMMCPE::PhotonSource::Unknown, "Unknown"},
                                                                                                      {CCMMCPE::PhotonSource::Scintillation, "Scintillation"},
                                                                                                      {CCMMCPE::PhotonSource::Cherenkov, "Cerenkov"},
                                                                                                      {CCMMCPE::PhotonSource::OpWLS, "OpWLS"}};
template <class Archive>
void CCMMCPE::save(Archive& ar, unsigned version) const {
    if (version > ccmmcpe_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMMCPE class.", version, ccmmcpe_version_);

    ar & make_nvp("parent_id",parent_id);
    ar & make_nvp("track_id",track_id);
    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("time",time);
    //ar & make_nvp("calculated_time",calculated_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("distance_uv",distance_uv);
    ar & make_nvp("original_wavelength",original_wavelength);
    ar & make_nvp("distance_visible",distance_visible);
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
        n_photons_per_wls = {};
        wls_loc = WLSLocationSeries();
    } else {
        if (ccmmcpe_version_ == 1) {
            size_t n_wls = 0;
            ar & make_nvp("n_wls",n_wls);
        }
        ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
        ar & make_nvp("wls_loc",wls_loc);
    }
    ar & make_nvp("time",time);
    //ar & make_nvp("calculated_time",calculated_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("distance_uv",distance_uv);
    ar & make_nvp("original_wavelength",original_wavelength);
    ar & make_nvp("distance_visible",distance_visible);
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
        << "\n  NPhotons Per WLS :" << n_photons_per_wls
        << "\n  WLS Location :" << wls_loc
        << "\n  Time :" << time
        //<< "\n  Calculated Time :" << calculated_time
        << "\n  Wavelength :" << wavelength
        << "\n  Distance UV :" << distance_uv
        << "\n  Original Wavelength :" << original_wavelength
        << "\n  Distance Visible :" << distance_visible
        //<< "\n  Calculated Distance UV :" << calculated_distance_uv
        //<< "\n  Calculated Distance Visible :" << calculated_distance_visible
        //<< "\n  Position :" << position
        //<< "\n  Direction :" << direction
        << "\n PhotonSource :" << photonSourceToProcessName.at(photon_source)
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(CCMMCPE);


