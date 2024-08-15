#include "simclasses/PhotonSummary.h"
#include <icetray/serialization.h>
#include <ostream>

const std::unordered_map<PhotonSummary::PhotonSource, std::string> PhotonSummary::photonSourceToProcessName = {{PhotonSummary::PhotonSource::Unknown, "Unknown"},
                                                                                                               {PhotonSummary::PhotonSource::Scintillation, "Scintillation"},
                                                                                                               {PhotonSummary::PhotonSource::Cerenkov, "Cerenkov"}};


template <class Archive>
void PhotonSummary::save(Archive& ar, unsigned version) const {
    if (version > photonsummary_version_)
        log_fatal("Attempting to save version %u from file but running version %u of PhotonSummary class.", version, photonsummary_version_);


    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    ar & make_nvp("g4_distance_visible",g4_distance_visible);
    ar & make_nvp("calculated_distance_uv",calculated_distance_uv);
    ar & make_nvp("calculated_distance_visible",calculated_distance_visible);
    ar & make_nvp("g4_time",g4_time);
    ar & make_nvp("calculated_time",calculated_time);
    ar & make_nvp("n_wls",n_wls);
    ar & make_nvp("photon_source",photon_source);
}

template <class Archive>
void PhotonSummary::load(Archive& ar, unsigned version) {
    if (version > photonsummary_version_)
        log_fatal("Attempting to read version %u from file but running version %u of PhotonSummary class.", version, photonsummary_version_);

    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    ar & make_nvp("g4_distance_visible",g4_distance_visible);
    ar & make_nvp("calculated_distance_uv",calculated_distance_uv);
    ar & make_nvp("calculated_distance_visible",calculated_distance_visible);
    ar & make_nvp("g4_time",g4_time);
    ar & make_nvp("calculated_time",calculated_time);
    ar & make_nvp("n_wls",n_wls);

    if (photonsummary_version_ > 0){
        ar & make_nvp("photon_source",photon_source);
    } else {
        photon_source = PhotonSummary::PhotonSource::Unknown;
    }
}


std::ostream& PhotonSummary::Print(std::ostream& os) const{
    os << "[ PhotonSummary::"
        << "\n  G4 Distance UV  :" << g4_distance_uv
        << "\n  G4 Distance Visible :" << g4_distance_visible
        << "\n  Calculated Distance UV  :" << calculated_distance_uv
        << "\n  Calculated Distance Visible :" << calculated_distance_visible
        << "\n  G4 Global Time :" << g4_time
        << "\n  Calculated Time :" << calculated_time
        << "\n  Number WLS :" << n_wls
        << "\n  Photon Source :" << photonSourceToProcessName.at(photon_source)
        << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& oss, PhotonSummary const & bcm) {
    return(bcm.Print(oss));
}

std::ostream& operator<<(std::ostream& oss, PhotonSummary & bcm) {
    return(bcm.Print(oss));
}

I3_SPLIT_SERIALIZABLE(PhotonSummary);
I3_SERIALIZABLE(PhotonSummarySeries);

