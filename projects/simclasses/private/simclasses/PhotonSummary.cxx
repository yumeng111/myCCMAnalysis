#include "simclasses/PhotonSummary.h"
#include <icetray/serialization.h>
#include <ostream>

const std::unordered_map<PhotonSummary::PhotonSource, std::string> PhotonSummary::photonSourceToProcessName = {{PhotonSummary::PhotonSource::Unknown, "Unknown"},
                                                                                                               {PhotonSummary::PhotonSource::Scintillation, "Scintillation"},
                                                                                                               {PhotonSummary::PhotonSource::OpWLS, "OpWLS"},
                                                                                                               {PhotonSummary::PhotonSource::Cherenkov, "Cerenkov"}};

template <class Archive>
void PhotonSummary::save(Archive& ar, unsigned version) const {
    if (version > photonsummary_version_)
        log_fatal("Attempting to save version %u from file but running version %u of PhotonSummary class.", version, photonsummary_version_);


    ar & make_nvp("distance_uv",distance_uv);
    ar & make_nvp("original_wavelength",original_wavelength);
    ar & make_nvp("distance_visible",distance_visible);
    ar & make_nvp("time",time);
    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("photon_source",photon_source);
    ar & make_nvp("current_process",current_process);
}

template <class Archive>
void PhotonSummary::load(Archive& ar, unsigned version) {
    if (version > photonsummary_version_)
        log_fatal("Attempting to read version %u from file but running version %u of PhotonSummary class.", version, photonsummary_version_);

    ar & make_nvp("distance_uv",distance_uv);
    ar & make_nvp("original_wavelength",original_wavelength);
    ar & make_nvp("distance_visible",distance_visible);
    ar & make_nvp("time",time);
    if(version <= 2) {
        size_t n_wls;
        ar & make_nvp("n_wls",n_wls);
    }
    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("photon_source",photon_source);
    if(version <= 2) {
        PhotonSource temp_parent;
        ar & make_nvp("temp_parent",temp_parent);
    }
    ar & make_nvp("current_process",current_process);
}


std::ostream& PhotonSummary::Print(std::ostream& os) const{
    os << "[ PhotonSummary::"
        << "\n  Distance UV  :" << distance_uv
        << "\n  Original Wavelength :" << original_wavelength
        << "\n  Distance Visible :" << distance_visible
        << "\n  Global Time :" << time
        << "\n  Number Photons Per WLS :" << n_photons_per_wls
        << "\n  WLS Location :" << wls_loc
        << "\n  Photon Source :" << photonSourceToProcessName.at(photon_source)
        << "\n  Current Process :" << photonSourceToProcessName.at(current_process)
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

