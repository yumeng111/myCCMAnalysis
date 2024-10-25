#include <simclasses/LightweightCCMMCPE.h>
#include <ostream>
#include <icetray/serialization.h>

I3_SERIALIZABLE(LightweightCCMMCPESeriesMap);
I3_SERIALIZABLE(LightweightCCMMCPESeries);

template <class Archive>
void LightweightCCMMCPE::save(Archive& ar, unsigned version) const {
    if (version > lightweightccmmcpe_version_)
        log_fatal("Attempting to save version %u from file but running version %u of LightweightCCMMCPE class.", version, lightweightccmmcpe_version_);

    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("g4_time",g4_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("g4_distance_uv",g4_distance_uv);
}

template <class Archive>
void LightweightCCMMCPE::load(Archive& ar, unsigned version) {
    if (version > lightweightccmmcpe_version_)
        log_fatal("Attempting to read version %u from file but running version %u of LightweightCCMMCPE class.", version, lightweightccmmcpe_version_);
    
    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("g4_time",g4_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("g4_distance_uv",g4_distance_uv);

}


std::ostream& operator<<(std::ostream& os, const LightweightCCMMCPE& pe) {
    return(pe.Print(os));
}

std::ostream& LightweightCCMMCPE::Print(std::ostream& os) const{
    os << "[ LightweightCCMMCPE::"
        << "\n  NPhotons Per WLS :" << n_photons_per_wls 
        << "\n  WLS Location :" << wls_loc 
        << "\n  G4 Time :" << g4_time
        << "\n  Wavelength :" << wavelength
        << "\n  G4 Distance UV :" << g4_distance_uv
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(LightweightCCMMCPE);


