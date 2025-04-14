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
    ar & make_nvp("time",time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("distance_uv",distance_uv);
}

template <class Archive>
void LightweightCCMMCPE::load(Archive& ar, unsigned version) {
    if (version > lightweightccmmcpe_version_)
        log_fatal("Attempting to read version %u from file but running version %u of LightweightCCMMCPE class.", version, lightweightccmmcpe_version_);
    
    ar & make_nvp("n_photons_per_wls",n_photons_per_wls);
    ar & make_nvp("wls_loc",wls_loc);
    ar & make_nvp("time",time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("distance_uv",distance_uv);

}


std::ostream& operator<<(std::ostream& os, const LightweightCCMMCPE& pe) {
    return(pe.Print(os));
}

std::ostream& LightweightCCMMCPE::Print(std::ostream& os) const{
    os << "[ LightweightCCMMCPE::"
        << "\n  NPhotons Per WLS :" << n_photons_per_wls 
        << "\n  WLS Location :" << wls_loc 
        << "\n  Time :" << time
        << "\n  Wavelength :" << wavelength
        << "\n  Distance UV :" << distance_uv
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(LightweightCCMMCPE);


