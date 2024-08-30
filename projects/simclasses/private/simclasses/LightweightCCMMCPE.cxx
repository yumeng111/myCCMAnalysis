#include <simclasses/LightweightCCMMCPE.h>
#include <ostream>
#include <icetray/serialization.h>

I3_SERIALIZABLE(LightweightCCMMCPESeriesMap);
I3_SERIALIZABLE(LightweightCCMMCPESeries);

const std::unordered_map<LightweightCCMMCPE::PhotonSource, std::string> LightweightCCMMCPE::photonSourceToProcessName = {{LightweightCCMMCPE::PhotonSource::Unknown, "Unknown"},
                                                                                                                         {LightweightCCMMCPE::PhotonSource::Scintillation, "Scintillation"},
                                                                                                                         {LightweightCCMMCPE::PhotonSource::Cerenkov, "Cerenkov"},
                                                                                                                         {LightweightCCMMCPE::PhotonSource::OpWLS, "OpWLS"}};
template <class Archive>
void LightweightCCMMCPE::save(Archive& ar, unsigned version) const {
    if (version > lightweightccmmcpe_version_)
        log_fatal("Attempting to save version %u from file but running version %u of LightweightCCMMCPE class.", version, lightweightccmmcpe_version_);

    ar & make_nvp("g4_time",g4_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    ar & make_nvp("photon_source",photon_source);
}

template <class Archive>
void LightweightCCMMCPE::load(Archive& ar, unsigned version) {
    if (version > lightweightccmmcpe_version_)
        log_fatal("Attempting to read version %u from file but running version %u of LightweightCCMMCPE class.", version, lightweightccmmcpe_version_);
    
    ar & make_nvp("g4_time",g4_time);
    ar & make_nvp("wavelength",wavelength);
    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    ar & make_nvp("photon_source",photon_source);

}


std::ostream& operator<<(std::ostream& os, const LightweightCCMMCPE& pe) {
    return(pe.Print(os));
}

std::ostream& LightweightCCMMCPE::Print(std::ostream& os) const{
    os << "[ LightweightCCMMCPE::"
        << "\n  G4 Time :" << g4_time
        << "\n  Wavelength :" << wavelength
        << "\n  G4 Distance UV :" << g4_distance_uv
        << "\n PhotonSource :" << photonSourceToProcessName.at(photon_source) 
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(LightweightCCMMCPE);


