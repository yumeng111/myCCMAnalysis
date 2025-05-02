
#include <simclasses/SampledRecoPulse.h>
#include <ostream>
#include <icetray/serialization.h>

I3_SERIALIZABLE(SampledRecoPulseSeriesMap);
I3_SERIALIZABLE(SampledRecoPulseSeries);

template <class Archive>
void SampledRecoPulse::save(Archive& ar, unsigned version) const {
    if (version > sampledrecopulse_version_)
        log_fatal("Attempting to save version %u from file but running version %u of SampledRecoPulse class.", version, sampledrecopulse_version_);

    ar & make_nvp("light_time",light_time);
    ar & make_nvp("afterpulse_time",afterpulse_time);
    ar & make_nvp("reco_time",reco_time);
    ar & make_nvp("amplitude",amplitude);
    ar & make_nvp("in_gev",in_gev);
}

template <class Archive>
void SampledRecoPulse::load(Archive& ar, unsigned version) {
    if (version > sampledrecopulse_version_)
        log_fatal("Attempting to read version %u from file but running version %u of SampledRecoPulse class.", version, sampledrecopulse_version_);

    ar & make_nvp("light_time",light_time);
    ar & make_nvp("afterpulse_time",afterpulse_time);
    ar & make_nvp("reco_time",reco_time);
    ar & make_nvp("amplitude",amplitude);
    ar & make_nvp("in_gev",in_gev);
}


std::ostream& operator<<(std::ostream& os, const SampledRecoPulse& pe) {
    return(pe.Print(os));
}

std::ostream& SampledRecoPulse::Print(std::ostream& os) const{
    os << "[ SampledRecoPulse::"
        << "\n  Light Time :" << light_time
        << "\n  Afterpulse Time :" << afterpulse_time
        << "\n  Reco Time :" << reco_time
        << "\n  Amplitude :" << amplitude
        << "\n  In GeV    :" << in_gev
        << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(SampledRecoPulse);


