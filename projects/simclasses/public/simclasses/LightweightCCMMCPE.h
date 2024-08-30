/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef LightweightCCMMCPE_H_INCLUDED
#define LightweightCCMMCPE_H_INCLUDED

#include "dataclasses/Utility.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"

#include <string>
#include <iostream>
#include <sstream>

#include <icetray/CCMPMTKey.h>

static const unsigned lightweightccmmcpe_version_ = 0;

/**
 * @brief LightweightCCMMCPE struct that stores the photon arrival time
 * (i.e.PE creation time), number of PE (for binning), and
 * the IDs of the particle that created this.
 */

struct LightweightCCMMCPE {
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Scintillation = 1,
        Cerenkov = 2,
        OpWLS = 3
    };

    // things we want to save about a photon hitting our pmts in simulation
    float g4_time;
    float wavelength; // wavelength of photon
    float g4_distance_uv; // g4_distances travelled as uv photon
    PhotonSource photon_source;

    SET_LOGGER("LightweightCCMMCPE");

    bool operator==(const LightweightCCMMCPE& rhs) const {
        return g4_time == rhs.g4_time
            && wavelength == rhs.wavelength
            && g4_distance_uv == rhs.g4_distance_uv
            && photon_source == rhs.photon_source;
    }

    LightweightCCMMCPE(const LightweightCCMMCPE&) = default;

    LightweightCCMMCPE(float g4_time_ = 0,
            float wavelength_ = 0,
            float g4_distance_uv_ = 0,
            PhotonSource photon_source_ = LightweightCCMMCPE::PhotonSource::Unknown):
        g4_time(g4_time_), wavelength(wavelength_),
        g4_distance_uv(g4_distance_uv_), photon_source(photon_source_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
    static const std::unordered_map<LightweightCCMMCPE::PhotonSource, std::string> photonSourceToProcessName;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

};

I3_CLASS_VERSION(LightweightCCMMCPE,lightweightccmmcpe_version_);

typedef I3Vector<LightweightCCMMCPE> LightweightCCMMCPESeries;
typedef I3Map<CCMPMTKey, LightweightCCMMCPESeries > LightweightCCMMCPESeriesMap;

std::ostream& operator<<(std::ostream&, const LightweightCCMMCPE&);

I3_POINTER_TYPEDEFS(LightweightCCMMCPE);
I3_POINTER_TYPEDEFS(LightweightCCMMCPESeries);
I3_POINTER_TYPEDEFS(LightweightCCMMCPESeriesMap);
#endif
