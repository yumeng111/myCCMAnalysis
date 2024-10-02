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
#include <simclasses/WLSLocation.h>

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

    // things we want to save about a photon hitting our pmts in simulation
    std::vector<size_t> n_photons_per_wls;
    WLSLocationSeries wls_loc; 
    float g4_time;
    float wavelength; // wavelength of photon
    float g4_distance_uv; // g4_distances travelled as uv photon

    SET_LOGGER("LightweightCCMMCPE");

    bool operator==(const LightweightCCMMCPE& rhs) const {
        return n_photons_per_wls == rhs.n_photons_per_wls
            && wls_loc == rhs.wls_loc
            && g4_time == rhs.g4_time
            && wavelength == rhs.wavelength
            && g4_distance_uv == rhs.g4_distance_uv;
    }

    LightweightCCMMCPE(const LightweightCCMMCPE&) = default;

    LightweightCCMMCPE(std::vector<size_t> n_photons_per_wls_ = {0},
                       WLSLocationSeries wls_loc_ = WLSLocationSeries(),
                       float g4_time_ = 0,
                       float wavelength_ = 0,
                       float g4_distance_uv_ = 0):
                       n_photons_per_wls(n_photons_per_wls_),
                       wls_loc(wls_loc_),
                       g4_time(g4_time_),
                       wavelength(wavelength_),
                       g4_distance_uv(g4_distance_uv_){
    }


    std::ostream& Print(std::ostream&) const;

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
