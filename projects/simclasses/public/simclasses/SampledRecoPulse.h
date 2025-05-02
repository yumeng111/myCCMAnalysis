/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef SampledRecoPulse_H_INCLUDED
#define SampledRecoPulse_H_INCLUDED

#include "dataclasses/I3Position.h"
#include "dataclasses/I3Direction.h"
#include "dataclasses/Utility.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include <simclasses/WLSLocation.h>

#include <string>
#include <iostream>
#include <sstream>

#include <icetray/CCMPMTKey.h>

static const unsigned sampledrecopulse_version_ = 0;

/**
 * @brief SampledRecoPulse struct that stores the photon arrival time
 * (i.e.PE creation time), number of PE (for binning), and
 * the IDs of the particle that created this.
 */

struct SampledRecoPulse {
    float light_time;
    float afterpulse_time;
    float reco_time;
    float amplitude;
    bool in_gev;

    SET_LOGGER("SampledRecoPulse");

    bool operator==(const SampledRecoPulse& rhs) const {
        return light_time == rhs.light_time
            && afterpulse_time == rhs.afterpulse_time
            && reco_time == rhs.reco_time
            && amplitude == rhs.amplitude
            && in_gev == rhs.in_gev
            ;
    }

    SampledRecoPulse(const SampledRecoPulse&) = default;

    SampledRecoPulse(float light_time_ = 0,
                     float afterpulse_time_ = 0,
                     float reco_time_ = 0,
                     float amplitude_ = 0,
                     bool in_gev_ = false):
                     light_time(light_time_),
                     afterpulse_time(afterpulse_time_),
                     reco_time(reco_time_),
                     amplitude(amplitude_),
                     in_gev(in_gev_) {
    }


    std::ostream& Print(std::ostream&) const;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

};

I3_CLASS_VERSION(SampledRecoPulse,sampledrecopulse_version_);

typedef I3Vector<SampledRecoPulse> SampledRecoPulseSeries;
typedef I3Map<CCMPMTKey, SampledRecoPulseSeries > SampledRecoPulseSeriesMap;

std::ostream& operator<<(std::ostream&, const SampledRecoPulse&);

I3_POINTER_TYPEDEFS(SampledRecoPulse);
I3_POINTER_TYPEDEFS(SampledRecoPulseSeries);
I3_POINTER_TYPEDEFS(SampledRecoPulseSeriesMap);
#endif
