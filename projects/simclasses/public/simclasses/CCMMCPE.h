/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef CCMMCPE_H_INCLUDED
#define CCMMCPE_H_INCLUDED

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

static const unsigned ccmmcpe_version_ = 2;

/**
 * @brief CCMMCPE struct that stores the photon arrival time
 * (i.e.PE creation time), number of PE (for binning), and
 * the IDs of the particle that created this.
 */

struct CCMMCPE {
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Scintillation = 1,
        Cherenkov = 2,
        OpWLS = 3
    };

    // things we want to save about a photon hitting our pmts in simulation
    size_t parent_id;
    size_t track_id;
    std::vector<size_t> n_photons_per_wls;
    WLSLocationSeries wls_loc;
    float time;
    //float calculated_time;
    float wavelength; // wavelength of photon
    float distance_uv; // distances travelled as uv photon
    float original_wavelength;
    float distance_visible; // distance travelled as visible photon
    //float calculated_distance_uv; // calculated_distances travelled as uv photon
    //float calculated_distance_visible; // calculated_distance travelled as visible photon
    //I3Position position; // hit position on PMT
    //I3Direction direction; // hit direction on PMT
    PhotonSource photon_source;

    SET_LOGGER("CCMMCPE");

    bool operator==(const CCMMCPE& rhs) const {
        return parent_id == rhs.parent_id
            && track_id == rhs.track_id
            && n_photons_per_wls == rhs.n_photons_per_wls
            && wls_loc == rhs.wls_loc
            && time == rhs.time
            //&& calculated_time == rhs.calculated_time
            && wavelength == rhs.wavelength
            && distance_uv == rhs.distance_uv
            && original_wavelength == rhs.original_wavelength
            && distance_visible == rhs.distance_visible
            //&& calculated_distance_uv == rhs.calculated_distance_uv
            //&& calculated_distance_visible == rhs.calculated_distance_visible
            //&& position == rhs.position
            //&& direction == rhs.direction
            && photon_source == rhs.photon_source;
    }

    CCMMCPE(const CCMMCPE&) = default;

    CCMMCPE(size_t parent_id_ = 0,
            size_t track_id_ = 0,
            std::vector<size_t> n_photons_per_wls_ = {0},
            WLSLocationSeries wls_loc_ = WLSLocationSeries(),
            float time_ = 0,
            //float calculated_time_ = 0,
            float wavelength_ = 0,
            float distance_uv_ = 0,
            float original_wavelength_ = 0,
            float distance_visible_ = 0,
            //float calculated_distance_uv_ = 0,
            //float calculated_distance_visible_ = 0,
            //I3Position position_ = I3Position(0.0, 0.0, 0.0),
            //I3Direction direction_ = I3Direction(0.0, 0.0, 0.0),
            PhotonSource photon_source_ = CCMMCPE::PhotonSource::Unknown):
            parent_id(parent_id_), track_id(track_id_),
            n_photons_per_wls(n_photons_per_wls_), wls_loc(wls_loc_),
            time(time_), //calculated_time(calculated_time_),
            wavelength(wavelength_),
            distance_uv(distance_uv_),
            original_wavelength(original_wavelength_),
            distance_visible(distance_visible_),
            //calculated_distance_uv(calculated_distance_uv_), calculated_distance_visible(calculated_distance_visible_),
            //position(position_), direction(direction_),
            photon_source(photon_source_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
    static const std::unordered_map<CCMMCPE::PhotonSource, std::string> photonSourceToProcessName;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

};

I3_CLASS_VERSION(CCMMCPE,ccmmcpe_version_);

typedef I3Vector<CCMMCPE> CCMMCPESeries;
typedef I3Map<CCMPMTKey, CCMMCPESeries > CCMMCPESeriesMap;

std::ostream& operator<<(std::ostream&, const CCMMCPE&);

I3_POINTER_TYPEDEFS(CCMMCPE);
I3_POINTER_TYPEDEFS(CCMMCPESeries);
I3_POINTER_TYPEDEFS(CCMMCPESeriesMap);
#endif
