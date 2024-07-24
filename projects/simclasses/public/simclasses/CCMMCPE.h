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

#include <string>
#include <iostream>
#include <sstream>

#include <icetray/CCMPMTKey.h>

static const unsigned ccmmcpe_version_ = 0;

/**
 * @brief CCMMCPE struct that stores the photon arrival time
 * (i.e.PE creation time), number of PE (for binning), and
 * the IDs of the particle that created this.
 */

struct CCMMCPE {
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Scintillation = 1,
        Cerenkov = 2,
        OpWLS = 3
    };

    // things we want to save about a photon hitting our pmts in simulation
    size_t parent_id;
    size_t track_id;
    float g4_time; 
    float calculated_time; 
    float wavelength; // wavelength of photon
    float g4_distance_uv; // g4_distances travelled as uv photon
    float g4_distance_visible; // g4_distance travelled as visible photon
    float calculated_distance_uv; // calculated_distances travelled as uv photon
    float calculated_distance_visible; // calculated_distance travelled as visible photon
    I3Position position; // hit position on PMT
    I3Direction direction; // hit direction on PMT
    PhotonSource photon_source; 

    SET_LOGGER("CCMMCPE");

    bool operator==(const CCMMCPE& rhs) const {
        return parent_id == rhs.parent_id 
            && track_id == rhs.track_id
            && g4_time == rhs.g4_time
            && calculated_time == rhs.calculated_time
            && wavelength == rhs.wavelength
            && g4_distance_uv == rhs.g4_distance_uv
            && g4_distance_visible == rhs.g4_distance_visible
            && calculated_distance_uv == rhs.calculated_distance_uv
            && calculated_distance_visible == rhs.calculated_distance_visible
            && position == rhs.position
            && direction == rhs.direction
            && photon_source == rhs.photon_source;
    }


  CCMMCPE(size_t parent_id_ = 0,
          size_t track_id_ = 0,
          float g4_time_ = 0,
          float calculated_time_ = 0,
          float wavelength_ = 0,
          float g4_distance_uv_ = 0,
          float g4_distance_visible_ = 0,
          float calculated_distance_uv_ = 0,
          float calculated_distance_visible_ = 0,
          I3Position position_ = I3Position(0.0, 0.0, 0.0),
          I3Direction direction_ = I3Direction(0.0, 0.0, 0.0),
          PhotonSource photon_source_ = CCMMCPE::PhotonSource::Unknown):
        parent_id(parent_id_), track_id(track_id_), g4_time(g4_time_), calculated_time(calculated_time_), wavelength(wavelength_), 
        g4_distance_uv(g4_distance_uv_), g4_distance_visible(g4_distance_visible_),
        calculated_distance_uv(calculated_distance_uv_), calculated_distance_visible(calculated_distance_visible_),
        position(position_), direction(direction_), photon_source(photon_source_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
        static const std::unordered_map<CCMMCPE::PhotonSource, std::string> photonSourceToProcessName;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

    //friend class icecube::serialization::access;
    //template <class Archive> void serialize(Archive & ar, const unsigned version) {
    //    if (version>ccmmcpe_version_)
    //        log_fatal("Attempting to read version %u from file but running version %u of CCMMCPE class.",
    //                version,ccmmcpe_version_);
    //    ar & make_nvp("parent_id",parent_id);
    //    ar & make_nvp("track_id",track_id);
    //    ar & make_nvp("g4_time",g4_time);
    //    ar & make_nvp("calculated_time",calculated_time);
    //    ar & make_nvp("wavelength",wavelength);
    //    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    //    ar & make_nvp("g4_distance_visible",g4_distance_visible);
    //    ar & make_nvp("calculated_distance_uv",calculated_distance_uv);
    //    ar & make_nvp("calculated_distance_visible",calculated_distance_visible);
    //    ar & make_nvp("position",position);
    //    ar & make_nvp("direction",direction);
    //    ar & make_nvp("photon_source",photon_source);
    //}

};

I3_CLASS_VERSION(CCMMCPE,ccmmcpe_version_);

typedef I3Vector<CCMMCPE> CCMMCPESeries;
typedef I3Map<CCMPMTKey, CCMMCPESeries > CCMMCPESeriesMap;

std::ostream& operator<<(std::ostream&, const CCMMCPE&);

I3_POINTER_TYPEDEFS(CCMMCPE);
I3_POINTER_TYPEDEFS(CCMMCPESeries);
I3_POINTER_TYPEDEFS(CCMMCPESeriesMap);
#endif
