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
    float time; // photon hit time
    float wavelength; // wavelength of photon
    I3Position position; // hit position on PMT
    I3Direction direction; // hit direction on PMT
    PhotonSource photon_source; 

    SET_LOGGER("CCMMCPE");

    bool operator==(const CCMMCPE& rhs) const {
        return parent_id == rhs.parent_id 
            && track_id == rhs.track_id
            && time == rhs.time
            && wavelength == rhs.wavelength
            && position == rhs.position
            && direction == rhs.direction
            && photon_source == rhs.photon_source;
    }


  CCMMCPE(size_t parent_id_ = 0,
          size_t track_id_ = 0,
          float time_ = 0,
          float wavelength_ = 0,
          I3Position position_ = I3Position(0.0, 0.0, 0.0),
          I3Direction direction_ = I3Direction(0.0, 0.0, 0.0),
          PhotonSource photon_source_ = CCMMCPE::PhotonSource::Unknown):
        parent_id(parent_id_), track_id(track_id_), time(time_), wavelength(wavelength_), position(position_), direction(direction_), photon_source(photon_source_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
        static const std::unordered_map<CCMMCPE::PhotonSource, std::string> photonSourceToProcessName;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>ccmmcpe_version_)
            log_fatal("Attempting to read version %u from file but running version %u of CCMMCPE class.",
                    version,ccmmcpe_version_);
        ar & make_nvp("parent_id",parent_id);
        ar & make_nvp("track_id",track_id);
        ar & make_nvp("time",time);
        ar & make_nvp("wavelength",wavelength);
        ar & make_nvp("position",position);
        ar & make_nvp("direction",direction);
        ar & make_nvp("photon_source",photon_source);
    }

};

I3_CLASS_VERSION(CCMMCPE,ccmmcpe_version_);

typedef I3Vector<CCMMCPE> CCMMCPESeries;
typedef I3Map<CCMPMTKey, CCMMCPESeries > CCMMCPESeriesMap;

std::ostream& operator<<(std::ostream&, const CCMMCPE&);

I3_POINTER_TYPEDEFS(CCMMCPE);
I3_POINTER_TYPEDEFS(CCMMCPESeries);
I3_POINTER_TYPEDEFS(CCMMCPESeriesMap);
#endif
