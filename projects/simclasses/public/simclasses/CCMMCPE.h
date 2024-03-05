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
        Cherenkov = 2
    };

    //static const std::map<PhotonSource, std::string> PhotonSourceNames = {
    //#define X(a, b) {CCMMCPE::PhotonSource::a , #a },
    //#undef X
    //};
    //static const std::map<PhotonSource, std::string> PhotonSourceNames = {{CCMMCPE::PhotonSource::Unknown, "Unknown"}, 
    //                                                                      {CCMMCPE::PhotonSource::Scintillation, "Scintillation"},
    //                                                                      {CCMMCPE::PhotonSource::Cherenkov, "Cherenkov"}}; 

    // things we want to save about a photon hitting our pmts in simulation
    float time; // photon hit time
    float wavelength; // wavelength of photon
    I3Position position; // hit position on PMT
    I3Direction direction; // hit direction on PMT
    PhotonSource photon_source; // true if photon produced from scintillation, false if photon produced via cherenkov radiation 

    SET_LOGGER("CCMMCPE");

    bool operator==(const CCMMCPE& rhs) const {
        return time == rhs.time
            && wavelength == rhs.wavelength
            && position == rhs.position
            && direction == rhs.direction
            && photon_source == rhs.photon_source;
    }


  CCMMCPE(float time_ = 0, float wavelength_ = 0, I3Position position_ = I3Position(0.0, 0.0, 0.0), I3Direction direction_ = I3Direction(0.0, 0.0, 0.0), PhotonSource photon_source_ = CCMMCPE::PhotonSource::Unknown):
        time(time_), wavelength(wavelength_), position(position_), direction(direction_), photon_source(photon_source_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>ccmmcpe_version_)
            log_fatal("Attempting to read version %u from file but running version %u of CCMMCPE class.",
                    version,ccmmcpe_version_);
        ar & make_nvp("time",time);
        ar & make_nvp("wavelength",wavelength);
        ar & make_nvp("position",position);
        ar & make_nvp("direction",direction);
        ar & make_nvp("photon_source",photon_source);
    }

};

I3_CLASS_VERSION(CCMMCPE,ccmmcpe_version_);

typedef std::vector<CCMMCPE> CCMMCPESeries;
typedef I3Map<CCMPMTKey, CCMMCPESeries > CCMMCPESeriesMap;

std::ostream& operator<<(std::ostream&, const CCMMCPE&);

I3_POINTER_TYPEDEFS(CCMMCPE);
I3_POINTER_TYPEDEFS(CCMMCPESeries);
I3_POINTER_TYPEDEFS(CCMMCPESeriesMap);
#endif
