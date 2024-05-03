#ifndef PhotonYieldSummary_H_INCLUDED
#define PhotonYieldSummary_H_INCLUDED

#include "dataclasses/Utility.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"

#include <string>
#include <iostream>
#include <sstream>

#include <icetray/CCMPMTKey.h>

static const unsigned photonyieldsummary_version_ = 0;

class PhotonYieldSummary : public I3FrameObject {
    public:
    
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Vertex = 1,
        TPBFoil = 2,
    };

    float time; // photon hit time
    float yield; // number of photon hits 
    PhotonSource photon_source; 

    SET_LOGGER("PhotonYieldSummary");

    bool operator==(const PhotonYieldSummary& rhs) const {
        return time == rhs.time
            && yield == rhs.yield
            && photon_source == rhs.photon_source;
    }


  PhotonYieldSummary(float time_ = 0,
          float yield_ = 0,
          PhotonSource photon_source_ = PhotonYieldSummary::PhotonSource::Unknown):
        time(time_), yield(yield_), photon_source(photon_source_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
        static const std::unordered_map<PhotonYieldSummary::PhotonSource, std::string> photonSourceToProcessName;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>photonyieldsummary_version_)
            log_fatal("Attempting to read version %u from file but running version %u of photonyieldsummary class.",
                    version,photonyieldsummary_version_);
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
        ar & make_nvp("time",time);
        ar & make_nvp("yield",yield);
        ar & make_nvp("photon_source",photon_source);
    }

};

I3_CLASS_VERSION(PhotonYieldSummary,photonyieldsummary_version_);

typedef I3Vector<PhotonYieldSummary> PhotonYieldSummarySeries;
typedef I3Map<CCMPMTKey, PhotonYieldSummarySeries > PhotonYieldSummarySeriesMap;

std::ostream& operator<<(std::ostream&, const PhotonYieldSummary&);

I3_POINTER_TYPEDEFS(PhotonYieldSummary);
I3_POINTER_TYPEDEFS(PhotonYieldSummarySeries);
I3_POINTER_TYPEDEFS(PhotonYieldSummarySeriesMap);
#endif
