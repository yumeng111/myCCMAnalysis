#ifndef PhotonSummary_H_INCLUDED
#define PhotonSummary_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include "dataclasses/I3Vector.h"

static const unsigned ccmmuonsummary_version_ = 0;

class PhotonSummary : public I3FrameObject {
    public:
    // things I want to save about a simulation photon:
    // distance travelled as UV photon
    // distance travelled as visible photon
    // number of WLS

    float distance_uv;
    float distance_visible;
    size_t n_wls;

    SET_LOGGER("PhotonSummary");

    bool operator==(const PhotonSummary& rhs) const {
        return distance_uv == rhs.distance_uv 
            && distance_visible == rhs.distance_visible
            && n_wls == rhs.n_wls;
    }


  PhotonSummary(float distance_uv_ = 0, float distance_visible_ = 0, size_t n_wls_ = 0): 
                distance_uv(distance_uv_), distance_visible(distance_visible_), n_wls(n_wls_) {
    }


    std::ostream& Print(std::ostream&) const;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>ccmmuonsummary_version_)
            log_fatal("Attempting to read version %u from file but running version %u of PhotonSummary class.",
                    version,ccmmuonsummary_version_);
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
        ar & make_nvp("distance_uv",distance_uv);
        ar & make_nvp("distance_visible",distance_visible);
        ar & make_nvp("n_wls",n_wls);
    }

};

I3_CLASS_VERSION(PhotonSummary,ccmmuonsummary_version_);

typedef I3Vector<PhotonSummary> PhotonSummarySeries;

std::ostream& operator<<(std::ostream&, const PhotonSummary&);
I3_POINTER_TYPEDEFS(PhotonSummary);
I3_POINTER_TYPEDEFS(PhotonSummarySeries);
I3_DEFAULT_NAME(PhotonSummary);
#endif
