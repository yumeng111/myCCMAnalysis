#ifndef PhotonSummary_H_INCLUDED
#define PhotonSummary_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Vector.h"

static const unsigned ccmmuonsummary_version_ = 0;

class PhotonSummary : public I3FrameObject {
    public:
    enum class CreationProcess : int8_t {
        Unknown = 0,
        Scintillation = 1,
        Cerenkov = 2,
        OpWLS = 3
    };

    // things I want to save about a simulation photon:
    // parent id
    // track id
    // time
    // position
    // creation process 
    size_t parent_id;
    size_t track_id;
    float time;
    I3Position position; 
    CreationProcess creation_process;

    SET_LOGGER("PhotonSummary");

    bool operator==(const PhotonSummary& rhs) const {
        return parent_id == rhs.parent_id 
            && track_id == rhs.track_id
            && time == rhs.time
            && position == rhs.position
            && creation_process == rhs.creation_process;
    }


  PhotonSummary(size_t parent_id_ = 0, size_t track_id_ = 0, float time_ = 0, I3Position position_ = I3Position(0.0, 0.0, 0.0), 
                 CreationProcess creation_process_ = PhotonSummary::CreationProcess::Unknown):
                parent_id(parent_id_), track_id(track_id_), time(time_), position(position_), creation_process(creation_process_) {
    }


    std::ostream& Print(std::ostream&) const;

    private:
        static const std::unordered_map<PhotonSummary::CreationProcess, std::string> CreationProcesstoName;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>ccmmuonsummary_version_)
            log_fatal("Attempting to read version %u from file but running version %u of PhotonSummary class.",
                    version,ccmmuonsummary_version_);
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
        ar & make_nvp("parent_id",parent_id);
        ar & make_nvp("track_id",track_id);
        ar & make_nvp("time",time);
        ar & make_nvp("position",position);
        ar & make_nvp("creation_process",creation_process);
    }

};

I3_CLASS_VERSION(PhotonSummary,ccmmuonsummary_version_);

typedef I3Vector<PhotonSummary> PhotonSummarySeries;

std::ostream& operator<<(std::ostream&, const PhotonSummary&);
I3_POINTER_TYPEDEFS(PhotonSummary);
I3_POINTER_TYPEDEFS(PhotonSummarySeries);
I3_DEFAULT_NAME(PhotonSummary);
#endif
