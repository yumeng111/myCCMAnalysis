#ifndef HESodiumEvent_H_INCLUDED
#define HESodiumEvent_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Vector.h>

static const unsigned hesodiumevent_version_ = 0;

class HESodiumEvent : public I3FrameObject {
    public:

    I3Position photon_vertex;
    I3Position electron_vertex;
    I3Position positron_vertex;

    SET_LOGGER("HESodiumEvent");

    bool operator==(const HESodiumEvent& rhs) const {
        return photon_vertex == rhs.photon_vertex
            && electron_vertex == rhs.electron_vertex
            && positron_vertex == rhs.positron_vertex;
    }


  HESodiumEvent(I3Position photon_vertex_ = I3Position(0.0, 0.0, 0.0),
                 I3Position electron_vertex_ = I3Position(0.0, 0.0, 0.0), 
                 I3Position positron_vertex_ = I3Position(0.0, 0.0, 0.0)):
                photon_vertex(photon_vertex_),
                electron_vertex(electron_vertex_),
                positron_vertex(positron_vertex_) {
    }


    std::ostream& Print(std::ostream&) const;
    
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version) {
        if (version>hesodiumevent_version_)
            log_fatal("Attempting to read version %u from file but running version %u of HESodiumEvent class.",
                    version,hesodiumevent_version_);
        ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
        ar & make_nvp("photon_vertex",photon_vertex);
        ar & make_nvp("electron_vertex",electron_vertex);
        ar & make_nvp("positron_vertex",positron_vertex);
    }

};

I3_CLASS_VERSION(HESodiumEvent,hesodiumevent_version_);
typedef I3Vector<HESodiumEvent> HESodiumEventSeries;

std::ostream& operator<<(std::ostream&, const HESodiumEvent&);
I3_POINTER_TYPEDEFS(HESodiumEvent);
I3_POINTER_TYPEDEFS(HESodiumEventSeries);
I3_DEFAULT_NAME(HESodiumEvent);
#endif
