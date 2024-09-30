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

static const unsigned photonsummary_version_ = 2;

class PhotonSummary : public I3FrameObject {
    public:
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Scintillation = 1,
        Cerenkov = 2,
        OpWLS = 3
    };

    float g4_distance_uv; // record distance as given by GetStepLength
    float g4_distance_visible;
    float calculated_distance_uv; // record distance as calculated using delta G4GlobalTime
    float calculated_distance_visible;
    float g4_time; // as as reported using g4 global time
    float calculated_time; // time as calculated using g4 distance
    size_t n_wls;
    std::vector<size_t> n_photons_per_wls;
    PhotonSource photon_source; 
    PhotonSource temp_parent; 
    PhotonSource current_process; 

    SET_LOGGER("PhotonSummary");

    bool operator==(const PhotonSummary& rhs) const {
        return g4_distance_uv == rhs.g4_distance_uv 
            && g4_distance_visible == rhs.g4_distance_visible
            && calculated_distance_uv == rhs.calculated_distance_uv 
            && calculated_distance_visible == rhs.calculated_distance_visible
            && g4_time == rhs.g4_time
            && calculated_time == rhs.calculated_time
            && n_wls == rhs.n_wls
            && n_photons_per_wls == rhs.n_photons_per_wls
            && photon_source == rhs.photon_source
            && temp_parent == rhs.temp_parent
            && current_process == rhs.current_process;
    }


  PhotonSummary(float g4_distance_uv_ = 0, float g4_distance_visible_ = 0, 
                float calculated_distance_uv_ = 0, float calculated_distance_visible_ = 0,
                float g4_time_ = 0, float calculated_time_ = 0, size_t n_wls_ = 0, std::vector<size_t> n_photons_per_wls_ = {0},
                PhotonSource photon_source_ = PhotonSummary::PhotonSource::Unknown, PhotonSource temp_parent_ = PhotonSummary::PhotonSource::Unknown,
                PhotonSource current_process_ = PhotonSummary::PhotonSource::Unknown): 
                g4_distance_uv(g4_distance_uv_), g4_distance_visible(g4_distance_visible_), 
                calculated_distance_uv(calculated_distance_uv_), calculated_distance_visible(calculated_distance_visible_),
                g4_time(g4_time_), calculated_time(calculated_time_), n_wls(n_wls_), n_photons_per_wls(n_photons_per_wls_), photon_source(photon_source_),
                temp_parent(temp_parent_), current_process(current_process_) {
    }


    std::ostream& Print(std::ostream&) const;
    
    private:
    static const std::unordered_map<PhotonSummary::PhotonSource, std::string> photonSourceToProcessName;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

    //friend class icecube::serialization::access;
    //template <class Archive> void serialize(Archive & ar, const unsigned version) {
    //    if (version>photonsummary_version_)
    //        log_fatal("Attempting to read version %u from file but running version %u of PhotonSummary class.",
    //                version,photonsummary_version_);
    //    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    //    ar & make_nvp("g4_distance_uv",g4_distance_uv);
    //    ar & make_nvp("g4_distance_visible",g4_distance_visible);
    //    ar & make_nvp("calculated_distance_uv",calculated_distance_uv);
    //    ar & make_nvp("calculated_distance_visible",calculated_distance_visible);
    //    ar & make_nvp("g4_time",g4_time);
    //    ar & make_nvp("calculated_time",calculated_time);
    //    ar & make_nvp("n_wls",n_wls);
    //}

};

I3_CLASS_VERSION(PhotonSummary,photonsummary_version_);

typedef I3Vector<PhotonSummary> PhotonSummarySeries;

std::ostream& operator<<(std::ostream&, const PhotonSummary&);
I3_POINTER_TYPEDEFS(PhotonSummary);
I3_POINTER_TYPEDEFS(PhotonSummarySeries);
I3_DEFAULT_NAME(PhotonSummary);
#endif
