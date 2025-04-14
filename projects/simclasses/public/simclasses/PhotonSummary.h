#ifndef PhotonSummary_H_INCLUDED
#define PhotonSummary_H_INCLUDED

#include <utility>
#include <vector>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>
#include <dataclasses/I3Vector.h>
#include <simclasses/WLSLocation.h>

static const unsigned photonsummary_version_ = 3;

class PhotonSummary : public I3FrameObject {
public:
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Scintillation = 1,
        Cherenkov = 2,
        OpWLS = 3
    };

    float distance_uv; // record distance as given by GetStepLength
    float original_wavelength;
    float distance_visible;
    float time; // as as reported using g4 global time
    std::vector<size_t> n_photons_per_wls;
    WLSLocationSeries wls_loc;
    PhotonSource photon_source;
    PhotonSource current_process;

    SET_LOGGER("PhotonSummary");

    bool operator==(const PhotonSummary& rhs) const {
        return distance_uv == rhs.distance_uv
            && original_wavelength == rhs.original_wavelength
            && distance_visible == rhs.distance_visible
            && time == rhs.time
            && n_photons_per_wls == rhs.n_photons_per_wls
            && wls_loc == rhs.wls_loc
            && photon_source == rhs.photon_source
            && current_process == rhs.current_process;
    }


    PhotonSummary(float distance_uv_ = 0,
            float original_wavelength_ = 0,
            float distance_visible_ = 0,
            float time_ = 0,
            std::vector<size_t> n_photons_per_wls_ = {},
            WLSLocationSeries wls_loc_ = WLSLocationSeries(),
            PhotonSource photon_source_ = PhotonSummary::PhotonSource::Unknown,
            PhotonSource current_process_ = PhotonSummary::PhotonSource::Unknown)
        :
            distance_uv(distance_uv_),
            original_wavelength(original_wavelength_),
            distance_visible(distance_visible_),
            time(time_),
            n_photons_per_wls(n_photons_per_wls_),
            wls_loc(wls_loc_),
            photon_source(photon_source_),
            current_process(current_process_) {}

    std::ostream& Print(std::ostream&) const;

private:
    static const std::unordered_map<PhotonSummary::PhotonSource, std::string> photonSourceToProcessName;

    friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();

};

I3_CLASS_VERSION(PhotonSummary,photonsummary_version_);

typedef I3Vector<PhotonSummary> PhotonSummarySeries;

std::ostream& operator<<(std::ostream&, const PhotonSummary&);
I3_POINTER_TYPEDEFS(PhotonSummary);
I3_POINTER_TYPEDEFS(PhotonSummarySeries);
I3_DEFAULT_NAME(PhotonSummary);
#endif
