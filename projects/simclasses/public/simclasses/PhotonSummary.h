#ifndef PhotonSummary_H_INCLUDED
#define PhotonSummary_H_INCLUDED

#include <utility>
#include <vector>
#include <tuple>
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

    float distance_uv = 0; // record distance as given by GetStepLength
    float original_wavelength = 0;
    float distance_visible = 0;
    float time = 0; // as as reported using g4 global time
    std::vector<size_t> n_photons_per_wls;
    WLSLocationSeries wls_loc;
    PhotonSource photon_source = PhotonSummary::PhotonSource::Unknown;
    PhotonSource current_process = PhotonSummary::PhotonSource::Unknown;

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

    PhotonSummary() = default;

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
