#include <icetray/serialization.h>
#include <dataclasses/status/CCMDetectorStatus.h>

CCMDetectorStatus::~CCMDetectorStatus() {}

template <class Archive>
void
CCMDetectorStatus::save(Archive& ar, unsigned version) const {
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTStatus", pmtStatus);
    ar & make_nvp("StartTime", startTime);
    ar & make_nvp("EndTime", endTime);
}

template <class Archive>
void
CCMDetectorStatus::load(Archive& ar, unsigned version) {
    if (version>ccmdetectorstatus_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMDetectorStatus class.",
                version,ccmdetectorstatus_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTStatus", pmtStatus);
    ar & make_nvp("StartTime", startTime);
    ar & make_nvp("EndTime", endTime);
}

I3_SPLIT_SERIALIZABLE(CCMDetectorStatus);
