#include <icetray/serialization.h>
#include <dataclasses/status/CCMPMTStatus.h>
#include "icetray/I3Units.h"

CCMPMTStatus::~CCMPMTStatus() {}

template <class Archive>
void CCMPMTStatus::serialize (Archive& ar, const unsigned version) {
    if (version>ccmpmtstatus_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMPMTStatus class.",version,ccmpmtstatus_version_);

    ar & make_nvp("hvPower", hvPower);
    ar & make_nvp("PMTHV", pmtHV);
}

I3_SERIALIZABLE(CCMPMTStatus);
