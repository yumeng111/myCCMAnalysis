#include <icetray/serialization.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <icetray/I3Units.h>
#include <dataclasses/external/CompareFloatingPoint.h>

NIMLogicPulse::~NIMLogicPulse() {}

template <class Archive> void NIMLogicPulse::serialize(Archive& ar, unsigned version) {
    if (version>nimlogicpulse_version_)
        log_fatal("Attempting to read version %u from file but running version %u of NIMLogicPulse class.", version, nimlogicpulse_version_);

    ar & make_nvp("NIMPulseTime", time_);
    ar & make_nvp("NIMPulseLength", length_);
}

using CompareFloatingPoint::Compare;

bool 
NIMLogicPulse::operator==(const NIMLogicPulse& rhs) const {
    return Compare(time_, rhs.time_) && Compare(length_, rhs.length_);
}

bool 
NIMLogicPulse::operator!=(const NIMLogicPulse& rhs) const {
    return !( *this == rhs );
}

std::ostream& NIMLogicPulse::Print(std::ostream& oss) const {
    oss << " length = " << GetNIMPulseLength() << " ns"
        << " time = " << GetNIMPulseTime() << " ns";
    return oss;
}

std::ostream& operator<<(std::ostream& oss, const NIMLogicPulse& t) {
    return(t.Print(oss));
}

I3_SERIALIZABLE(NIMLogicPulse);
I3_SERIALIZABLE(I3VectorNIMLogicPulse);
