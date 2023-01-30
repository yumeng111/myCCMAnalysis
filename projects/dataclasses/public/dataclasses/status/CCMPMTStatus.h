#ifndef CCMPMTStatus_H_INCLUDED
#define CCMPMTStatus_H_INCLUDED

#include "dataclasses/Utility.h"
#include <icetray/CCMPMTKey.h>
#include <serialization/version.hpp>

static const unsigned ccmpmtstatus_version_ = 6;

class CCMPMTStatus {
public:
    CCMPMTStatus():
        hvPower(OnOff::Unknown),
        pmtHV(NAN)
{};

    ~CCMPMTStatus();

    bool operator==(const CCMPMTStatus& rhs) const {
        return (hvPower == rhs.hvPower && pmtHV == rhs.pmtHV);
    }

    bool operator!=(const CCMPMTStatus& rhs) const {
        return !operator==(rhs);
    }

    /**
     * There is also provision to turn on or off various settings in the 
     * PMT
     */
    enum OnOff {Unknown = -1, Off = 0, On = 1};

    OnOff hvPower;

    /**
     * Real operating PMT HV (Volts)
     */
    double pmtHV;

private:
    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive & ar, const unsigned version);
};

typedef std::map<CCMPMTKey, CCMPMTStatus> CCMPMTStatusMap;
I3_POINTER_TYPEDEFS(CCMPMTStatusMap);

I3_POINTER_TYPEDEFS(CCMPMTStatus);
I3_CLASS_VERSION(CCMPMTStatus, ccmpmtstatus_version_);

#endif //CCMPMTStatus_H_INCLUDED

