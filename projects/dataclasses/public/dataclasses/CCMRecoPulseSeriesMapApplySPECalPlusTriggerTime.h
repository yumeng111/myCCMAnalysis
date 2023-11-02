/**
 *  $Id$
 *
 *  Copyright (C) 2016
 *  Claudio Kopper <ckopper@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#ifndef dataclasses_CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime_H
#define dataclasses_CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime_H

#include <functional>
#include <string>
#include "icetray/I3FrameObject.h"
#include "icetray/OMKey.h"
#include "icetray/I3Frame.h"
#include "icetray/serialization.h"
#include "dataclasses/physics/CCMRecoPulse.h"

static const unsigned ccmrecopulseseriesmapapplyspecalplustriggertime_version_ = 0;

class CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime : public I3FrameObject {
public:
    /*
     * Use the SPE correction from CCMCalibration to shift the pulse amplitudes.
     * Assumes the input pulses are not shifted yet.
     * Use the NIM pulses to shift the pulse times to be relative to the trigger time.
     */
    CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime(
            std::string const & pulses_key,
            std::string const & calibration_key,
            std::string const & nim_pulses_key,
            std::string const & geometry_key
    );
    CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime();

    std::ostream& Print(std::ostream&) const override;

    CCMRecoPulseSeriesMapConstPtr Apply(const I3Frame&) const;
    std::string GetPulsesSource() const { return pulses_key_; }
    std::string GetCalibrationSource() const { return calibration_key_; }
    std::string GetNIMPulsesSource() const { return nim_pulses_key_; }
    std::string GetGeometrySource() const { return geometry_key_; }

    bool operator==(const CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime&) const;
    bool operator!=(const CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime&) const;
private:
    std::string pulses_key_;
    std::string calibration_key_;
    std::string nim_pulses_key_;
    std::string geometry_key_;
    mutable CCMRecoPulseSeriesMapPtr shifted_;

    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive& ar, unsigned version);

    SET_LOGGER("CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime");
};

std::ostream& operator<<(std::ostream&, const CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime&);

I3_CLASS_VERSION(CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime,
        ccmrecopulseseriesmapapplyspecalplustriggertime_version_);
I3_POINTER_TYPEDEFS(CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime);

#endif // dataclasses_CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime_H
