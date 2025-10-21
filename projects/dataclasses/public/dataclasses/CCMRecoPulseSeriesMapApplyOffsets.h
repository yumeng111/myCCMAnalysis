/**
 *  $Id$
 *
 *  Copyright (C) 2016
 *  Claudio Kopper <ckopper@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#ifndef dataclasses_CCMRecoPulseSeriesMapApplyOffsets_H
#define dataclasses_CCMRecoPulseSeriesMapApplyOffsets_H

#include <string>
#include "icetray/I3FrameObject.h"
#include "icetray/I3Frame.h"
#include "icetray/serialization.h"
#include "dataclasses/physics/CCMRecoPulse.h"

static const unsigned ccmrecopulseseriesmapapplyoffsets_version_ = 0;

class CCMRecoPulseSeriesMapApplyOffsets : public I3FrameObject {
public:
    /*
     * Apply per-PMT time offsets to a CCMRecoPulseSeriesMap.
     */
    CCMRecoPulseSeriesMapApplyOffsets(
            std::string const & pulses_key,
            std::string const & original_offsets_key,
            std::string const & target_offsets_key,
            double const & mod = 0.0
    );
    CCMRecoPulseSeriesMapApplyOffsets();

    std::ostream& Print(std::ostream&) const override;

    CCMRecoPulseSeriesMapConstPtr Apply(const I3Frame&) const;
    std::string GetPulsesSource() const { return pulses_key_; }
    std::string GetOriginalOffsetsSource() const { return original_offsets_key_; }
    std::string GetTargetOffsetsSource() const { return target_offsets_key_; }
    double GetMod() const { return mod_; }

    bool operator==(const CCMRecoPulseSeriesMapApplyOffsets&) const;
    bool operator!=(const CCMRecoPulseSeriesMapApplyOffsets&) const;
private:
    std::string pulses_key_;
    std::string original_offsets_key_;
    std::string target_offsets_key_;
    double mod_;
    mutable CCMRecoPulseSeriesMapPtr shifted_;

    friend class icecube::serialization::access;
    template <class Archive> void serialize(Archive& ar, unsigned version);

    SET_LOGGER("CCMRecoPulseSeriesMapApplyOffsets");
};

std::ostream& operator<<(std::ostream&, const CCMRecoPulseSeriesMapApplyOffsets&);

I3_CLASS_VERSION(CCMRecoPulseSeriesMapApplyOffsets,
        ccmrecopulseseriesmapapplyoffsets_version_);
I3_POINTER_TYPEDEFS(CCMRecoPulseSeriesMapApplyOffsets);

#endif // dataclasses_CCMRecoPulseSeriesMapApplyOffsets_H
