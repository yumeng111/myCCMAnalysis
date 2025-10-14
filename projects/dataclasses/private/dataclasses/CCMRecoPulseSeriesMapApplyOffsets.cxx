/**
 *  $Id$
 *
 *  Copyright (C) 2016
 *  Claudio Kopper <ckopper@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#include <string>
#include <vector>

#include "icetray/CCMPMTKey.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/CCMRecoPulseSeriesMapApplyOffsets.h"
#include "dataclasses/physics/CCMRecoPulse.h"

CCMRecoPulseSeriesMapApplyOffsets::CCMRecoPulseSeriesMapApplyOffsets()
    : pulses_key_(), offsets_key_(), shifted_() {}

CCMRecoPulseSeriesMapApplyOffsets::CCMRecoPulseSeriesMapApplyOffsets(
        std::string const & pulses_key,
        std::string const & offsets_key
        )
    : pulses_key_(pulses_key), offsets_key_(offsets_key), shifted_() {}

CCMRecoPulseSeriesMapConstPtr
CCMRecoPulseSeriesMapApplyOffsets::Apply(const I3Frame &frame) const {
    typedef CCMRecoPulseSeriesMap Map;
    typedef boost::shared_ptr<const Map> MapConstPtr;
    typedef Map::value_type Pair;
    typedef Pair::second_type Series;
    typedef Series::value_type Element;

    if (shifted_)
        return shifted_;

    I3MapPMTKeyDoubleConstPtr offsets = frame.Get<I3MapPMTKeyDoubleConstPtr>(offsets_key_);
    if (!offsets)
        log_fatal("Couldn't find '%s' in the frame!",
                offsets_key_.c_str());

    MapConstPtr in_pulses = frame.Get<MapConstPtr>(pulses_key_);
    if (!in_pulses)
        log_fatal("Couldn't find '%s' in the frame!",
                pulses_key_.c_str());

    CCMRecoPulseSeriesMapPtr shifted = boost::make_shared<Map>();

    for(Pair const & pair : *in_pulses) {
        if(pair.second.size() == 0) {
            shifted->insert(Pair(pair.first, Series()));
            continue;
        }
        // retrieve the calibration for this PMT
        std::map<CCMPMTKey, double>::const_iterator offset =
            offsets->find(pair.first);
        if (offset == offsets->end())
            log_fatal("Could not find PMT (%i/%u) in '%s'",
                    pair.first.GetRegion(), pair.first.GetSensor(),
                    offsets_key_.c_str());

        double const & time_correction = offset->second;

        // insert an entry for this DOM into the output map
        // and retrieve a reference
        Series &shiftedvec = (*shifted)[pair.first];

        for(Element const & element : pair.second) {
            // plain copy of the pulse first
            shiftedvec.push_back(element);

            // retrieve a reference to the new element
            Element &shifted_element = shiftedvec.back();

            shifted_element.SetTime(element.GetTime() + time_correction);
        }
    }

    // save in cache
    shifted_ = shifted;

    return shifted_;
}

bool
CCMRecoPulseSeriesMapApplyOffsets::operator==(
        const CCMRecoPulseSeriesMapApplyOffsets& other) const
{
    return (pulses_key_ == other.pulses_key_) &&
        (offsets_key_ == other.offsets_key_);
}

bool
CCMRecoPulseSeriesMapApplyOffsets::operator!=(
        const CCMRecoPulseSeriesMapApplyOffsets& other) const
{
    return (pulses_key_ != other.pulses_key_) ||
        (offsets_key_ != other.offsets_key_);
}

std::ostream& CCMRecoPulseSeriesMapApplyOffsets::Print(std::ostream& os) const
{
    os << "Pulses: " << pulses_key_ << " Offsets: " << offsets_key_;
    return os;
}

std::ostream& operator<<(std::ostream& os, const CCMRecoPulseSeriesMapApplyOffsets& corr)
{
    return(corr.Print(os));
}

template <class Archive>
    void
CCMRecoPulseSeriesMapApplyOffsets::serialize(Archive& ar, unsigned version)
{
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PulsesKey", pulses_key_);
    ar & make_nvp("OffsetsKey", offsets_key_);
}

I3_SERIALIZABLE(CCMRecoPulseSeriesMapApplyOffsets);
