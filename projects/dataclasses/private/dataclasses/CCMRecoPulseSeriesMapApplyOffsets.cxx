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

inline double Mod(double x, double mod) {
    double result = std::fmod(x, mod);
    result += (result < 0) ? mod : 0;
    return result;
}

inline double SnapToMod(double sum, double original_offset, double target, double mod) {
    double sum_remainder = Mod(sum, mod);
    double original_remainder = Mod(original_offset, mod);
    // First snap the result to the original offset's modulo class
    // This ensures we remove any finer timing differences
    double snapped = sum - sum_remainder + original_remainder;
    snapped += (sum_remainder > original_remainder) ? 0 : -mod;

    // Now adjust the result to match the target's modulo class with minimal change
    double target_remainder = Mod(target, mod);
    double delta = target_remainder - original_remainder;
    if(std::abs(delta) > std::abs(delta - mod)) {
        delta = delta - mod;
    }
    snapped += delta;
    return snapped;
}

inline void SnapPulsesToMod(std::vector<CCMRecoPulse> & pulses, double shift, double original_offset, double target, double mod, bool zero_out_width=false) {
    // Pre-compute which direction we need to shift in order to match the target modulo class
    double snapped_remainder = Mod(original_offset + shift, mod);
    double target_remainder = Mod(target, mod);
    double delta = target_remainder - snapped_remainder;
    if(std::abs(delta) > std::abs(delta - mod)) {
        delta = delta - mod;
    }

    double original_remainder = Mod(original_offset, mod);

    for(CCMRecoPulse & pulse : pulses) {
        // First snap the pulse to the original offset's modulo class
        double sum = pulse.GetTime() + shift;
        double sum_remainder = Mod(sum, mod);
        double snapped = sum - sum_remainder + original_remainder;
        snapped += (sum_remainder > original_remainder) ? 0 : -mod;

        // Now adjust the pulse to match the target's modulo class
        snapped += delta;
        pulse.SetTime(snapped);
        pulse.SetWidth(pulse.GetWidth() * (1 - zero_out_width));
    }
}

CCMRecoPulseSeriesMapApplyOffsets::CCMRecoPulseSeriesMapApplyOffsets()
    : pulses_key_(), original_offsets_key_(), target_offsets_key_(), mod_(0.0), shifted_() {}

CCMRecoPulseSeriesMapApplyOffsets::CCMRecoPulseSeriesMapApplyOffsets(
        std::string const & pulses_key,
        std::string const & original_offsets_key,
        std::string const & target_offsets_key,
        double const & mod
        )
    : pulses_key_(pulses_key), original_offsets_key_(original_offsets_key), target_offsets_key_(target_offsets_key), mod_(mod), shifted_() {}

CCMRecoPulseSeriesMapConstPtr
CCMRecoPulseSeriesMapApplyOffsets::Apply(const I3Frame &frame) const {
    if (shifted_)
        return shifted_;

    I3MapPMTKeyDoubleConstPtr original_offsets = frame.Get<I3MapPMTKeyDoubleConstPtr>(original_offsets_key_);
    if (!original_offsets)
        log_fatal("Couldn't find '%s' in the frame!",
                original_offsets_key_.c_str());
    I3MapPMTKeyDoubleConstPtr target_offsets = frame.Get<I3MapPMTKeyDoubleConstPtr>(target_offsets_key_);
    if (!target_offsets)
        log_fatal("Couldn't find '%s' in the frame!",
                target_offsets_key_.c_str());

    CCMRecoPulseSeriesMapConstPtr in_pulses = frame.Get<CCMRecoPulseSeriesMapConstPtr>(pulses_key_);
    if (!in_pulses)
        log_fatal("Couldn't find '%s' in the frame!",
                pulses_key_.c_str());

    CCMRecoPulseSeriesMapPtr shifted = boost::make_shared<CCMRecoPulseSeriesMap>();

    if(mod_ == 0.0) {
        for(auto const & [key, series] : *in_pulses) {
            if(series.size() == 0) {
                shifted->insert({key, CCMRecoPulseSeries()});
                continue;
            }
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, double>::const_iterator original_offset =
                original_offsets->find(key);
            if (original_offset == original_offsets->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        key.GetRegion(), key.GetSensor(),
                        original_offsets_key_.c_str());

            std::map<CCMPMTKey, double>::const_iterator target_offset =
                target_offsets->find(key);
            if (target_offset == target_offsets->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        key.GetRegion(), key.GetSensor(),
                        target_offsets_key_.c_str());

            double const & time_correction = target_offset->second - original_offset->second;

            // insert an entry for this DOM into the output map
            // and retrieve a reference
            CCMRecoPulseSeries & shiftedvec = (*shifted)[key];

            for(CCMRecoPulse const & pulse : series) {
                // plain copy of the pulse first
                shiftedvec.push_back(pulse);

                // retrieve a reference to the new element
                CCMRecoPulse & shifted_pulse = shiftedvec.back();

                shifted_pulse.SetTime(pulse.GetTime() + time_correction);
            }
        }
    } else {
        for(auto const & [key, series] : *in_pulses) {
            if(series.size() == 0) {
                shifted->insert({key, CCMRecoPulseSeries()});
                continue;
            }
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, double>::const_iterator original_offset =
                original_offsets->find(key);
            if (original_offset == original_offsets->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        key.GetRegion(), key.GetSensor(),
                        original_offsets_key_.c_str());

            std::map<CCMPMTKey, double>::const_iterator target_offset =
                target_offsets->find(key);
            if (target_offset == target_offsets->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        key.GetRegion(), key.GetSensor(),
                        target_offsets_key_.c_str());

            // insert an entry for this DOM into the output map
            // and retrieve a reference
            CCMRecoPulseSeries & shifted_vec = (*shifted)[key];
            shifted_vec = series; // copy all pulses

            SnapPulsesToMod(shifted_vec, 0.0, original_offset->second, target_offset->second, mod_, false);
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
        (original_offsets_key_ == other.original_offsets_key_) &&
        (target_offsets_key_ == other.target_offsets_key_) &&
        (mod_ == other.mod_);
}

bool
CCMRecoPulseSeriesMapApplyOffsets::operator!=(
        const CCMRecoPulseSeriesMapApplyOffsets& other) const
{
    return (pulses_key_ != other.pulses_key_) ||
        (original_offsets_key_ != other.original_offsets_key_) ||
        (target_offsets_key_ != other.target_offsets_key_) ||
        (mod_ != other.mod_);
}

std::ostream& CCMRecoPulseSeriesMapApplyOffsets::Print(std::ostream& os) const
{
    os << "Pulses: " << pulses_key_ << " OriginalOffsets: " << original_offsets_key_
       << " TargetOffsets: " << target_offsets_key_ << " Mod: " << mod_;
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
    ar & make_nvp("OriginalOffsetsKey", original_offsets_key_);
    ar & make_nvp("TargetOffsetsKey", target_offsets_key_);
    ar & make_nvp("Mod", mod_);
}

I3_SERIALIZABLE(CCMRecoPulseSeriesMapApplyOffsets);
