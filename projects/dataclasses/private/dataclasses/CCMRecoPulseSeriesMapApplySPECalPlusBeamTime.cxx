/**
 *  $Id$
 *
 *  Copyright (C) 2016
 *  Claudio Kopper <ckopper@icecube.wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#include "dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusBeamTime.h"
#include "dataclasses/physics/CCMRecoPulse.h"
#include "dataclasses/physics/CCMBCMSummary.h"
#include "dataclasses/geometry/CCMGeometry.h"
#include "dataclasses/physics/NIMLogicPulse.h"
#include "dataclasses/calibration/CCMCalibration.h"
#include "boost/make_shared.hpp"
#include "boost/foreach.hpp"

CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::CCMRecoPulseSeriesMapApplySPECalPlusBeamTime() 
    : pulses_key_(), calibration_key_(), shifted_() {}

CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::CCMRecoPulseSeriesMapApplySPECalPlusBeamTime(
        std::string const & pulses_key,
        std::string const & calibration_key,
        std::string const & nim_pulses_key,
        std::string const & geometry_key,
        std::string const & bcm_summary_key
        ) 
    : pulses_key_(pulses_key), calibration_key_(calibration_key), nim_pulses_key_(nim_pulses_key), geometry_key_(geometry_key), bcm_summary_key_(bcm_summary_key), shifted_() {}

CCMRecoPulseSeriesMapConstPtr
CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::Apply(const I3Frame &frame) const {
    typedef CCMRecoPulseSeriesMap Map;
    typedef boost::shared_ptr<const Map> MapConstPtr;
    typedef Map::value_type Pair;
    typedef Pair::second_type Series;
    typedef Series::value_type Element;

    if (shifted_)
        return shifted_;

    CCMBCMSummaryConstPtr bcm = frame.Get<CCMBCMSummaryConstPtr>(bcm_summary_key_);
    if(!bcm)
        log_fatal("Couldn't find '%s' in the frame!",
                bcm_summary_key_.c_str());

    CCMGeometryConstPtr geo = frame.Get<CCMGeometryConstPtr>(geometry_key_);
    if (!geo)
        log_fatal("Couldn't find '%s' in the frame!",
                geometry_key_.c_str());

    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame.Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses)
        log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());

    CCMCalibrationConstPtr calibration = frame.Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration)
        log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());

    MapConstPtr in_pulses = frame.Get<MapConstPtr>(pulses_key_);
    if (!in_pulses)
        log_fatal("Couldn't find '%s' in the frame!",
                pulses_key_.c_str());

    CCMRecoPulseSeriesMapPtr shifted = boost::make_shared<Map>();

    // Find the BCM nim pulse
    CCMPMTKey bcm_pmt_key = CCMPMTKey(10, 1);
    std::map<CCMPMTKey, CCMTriggerKey>::const_iterator bcm_trigger_key = 
        geo->trigger_copy_map.find(bcm_pmt_key);
    if(bcm_trigger_key == geo->trigger_copy_map.end())
        log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                geometry_key_.c_str());
    NIMLogicPulseSeriesMap::const_iterator bcm_trigger_nim_pulses =
        nim_pulses->find(bcm_trigger_key->second);
    if(bcm_trigger_nim_pulses == nim_pulses->end())
        log_fatal("Could not find PMT (%i/%u) in '%s'",
                bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                nim_pulses_key_.c_str());

    double bcm_nim_pulse_time = 0.0;
    double bcm_max_nim_pulse_length = 0.0;
    for(size_t j=0; j<bcm_trigger_nim_pulses->second.size(); ++j) {
        double length = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseLength();
        if(length > bcm_max_nim_pulse_length) {
            bcm_nim_pulse_time = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseTime();
            bcm_max_nim_pulse_length = length;
        }
    }

    for(Pair const & pair : *in_pulses) {
        if(pair.second.size() == 0) {
            shifted->insert(Pair(pair.first, Series()));
            continue;
        }
        // retrieve the calibration for this PMT
        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
            calibration->pmtCal.find(pair.first);
        if (calib == calibration->pmtCal.end()) 
            log_fatal("Could not find PMT (%i/%u) in '%s'",
                    pair.first.GetRegion(), pair.first.GetSensor(),
                    calibration_key_.c_str());
        CCMPMTCalibration const & pmt_calibration = calib->second;

        // retrieve the board trigger nim pulse time for this PMT
        std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key = 
            geo->trigger_copy_map.find(pair.first);
        if(trigger_key == geo->trigger_copy_map.end())
            log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                    pair.first.GetRegion(), pair.first.GetSensor(),
                    geometry_key_.c_str());
        NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
            nim_pulses->find(trigger_key->second);
        if(pmt_trigger_nim_pulses == nim_pulses->end())
            log_fatal("Could not find PMT (%i/%u) in '%s'",
                    pair.first.GetRegion(), pair.first.GetSensor(),
                    nim_pulses_key_.c_str());
        double nim_pulse_time = 0.0;
        double max_nim_pulse_length = 0.0;
        for(size_t j=0; j<pmt_trigger_nim_pulses->second.size(); ++j) {
            double length = pmt_trigger_nim_pulses->second.at(j).GetNIMPulseLength();
            if(length > max_nim_pulse_length) {
                nim_pulse_time = pmt_trigger_nim_pulses->second.at(j).GetNIMPulseTime();
                max_nim_pulse_length = length;
            }
        }

        // get the BCM time
        double bcm_time = bcm->bcm_start_time;

        // Load the SPE corrections
        double SPECorrection = 1.0;
        if((not std::isnan(pmt_calibration.GetMeanPMTCharge())) and pmt_calibration.GetMeanPMTCharge() > 0) {
            SPECorrection = 1.0 / pmt_calibration.GetMeanPMTCharge();
        }

        double time_correction = 0.0;
        if(not std::isnan(nim_pulse_time)) {
            time_correction -= nim_pulse_time;
        }
        if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
            time_correction -= pmt_calibration.GetPMTDeltaT();
        }
        if(not std::isnan(bcm_time)) {
            time_correction -= bcm_time;
        }
        if(not std::isnan(bcm_nim_pulse_time)) {
            time_correction -= (-bcm_nim_pulse_time);
        }

        // insert an entry for this DOM into the output map
        // and retrieve a reference
        Series &shiftedvec = (*shifted)[pair.first];

        for(Element const & element : pair.second) {
            // plain copy of the pulse first
            shiftedvec.push_back(element);

            // retrieve a reference to the new element
            Element &shifted_element = shiftedvec.back();

            // now set the shifted charge
            shifted_element.SetCharge(element.GetCharge() * SPECorrection);
            shifted_element.SetTime(element.GetTime() + time_correction);
        }
    }

    // save in cache
    shifted_ = shifted;

    return shifted_;
}

bool
CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::operator==(
        const CCMRecoPulseSeriesMapApplySPECalPlusBeamTime& other) const
{
    return (pulses_key_ == other.pulses_key_) && 
        (calibration_key_ == other.calibration_key_) &&
        (nim_pulses_key_ == other.nim_pulses_key_) &&
        (geometry_key_ == other.geometry_key_) &&
        (bcm_summary_key_ == other.bcm_summary_key_);
}

bool
CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::operator!=(
        const CCMRecoPulseSeriesMapApplySPECalPlusBeamTime& other) const
{
    return (pulses_key_ != other.pulses_key_) ||
        (calibration_key_ != other.calibration_key_) ||
        (nim_pulses_key_ != other.nim_pulses_key_) ||
        (geometry_key_ != other.geometry_key_) ||
        (bcm_summary_key_ != other.bcm_summary_key_);
}

std::ostream& CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::Print(std::ostream& os) const
{
    os << "Pulses: " << pulses_key_ << " Calibration: " << calibration_key_ << " NIMPulses: " << nim_pulses_key_ << " Geometry: " << geometry_key_ << " BCMSummary: " << bcm_summary_key_;
    return os;
}

std::ostream& operator<<(std::ostream& os, const CCMRecoPulseSeriesMapApplySPECalPlusBeamTime& corr)
{
    return(corr.Print(os));
}

template <class Archive>
    void
CCMRecoPulseSeriesMapApplySPECalPlusBeamTime::serialize(Archive& ar, unsigned version)
{
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PulsesKey", pulses_key_);
    ar & make_nvp("CalibrationKey", calibration_key_);
    ar & make_nvp("NIMPulsesKey", nim_pulses_key_);
    ar & make_nvp("GeometryKey", geometry_key_);
    ar & make_nvp("BCMSummaryKey", bcm_summary_key_);
}

I3_SERIALIZABLE(CCMRecoPulseSeriesMapApplySPECalPlusBeamTime);
