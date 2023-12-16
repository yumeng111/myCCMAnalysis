#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <string>
#include <random>
#include <iostream>
#include <algorithm>

#include <icetray/I3Frame.h>
#include <icetray/I3ConditionalModule.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/calibration/CCMCalibration.h>

#include "daqtools/OnlineRobustStats.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"

class AccumulateIndividualPMTWaveforms : public I3ConditionalModule {
    std::string geometry_key_;
    std::string nim_pulses_key_;
    std::string bcm_summary_key_;
    std::string calibration_key_;
    std::string waveforms_key_;
    std::string output_prefix_;
    bool consume_frames_;
    bool correct_nim_pulse_time_;
    bool correct_electron_transit_time_;
    I3Frame::Stream output_frame_type_;

    bool geo_seen_ = false;
    CCMGeometry geo_;
    I3Vector<CCMPMTKey> allowed_pmt_keys_;
    std::vector<CCMPMTKey> pmt_keys_;
    std::map<CCMPMTKey, WaveformAccumulator> accumulated_waveforms_;
    bool trigger_reference_time_;
    bool bcm_reference_time_;
    std::string reference_time_key_;

    enum class ReferenceTimeType {TriggerTime, BCMTime, UserSpecifiedTime};
    ReferenceTimeType chosen_time_reference_;

public:
    AccumulateIndividualPMTWaveforms(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    void DumpWaveforms();

    void ProcessFrame(I3FramePtr frame);
    std::map<CCMPMTKey, int32_t> GetReferenceIndices(I3FramePtr frame);
    std::map<CCMPMTKey, int32_t> GetReferenceIndicesTrigger(I3FramePtr frame);
    std::map<CCMPMTKey, int32_t> GetReferenceIndicesBCM(I3FramePtr frame);
    std::map<CCMPMTKey, int32_t> GetReferenceIndicesUser(I3FramePtr frame);
};

I3_MODULE(AccumulateIndividualPMTWaveforms);

std::map<CCMPMTKey, int32_t> AccumulateIndividualPMTWaveforms::GetReferenceIndicesTrigger(I3FramePtr frame) {
    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame->Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses and correct_nim_pulse_time_)
        log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());

    CCMCalibrationConstPtr calibration = frame->Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration and correct_electron_transit_time_)
        log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());

    std::map<CCMPMTKey, int32_t> output_indices;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        double time_correction = 0.0;
        CCMOMGeo const & pmt = geo_.pmt_geo.at(pmt_key);
        CCMPMTType const & pmt_type = pmt.omtype;
        bool is_pmt = pmt_type == CCMPMTType::CCM8inCoated or pmt_type == CCMPMTType::CCM8inUncoated or pmt_type == CCMPMTType::CCM1in;

        if(correct_electron_transit_time_ and is_pmt) {
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
                calibration->pmtCal.find(pmt_key);
            if (calib == calibration->pmtCal.end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        calibration_key_.c_str());
            CCMPMTCalibration const & pmt_calibration = calib->second;
            if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
                time_correction -= pmt_calibration.GetPMTDeltaT();
            }
        }

        if(correct_nim_pulse_time_) {
            // retrieve the board trigger nim pulse time for this PMT
            std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key =
                geo_.trigger_copy_map.find(pmt_key);
            if(trigger_key == geo_.trigger_copy_map.end())
                log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        geometry_key_.c_str());
            NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
                nim_pulses->find(trigger_key->second);
            if(pmt_trigger_nim_pulses == nim_pulses->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
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
            if(not std::isnan(nim_pulse_time)) {
                time_correction -= nim_pulse_time;
            }
        }
        output_indices[pmt_key] = (int32_t)(-time_correction / 2.0);
    }

    return output_indices;
}

std::map<CCMPMTKey, int32_t> AccumulateIndividualPMTWaveforms::GetReferenceIndicesBCM(I3FramePtr frame) {
    CCMBCMSummaryConstPtr bcm = frame->Get<CCMBCMSummaryConstPtr>(bcm_summary_key_);
    if(!bcm)
        log_fatal("Couldn't find '%s' in the frame!",
                bcm_summary_key_.c_str());

    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame->Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses and correct_nim_pulse_time_)
        log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());

    CCMCalibrationConstPtr calibration = frame->Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration and correct_electron_transit_time_)
        log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());

    // Find the BCM nim pulse
    CCMPMTKey bcm_pmt_key = CCMPMTKey(10, 1);
    std::map<CCMPMTKey, CCMTriggerKey>::const_iterator bcm_trigger_key =
        geo_.trigger_copy_map.find(bcm_pmt_key);
    if(bcm_trigger_key == geo_.trigger_copy_map.end())
        log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                geometry_key_.c_str());

    NIMLogicPulseSeriesMap::const_iterator bcm_trigger_nim_pulses;
    if(correct_nim_pulse_time_) {
        bcm_trigger_nim_pulses = nim_pulses->find(bcm_trigger_key->second);
        if(bcm_trigger_nim_pulses == nim_pulses->end())
            log_fatal("Could not find PMT (%i/%u) in '%s'",
                    bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                    nim_pulses_key_.c_str());
    }

    double bcm_nim_pulse_time = 0.0;
    double bcm_max_nim_pulse_length = 0.0;
    if(correct_nim_pulse_time_) {
        for(size_t j=0; j<bcm_trigger_nim_pulses->second.size(); ++j) {
            double length = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseLength();
            if(length > bcm_max_nim_pulse_length) {
                bcm_nim_pulse_time = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseTime();
                bcm_max_nim_pulse_length = length;
            }
        }
    }

    // get the BCM time
    double bcm_time = bcm->bcm_start_time;

    std::map<CCMPMTKey, int32_t> output_indices;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        double time_correction = 0.0;

        CCMOMGeo const & pmt = geo_.pmt_geo.at(pmt_key);
        CCMPMTType const & pmt_type = pmt.omtype;
        bool is_pmt = pmt_type == CCMPMTType::CCM8inCoated or pmt_type == CCMPMTType::CCM8inUncoated or pmt_type == CCMPMTType::CCM1in;

        if(correct_electron_transit_time_ and is_pmt) {
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
                calibration->pmtCal.find(pmt_key);
            if (calib == calibration->pmtCal.end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        calibration_key_.c_str());
            CCMPMTCalibration const & pmt_calibration = calib->second;
            if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
                time_correction -= pmt_calibration.GetPMTDeltaT();
            }
        }

        if(correct_nim_pulse_time_) {
            // retrieve the board trigger nim pulse time for this PMT
            std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key =
                geo_.trigger_copy_map.find(pmt_key);
            if(trigger_key == geo_.trigger_copy_map.end())
                log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        geometry_key_.c_str());
            NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
                nim_pulses->find(trigger_key->second);
            if(pmt_trigger_nim_pulses == nim_pulses->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
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
            if(not std::isnan(nim_pulse_time)) {
                time_correction -= nim_pulse_time;
            }
            if(not std::isnan(bcm_nim_pulse_time)) {
                time_correction -= (-bcm_nim_pulse_time);
            }
        }

        if(not std::isnan(bcm_time)) {
            time_correction -= bcm_time;
        }

        output_indices[pmt_key] = (int32_t)(-time_correction / 2.0);
    }

    return output_indices;
}

std::map<CCMPMTKey, int32_t> AccumulateIndividualPMTWaveforms::GetReferenceIndicesUser(I3FramePtr frame) {
    boost::shared_ptr<I3Map<CCMPMTKey, int32_t> const> reference_times = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, int32_t> const>>(reference_time_key_);
    if(reference_times == nullptr) {
        log_fatal("Couldn't find '%s' in the frame!",
                reference_time_key_.c_str());
    }

    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame->Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses and correct_nim_pulse_time_)
        log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());

    CCMCalibrationConstPtr calibration = frame->Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration and correct_electron_transit_time_)
        log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());

    std::map<CCMPMTKey, int32_t> output_indices;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        double time_correction = 0.0;

        CCMOMGeo const & pmt = geo_.pmt_geo.at(pmt_key);
        CCMPMTType const & pmt_type = pmt.omtype;
        bool is_pmt = pmt_type == CCMPMTType::CCM8inCoated or pmt_type == CCMPMTType::CCM8inUncoated or pmt_type == CCMPMTType::CCM1in;

        if(correct_electron_transit_time_ and is_pmt) {
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
                calibration->pmtCal.find(pmt_key);
            if (calib == calibration->pmtCal.end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        calibration_key_.c_str());
            CCMPMTCalibration const & pmt_calibration = calib->second;
            if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
                time_correction -= pmt_calibration.GetPMTDeltaT();
            }
        }

        if(correct_nim_pulse_time_) {
            // retrieve the board trigger nim pulse time for this PMT
            std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key =
                geo_.trigger_copy_map.find(pmt_key);
            if(trigger_key == geo_.trigger_copy_map.end())
                log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        geometry_key_.c_str());
            NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
                nim_pulses->find(trigger_key->second);
            if(pmt_trigger_nim_pulses == nim_pulses->end())
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
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
            if(not std::isnan(nim_pulse_time)) {
                time_correction -= nim_pulse_time;
            }
        }

        if(reference_times and reference_times->count(pmt_key)) {
            time_correction -= (*reference_times).at(pmt_key);
        }

        output_indices[pmt_key] = (int32_t)(-time_correction / 2.0);
    }

    return output_indices;
}

std::map<CCMPMTKey, int32_t> AccumulateIndividualPMTWaveforms::GetReferenceIndices(I3FramePtr frame) {
    switch(chosen_time_reference_) {
        case ReferenceTimeType::TriggerTime:
            return GetReferenceIndicesTrigger(frame);
        case ReferenceTimeType::BCMTime:
            return GetReferenceIndicesBCM(frame);
        case ReferenceTimeType::UserSpecifiedTime:
            return GetReferenceIndicesUser(frame);
        default:
            break;
    };
    log_fatal("Received enum type that is not initialized properly!");
    return std::map<CCMPMTKey, int32_t>();
}

void AccumulateIndividualPMTWaveforms::ProcessFrame(I3FramePtr frame) {
    std::map<CCMPMTKey, int32_t> reference_indices = GetReferenceIndices(frame);
    boost::shared_ptr<CCMWaveformUInt16Series const> waveform_raw = frame->Get<boost::shared_ptr<CCMWaveformUInt16Series const>>(waveforms_key_);
    boost::shared_ptr<CCMWaveformDoubleSeries const> waveform_cal = frame->Get<boost::shared_ptr<CCMWaveformDoubleSeries const>>(waveforms_key_);
    if(!waveform_raw and !waveform_cal)
        log_fatal("Couldn't find '%s' in the frame!",
                waveforms_key_.c_str());

    I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map = geo_.pmt_channel_map;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        if(pmt_channel_map.count(pmt_key) == 0) {
            log_fatal("Could not find PMT (%i/%u) in '%s'",
                    pmt_key.GetRegion(), pmt_key.GetSensor(),
                    geometry_key_.c_str());
        }
        uint32_t channel = pmt_channel_map.at(pmt_key);
        if(waveform_raw) {
            if(channel >= waveform_raw->size()) {
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        waveforms_key_.c_str());
            }
            CCMWaveformUInt16 const & waveform = waveform_raw->at(channel);
            std::vector<double> wf(waveform.GetWaveform().begin(), waveform.GetWaveform().end());
            accumulated_waveforms_[pmt_key].AddWaveform(wf, reference_indices[pmt_key]);
        } else if(waveform_cal) {
            if(channel >= waveform_cal->size()) {
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        waveforms_key_.c_str());
            }
            CCMWaveformDouble const & waveform = waveform_cal->at(channel);
            accumulated_waveforms_[pmt_key].AddWaveform(waveform.GetWaveform(), reference_indices[pmt_key]);
        }
    }
}

AccumulateIndividualPMTWaveforms::AccumulateIndividualPMTWaveforms(const I3Context& context) : I3ConditionalModule(context) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMCalibrationName", "Key for CCMCalibration", std::string(I3DefaultName<CCMCalibration>::value()));
    AddParameter("NIMPulsesName", "Key for NIMPulses", std::string("NIMPulses"));
    AddParameter("BCMSummaryName", "Key for BCMSummary", std::string("BCMSummary"));
    AddParameter("WaveformsKey", "Key for Waveforms", std::string("CCMWaveforms"));
    AddParameter("OutputPrefix", "Prefix for the module output", std::string("AccumulatedWaveforms"));
    AddParameter("ConsumeFrames", "Consume frames used as input?", bool(true));
    AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>());
    AddParameter("TriggerReferenceTime", "Use the trigger time as a reference time? This is the default", bool(false));
    AddParameter("BCMReferenceTime", "Use the Beam Current Monitor start time as a reference time?", bool(false));
    AddParameter("ReferenceTimeKey", "Name of a frame key that contains an I3Map<CCMPMTKey, int32_t> for the PMT reference indices", std::string(""));
    AddParameter("CorrectNIMPulseTime", "Correct for channel 15 NIM pulse arrival time?", bool(true));
    AddParameter("CorrectElectronTransitTime", "Correct for PMT electron transit time?", bool(true));
    AddParameter("OutputFrameType", "The type of frame to use in the ouptut. Default: DAQ", I3Frame::DAQ);
}

void AccumulateIndividualPMTWaveforms::Configure() {
    GetParameter("CCMGeometryName", geometry_key_);
    GetParameter("CCMCalibrationName", calibration_key_);
    GetParameter("NIMPulsesName", nim_pulses_key_);
    GetParameter("BCMSummaryName", bcm_summary_key_);
    GetParameter("WaveformsKey", waveforms_key_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("ConsumeFrames", consume_frames_);
    GetParameter("PMTKeys", allowed_pmt_keys_);
    GetParameter("TriggerReferenceTime", trigger_reference_time_);
    GetParameter("BCMReferenceTime", bcm_reference_time_);
    GetParameter("ReferenceTimeKey", reference_time_key_);
    GetParameter("CorrectNIMPulseTime", correct_nim_pulse_time_);
    GetParameter("CorrectElectronTransitTime", correct_electron_transit_time_);
    GetParameter("OutputFrameType", output_frame_type_);

    // Preemtively sort the allowed_pmt_keys_ so they're ready for the Geometry function
    if(allowed_pmt_keys_.size() > 0)
        std::sort(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end());

    unsigned int n_options = 0;
    n_options += (unsigned int)(trigger_reference_time_);
    n_options += (unsigned int)(bcm_reference_time_);
    n_options += (unsigned int)(reference_time_key_ != "");

    if(n_options == 0) {
        chosen_time_reference_ = ReferenceTimeType::TriggerTime;
    } else if(n_options == 1) {
        if(trigger_reference_time_) {
            chosen_time_reference_ = ReferenceTimeType::TriggerTime;
        } else if (bcm_reference_time_) {
            chosen_time_reference_ = ReferenceTimeType::BCMTime;
        } else {
            chosen_time_reference_ = ReferenceTimeType::UserSpecifiedTime;
        }
    } else {
        std::stringstream ss;
        ss << "More than one reference time option specified! Please select only one. Received ";
        ss << "TriggerReferenceTime == " << (trigger_reference_time_ ? "true" : "false") << ", ";
        ss << "BCMReferenceTime == " << (bcm_reference_time_ ? "true" : "false") << ", ";
        ss << "ReferenceTimeKey == \"" << reference_time_key_ << "\"";
        log_fatal(ss.str().c_str());
    }
}

void AccumulateIndividualPMTWaveforms::Geometry(I3FramePtr frame) {
    // Dump the accumulated waveforms if we get a new geometry
    if(geo_seen_) {
        DumpWaveforms();
    }

    // Assumes allowed_pmt_keys_ is already sorted
    CCMGeometryConstPtr geo = frame->Get<CCMGeometryConstPtr>(geometry_key_);
    if (!geo)
        log_fatal("Couldn't find '%s' in the frame!",
                geometry_key_.c_str());
    geo_ = *geo;
    geo_seen_ = true;

    // Copy and sort the pmt keys from the geometry
    std::vector<CCMPMTKey> geo_keys; geo_keys.reserve(geo_.pmt_channel_map.size());
    for(std::pair<CCMPMTKey const, uint32_t> const & p : geo_.pmt_channel_map)
        geo_keys.push_back(p.first);
    std::sort(geo_keys.begin(), geo_keys.end());

    // Assume we're doing all PMTs if none are specified
    if(allowed_pmt_keys_.size() == 0) {
        pmt_keys_ = geo_keys;
        PushFrame(frame);
        return;
    }

    // Clear the final pmt key list
    pmt_keys_.clear();
    // Fill the final pmt key list with the intersection of specified pmt keys and those available in the geometry
    // If allowed_pmt_keys_ is not sorted then this will fail horribly
    std::set_intersection(geo_keys.begin(), geo_keys.end(), allowed_pmt_keys_.begin(), allowed_pmt_keys_.end(), std::back_inserter(pmt_keys_));

    if(pmt_keys_.size() < allowed_pmt_keys_.size()) {
        std::vector<CCMPMTKey> missing_keys;
        std::set_difference(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end(), pmt_keys_.begin(), pmt_keys_.end(), std::back_inserter(missing_keys));
        std::stringstream ss;
        ss << "Some specified CCMPMTKeys are not present in the geometry:";
        for(CCMPMTKey const & key : missing_keys)
            ss << " " << key;
        log_warn(ss.str().c_str());
    }

    // Initialize the accumulated waveforms
    accumulated_waveforms_.clear();
    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        if(accumulated_waveforms_.count(pmt_key) == 0) {
            accumulated_waveforms_[pmt_key] = WaveformAccumulator();
        }
    }

    PushFrame(frame);
}

void AccumulateIndividualPMTWaveforms::DAQ(I3FramePtr frame) {
    ProcessFrame(frame);
    if(not consume_frames_)
        PushFrame(frame);
}

void AccumulateIndividualPMTWaveforms::DumpWaveforms() {
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> p_samples = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<uint32_t>>> p_counts = boost::make_shared<I3Map<CCMPMTKey, std::vector<uint32_t>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, int32_t>> p_fixed_positions = boost::make_shared<I3Map<CCMPMTKey, int32_t>>();

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        std::deque<double> deque_samples = accumulated_waveforms_[pmt_key].GetSummedWaveform();
        std::deque<unsigned int> deque_counts = accumulated_waveforms_[pmt_key].GetCounts();
        int32_t fixed_position = accumulated_waveforms_[pmt_key].GetFixedPosition();

        p_samples->insert(std::make_pair<CCMPMTKey, std::vector<double>>(CCMPMTKey(pmt_key), std::vector<double>(deque_samples.begin(), deque_samples.end())));
        p_counts->insert(std::make_pair<CCMPMTKey, std::vector<uint32_t>>(CCMPMTKey(pmt_key), std::vector<uint32_t>(deque_counts.begin(), deque_counts.end())));
        p_fixed_positions->insert(std::make_pair<CCMPMTKey, int32_t>(CCMPMTKey(pmt_key), std::move(fixed_position)));
    }

    I3FramePtr frame = boost::make_shared<I3Frame>(output_frame_type_);

    frame->Put((output_prefix_).c_str(), p_samples);
    frame->Put((output_prefix_ + "Counts").c_str(), p_counts);
    frame->Put((output_prefix_ + "FixedPosition").c_str(), p_fixed_positions);
    PushFrame(frame);
}

void AccumulateIndividualPMTWaveforms::Finish() {
    DumpWaveforms();
    Flush();
}
