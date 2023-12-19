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
#include <dataclasses/calibration/BaselineEstimate.h>
#include <CCMAnalysis/CCMBinary/BinaryFormat.h>

#include "daqtools/OnlineRobustStats.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"

class AccumulateIndividualChannelWaveforms : public I3ConditionalModule {
    struct Corrections {
        bool nim_pulse_time = false;
        bool bcm_nim_pulse_time = false;
        bool electron_transit_time = false;
        bool bcm_start_time = false;
        bool user_time = false;
        bool baseline = false;
        bool operator==(Corrections const & o) const {
            return std::tie(
                    nim_pulse_time,
                    bcm_nim_pulse_time,
                    electron_transit_time,
                    bcm_start_time,
                    user_time,
                    baseline
                    ) == std::tie(
                    o.nim_pulse_time,
                    o.bcm_nim_pulse_time,
                    o.electron_transit_time,
                    o.bcm_start_time,
                    o.user_time,
                    o.baseline
                    );
        }
        bool operator!=(Corrections const & o) const {
            return not (*this == o);
        }
    };

    std::string geometry_key_;
    std::string daq_config_key_;
    std::string nim_pulses_key_;
    std::string bcm_summary_key_;
    std::string calibration_key_;
    std::string waveforms_key_;
    std::string output_prefix_;
    std::string baseline_estimates_key_;
    bool consume_frames_;
    bool correct_nim_pulse_time_;
    bool correct_electron_transit_time_;
    bool correct_baseline_and_invert_raw_waveforms_;
    bool allow_missing_;
    bool skip_missing_;
    I3Frame::Stream output_frame_type_;

    bool geo_seen_ = false;
    CCMGeometry geo_;
    I3Vector<uint32_t> allowed_channels_;
    std::vector<uint32_t> channels_;
    std::map<uint32_t, WaveformAccumulator> accumulated_waveforms_;
    bool trigger_reference_time_;
    bool bcm_reference_time_;
    std::string reference_time_key_;

    std::map<uint32_t, Corrections> corrections_;
    std::map<uint32_t, CCMPMTKey> channel_pmt_map_;
    std::map<uint32_t, CCMTriggerKey> channel_trigger_map_;

    enum class ReferenceTimeType {TriggerTime, BCMTime, UserSpecifiedTime};
    ReferenceTimeType chosen_time_reference_;

public:
    AccumulateIndividualChannelWaveforms(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    void DumpWaveforms();

    void ProcessFrame(I3FramePtr frame);
    bool ComputeReferenceIndices(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections);
    bool ComputeReferenceIndicesTrigger(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections);
    bool ComputeReferenceIndicesBCM(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections);
    bool ComputeReferenceIndicesUser(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections);
};

I3_MODULE(AccumulateIndividualChannelWaveforms);

bool AccumulateIndividualChannelWaveforms::ComputeReferenceIndicesTrigger(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections) {
    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame->Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses and correct_nim_pulse_time_) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                    nim_pulses_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());
    }

    CCMCalibrationConstPtr calibration = frame->Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration and correct_electron_transit_time_) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                calibration_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());
    }

    for(uint32_t const & channel : channels_) {
        Corrections pmt_corrections;

        double time_correction = 0.0;
        bool is_pmt = false;
        CCMPMTKey pmt_key;
        if(channel_pmt_map_.count(channel)) {
            pmt_key = channel_pmt_map_.at(channel);
            CCMOMGeo const & pmt = geo_.pmt_geo.at(pmt_key);
            CCMPMTType const & pmt_type = pmt.omtype;
            is_pmt = pmt_type == CCMPMTType::CCM8inCoated or pmt_type == CCMPMTType::CCM8inUncoated or pmt_type == CCMPMTType::CCM1in;
        }

        if(correct_electron_transit_time_ and is_pmt) {
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
                calibration->pmtCal.find(pmt_key);
            if (calib == calibration->pmtCal.end()) {
                if(not allow_missing_) {
                    log_fatal("Could not find PMT (%i/%u) in '%s'",
                            pmt_key.GetRegion(), pmt_key.GetSensor(),
                            calibration_key_.c_str());
                }
            } else {
                CCMPMTCalibration const & pmt_calibration = calib->second;
                if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
                    time_correction -= pmt_calibration.GetPMTDeltaT();
                    pmt_corrections.electron_transit_time = true;
                }
            }
        }

        if(correct_nim_pulse_time_) {
            // retrieve the board trigger nim pulse time for this PMT
            std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key =
                geo_.trigger_copy_map.find(pmt_key);
            if(trigger_key == geo_.trigger_copy_map.end()) {
                if(not allow_missing_) {
                    log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                            pmt_key.GetRegion(), pmt_key.GetSensor(),
                            geometry_key_.c_str());
                }
            } else {
                NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
                    nim_pulses->find(trigger_key->second);
                if(pmt_trigger_nim_pulses == nim_pulses->end()) {
                    if(not allow_missing_) {
                        log_fatal("Could not find PMT (%i/%u) in '%s'",
                                pmt_key.GetRegion(), pmt_key.GetSensor(),
                                nim_pulses_key_.c_str());
                    }
                } else {
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
                        pmt_corrections.nim_pulse_time = true;
                    }
                }
            }
        }
        output_indices[pmt_key] = (int32_t)(-time_correction / 2.0);
        output_corrections[pmt_key] = pmt_corrections;
    }

    return true;
}

bool AccumulateIndividualChannelWaveforms::ComputeReferenceIndicesBCM(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections) {
    CCMBCMSummaryConstPtr bcm = frame->Get<CCMBCMSummaryConstPtr>(bcm_summary_key_);
    if(!bcm) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                    bcm_summary_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                bcm_summary_key_.c_str());
    }

    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame->Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses and correct_nim_pulse_time_) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                    nim_pulses_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());
    }

    CCMCalibrationConstPtr calibration = frame->Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration and correct_electron_transit_time_) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                calibration_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());
    }

    // Find the BCM nim pulse
    CCMPMTKey bcm_pmt_key = CCMPMTKey(10, 1);
    std::map<CCMPMTKey, CCMTriggerKey>::const_iterator bcm_trigger_key =
        geo_.trigger_copy_map.find(bcm_pmt_key);
    if(bcm_trigger_key == geo_.trigger_copy_map.end()) {
        if(skip_missing_) {
            log_warn("Could not find PMT (%i/%u) in '%s'.trigger_copy_map! Skipping frame...",
                    bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                    geometry_key_.c_str());
            return false;
        } else
            log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                    bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                    geometry_key_.c_str());
    }

    double bcm_nim_pulse_time = 0.0;
    double bcm_max_nim_pulse_length = 0.0;

    NIMLogicPulseSeriesMap::const_iterator bcm_trigger_nim_pulses;
    if(correct_nim_pulse_time_) {
        bcm_trigger_nim_pulses = nim_pulses->find(bcm_trigger_key->second);
        if(bcm_trigger_nim_pulses == nim_pulses->end()) {
            if(skip_missing_) {
                log_warn("Could not find PMT (%i/%u) in '%s'! Skipping frame...",
                        bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                        nim_pulses_key_.c_str());
                return false;
            } else
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                        nim_pulses_key_.c_str());
        }
    } else {
        if(correct_nim_pulse_time_) {
            for(size_t j=0; j<bcm_trigger_nim_pulses->second.size(); ++j) {
                double length = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseLength();
                if(length > bcm_max_nim_pulse_length) {
                    bcm_nim_pulse_time = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseTime();
                    bcm_max_nim_pulse_length = length;
                }
            }
        }
    }

    // get the BCM time
    double bcm_time = bcm->bcm_start_time;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        Corrections pmt_corrections;
        double time_correction = 0.0;

        CCMOMGeo const & pmt = geo_.pmt_geo.at(pmt_key);
        CCMPMTType const & pmt_type = pmt.omtype;
        bool is_pmt = pmt_type == CCMPMTType::CCM8inCoated or pmt_type == CCMPMTType::CCM8inUncoated or pmt_type == CCMPMTType::CCM1in;

        if(correct_electron_transit_time_ and is_pmt) {
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
                calibration->pmtCal.find(pmt_key);
            if (calib == calibration->pmtCal.end()) {
                if(not allow_missing_) {
                    log_fatal("Could not find PMT (%i/%u) in '%s'",
                            pmt_key.GetRegion(), pmt_key.GetSensor(),
                            calibration_key_.c_str());
                }
            } else {
                CCMPMTCalibration const & pmt_calibration = calib->second;
                if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
                    time_correction -= pmt_calibration.GetPMTDeltaT();
                    pmt_corrections.electron_transit_time = true;
                }
            }
        }

        if(correct_nim_pulse_time_) {
            // retrieve the board trigger nim pulse time for this PMT
            std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key =
                geo_.trigger_copy_map.find(pmt_key);
            if(trigger_key == geo_.trigger_copy_map.end()) {
                if(not allow_missing_) {
                    log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                            pmt_key.GetRegion(), pmt_key.GetSensor(),
                            geometry_key_.c_str());
                }
            } else {
                NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
                    nim_pulses->find(trigger_key->second);
                if(pmt_trigger_nim_pulses == nim_pulses->end()) {
                    if(not allow_missing_) {
                        log_fatal("Could not find PMT (%i/%u) in '%s'",
                                pmt_key.GetRegion(), pmt_key.GetSensor(),
                                nim_pulses_key_.c_str());
                    }
                } else {
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
                        pmt_corrections.nim_pulse_time = true;
                    }
                    if(not std::isnan(bcm_nim_pulse_time)) {
                        time_correction -= (-bcm_nim_pulse_time);
                        pmt_corrections.bcm_nim_pulse_time = true;
                    }
                }
            }
        }

        if(not std::isnan(bcm_time)) {
            time_correction -= bcm_time;
            pmt_corrections.bcm_start_time = true;
        }

        output_indices[pmt_key] = (int32_t)(-time_correction / 2.0);
        output_corrections[pmt_key] = pmt_corrections;
    }

    return true;
}

bool AccumulateIndividualChannelWaveforms::ComputeReferenceIndicesUser(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections) {
    boost::shared_ptr<I3Map<CCMPMTKey, int32_t> const> reference_times = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, int32_t> const>>(reference_time_key_);
    if(reference_times == nullptr) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                reference_time_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                reference_time_key_.c_str());
    }

    CCMBCMSummaryConstPtr bcm = frame->Get<CCMBCMSummaryConstPtr>(bcm_summary_key_);
    if(!bcm) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                bcm_summary_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                bcm_summary_key_.c_str());
    }

    NIMLogicPulseSeriesMapConstPtr nim_pulses = frame->Get<NIMLogicPulseSeriesMapConstPtr>(nim_pulses_key_);
    if (!nim_pulses and correct_nim_pulse_time_) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                nim_pulses_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                nim_pulses_key_.c_str());
    }

    CCMCalibrationConstPtr calibration = frame->Get<CCMCalibrationConstPtr>(calibration_key_);
    if (!calibration and correct_electron_transit_time_) {
        if(skip_missing_) {
            log_warn("Couldn't find '%s' in the frame! Skipping frame...",
                calibration_key_.c_str());
            return false;
        } else
            log_fatal("Couldn't find '%s' in the frame!",
                calibration_key_.c_str());
    }

    // Find the BCM nim pulse
    CCMPMTKey bcm_pmt_key = CCMPMTKey(10, 1);
    std::map<CCMPMTKey, CCMTriggerKey>::const_iterator bcm_trigger_key =
        geo_.trigger_copy_map.find(bcm_pmt_key);
    if(bcm_trigger_key == geo_.trigger_copy_map.end()) {
        if(skip_missing_) {
            log_warn("Could not find PMT (%i/%u) in '%s'.trigger_copy_map! Skipping frame...",
                    bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                    geometry_key_.c_str());
            return false;
        } else
            log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                    bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                    geometry_key_.c_str());
    }

    double bcm_nim_pulse_time = 0.0;
    double bcm_max_nim_pulse_length = 0.0;

    NIMLogicPulseSeriesMap::const_iterator bcm_trigger_nim_pulses;
    if(correct_nim_pulse_time_) {
        bcm_trigger_nim_pulses = nim_pulses->find(bcm_trigger_key->second);
        if(bcm_trigger_nim_pulses == nim_pulses->end()) {
            if(skip_missing_) {
                log_warn("Could not find PMT (%i/%u) in '%s'! Skipping frame...",
                        bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                        nim_pulses_key_.c_str());
                return false;
            } else
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        bcm_pmt_key.GetRegion(), bcm_pmt_key.GetSensor(),
                        nim_pulses_key_.c_str());
        }
    } else {
        if(correct_nim_pulse_time_) {
            for(size_t j=0; j<bcm_trigger_nim_pulses->second.size(); ++j) {
                double length = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseLength();
                if(length > bcm_max_nim_pulse_length) {
                    bcm_nim_pulse_time = bcm_trigger_nim_pulses->second.at(j).GetNIMPulseTime();
                    bcm_max_nim_pulse_length = length;
                }
            }
        }
    }

    // get the BCM time
    double bcm_time = bcm->bcm_start_time;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        Corrections pmt_corrections;
        double time_correction = 0.0;

        CCMOMGeo const & pmt = geo_.pmt_geo.at(pmt_key);
        CCMPMTType const & pmt_type = pmt.omtype;
        bool is_pmt = pmt_type == CCMPMTType::CCM8inCoated or pmt_type == CCMPMTType::CCM8inUncoated or pmt_type == CCMPMTType::CCM1in;

        if(correct_electron_transit_time_ and is_pmt) {
            // retrieve the calibration for this PMT
            std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib =
                calibration->pmtCal.find(pmt_key);
            if (calib == calibration->pmtCal.end()) {
                if(not allow_missing_)
                    log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        calibration_key_.c_str());
            } else {
                CCMPMTCalibration const & pmt_calibration = calib->second;
                if(not std::isnan(pmt_calibration.GetPMTDeltaT())) {
                    time_correction -= pmt_calibration.GetPMTDeltaT();
                    pmt_corrections.electron_transit_time = true;
                }
            }
        }

        if(correct_nim_pulse_time_) {
            // retrieve the board trigger nim pulse time for this PMT
            std::map<CCMPMTKey, CCMTriggerKey>::const_iterator trigger_key =
                geo_.trigger_copy_map.find(pmt_key);
            if(trigger_key == geo_.trigger_copy_map.end()) {
                if(not allow_missing_)
                    log_fatal("Could not find PMT (%i/%u) in '%s'.trigger_copy_map",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        geometry_key_.c_str());
            } else {
                NIMLogicPulseSeriesMap::const_iterator pmt_trigger_nim_pulses =
                    nim_pulses->find(trigger_key->second);
                if(pmt_trigger_nim_pulses == nim_pulses->end()) {
                    if(not allow_missing_)
                        log_fatal("Could not find PMT (%i/%u) in '%s'",
                            pmt_key.GetRegion(), pmt_key.GetSensor(),
                            nim_pulses_key_.c_str());
                } else {
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
                        pmt_corrections.nim_pulse_time = true;
                    }
                    if(not std::isnan(bcm_nim_pulse_time)) {
                        time_correction -= (-bcm_nim_pulse_time);
                        pmt_corrections.bcm_nim_pulse_time = true;
                    }
                }
            }
        }

        if(not std::isnan(bcm_time)) {
            time_correction -= bcm_time;
            pmt_corrections.bcm_start_time = true;
        }

        if(reference_times and reference_times->count(pmt_key)) {
            time_correction -= (*reference_times).at(pmt_key);
            pmt_corrections.user_time = true;
        }

        output_indices[pmt_key] = (int32_t)(-time_correction / 2.0);
        output_corrections[pmt_key] = pmt_corrections;
    }

    return true;
}

bool AccumulateIndividualChannelWaveforms::ComputeReferenceIndices(I3FramePtr frame, std::map<uint32_t, int32_t> & output_indices, std::map<uint32_t, Corrections> & output_corrections) {
    switch(chosen_time_reference_) {
        case ReferenceTimeType::TriggerTime:
            return ComputeReferenceIndicesTrigger(frame, output_indices, output_corrections);
        case ReferenceTimeType::BCMTime:
            return ComputeReferenceIndicesBCM(frame, output_indices, output_corrections);
        case ReferenceTimeType::UserSpecifiedTime:
            return ComputeReferenceIndicesUser(frame, output_indices, output_corrections);
        default:
            break;
    };
    log_fatal("Received enum type that is not initialized properly!");
    return false;
}

void AccumulateIndividualChannelWaveforms::ProcessFrame(I3FramePtr frame) {
    std::map<uint32_t, int32_t> reference_indices;
    std::map<uint32_t, Corrections> corrections;
    bool computed_reference_indices = ComputeReferenceIndices(frame, reference_indices, corrections);
    if(not computed_reference_indices and skip_missing_)
        return PushFrame(frame);
    boost::shared_ptr<CCMWaveformUInt16Series const> waveform_raw = frame->Get<boost::shared_ptr<CCMWaveformUInt16Series const>>(waveforms_key_);
    boost::shared_ptr<CCMWaveformDoubleSeries const> waveform_cal = frame->Get<boost::shared_ptr<CCMWaveformDoubleSeries const>>(waveforms_key_);
    if(!waveform_raw and !waveform_cal)
        log_fatal("Couldn't find '%s' in the frame!",
                waveforms_key_.c_str());

    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const> baseline_estimates_map = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const>>(baseline_estimates_key_);
    boost::shared_ptr<I3Vector<uint32_t, BaselineEstimate> const> baseline_estimates_vec = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const>>(baseline_estimates_key_);

    I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map = geo_.pmt_channel_map;

    for(CCMPMTKey const & pmt_key : pmt_keys_) {
        if(pmt_channel_map.count(pmt_key) == 0) {
            log_fatal("Could not find PMT (%i/%u) in '%s'",
                    pmt_key.GetRegion(), pmt_key.GetSensor(),
                    geometry_key_.c_str());
        }
        uint32_t channel = pmt_channel_map.at(pmt_key);
        Corrections & pmt_corrections = corrections.at(pmt_key);
        if(waveform_raw) {
            if(channel >= waveform_raw->size()) {
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        waveforms_key_.c_str());
            }
            CCMWaveformUInt16 const & waveform = waveform_raw->at(channel);
            std::vector<double> wf(waveform.GetWaveform().begin(), waveform.GetWaveform().end());
            if(correct_baseline_and_invert_raw_waveforms_ and baseline_estimates and baseline_estimates->count(pmt_key)) {
                BaselineEstimate const & baseline_estimate = baseline_estimates->at(pmt_key);
                if(not std::isnan(baseline_estimate.baseline)) {
                    for(double & sample : wf) {
                        sample = -(sample + baseline_estimate.baseline);
                    }
                    pmt_corrections.baseline = true;
                }
            }
            if(allow_missing_) {
                if(corrections_.count(pmt_key)) {
                    Corrections const & previous_corrections = corrections_.at(pmt_key);
                    if(previous_corrections != pmt_corrections) {
                        log_fatal("Previous corrections for PMT (%i/%u) do not match current corrections! This should never happen!",
                                pmt_key.GetRegion(), pmt_key.GetSensor());
                    }
                } else {
                    corrections_[pmt_key] = pmt_corrections;
                }
            }
            accumulated_waveforms_[pmt_key].AddWaveform(wf, reference_indices[pmt_key]);
        } else if(waveform_cal) {
            if(channel >= waveform_cal->size()) {
                log_fatal("Could not find PMT (%i/%u) in '%s'",
                        pmt_key.GetRegion(), pmt_key.GetSensor(),
                        waveforms_key_.c_str());
            }
            CCMWaveformDouble const & waveform = waveform_cal->at(channel);
            if(allow_missing_) {
                if(corrections_.count(pmt_key)) {
                    Corrections const & previous_corrections = corrections_.at(pmt_key);
                    if(previous_corrections != pmt_corrections) {
                        log_fatal("Previous corrections for PMT (%i/%u) do not match current corrections! This should never happen!",
                                pmt_key.GetRegion(), pmt_key.GetSensor());
                    }
                } else {
                    corrections_[pmt_key] = pmt_corrections;
                }
            }
            accumulated_waveforms_[pmt_key].AddWaveform(waveform.GetWaveform(), reference_indices[pmt_key]);
        }
    }
}

AccumulateIndividualChannelWaveforms::AccumulateIndividualChannelWaveforms(const I3Context& context) : I3ConditionalModule(context) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string("CCMDAQConfig"));
    AddParameter("CCMCalibrationName", "Key for CCMCalibration", std::string(I3DefaultName<CCMCalibration>::value()));
    AddParameter("NIMPulsesName", "Key for NIMPulses", std::string("NIMPulses"));
    AddParameter("BCMSummaryName", "Key for BCMSummary", std::string("BCMSummary"));
    AddParameter("BaselineEstimatesName", "Key for BaselineEstimates", std::string("BaselineEstimates"));
    AddParameter("WaveformsKey", "Key for Waveforms", std::string("CCMWaveforms"));
    AddParameter("OutputPrefix", "Prefix for the module output", std::string("AccumulatedWaveforms"));
    AddParameter("ConsumeFrames", "Consume frames used as input?", bool(true));
    AddParameter("Channels", "Channels to run over", I3Vector<uint32_t>());
    AddParameter("TriggerReferenceTime", "Use the trigger time as a reference time? This is the default", bool(false));
    AddParameter("BCMReferenceTime", "Use the Beam Current Monitor start time as a reference time?", bool(false));
    AddParameter("ReferenceTimeKey", "Name of a frame key that contains an I3Map<uint32_t, int32_t> for the PMT reference indices", std::string(""));
    AddParameter("CorrectNIMPulseTime", "Correct for channel 15 NIM pulse arrival time?", bool(true));
    AddParameter("CorrectElectronTransitTime", "Correct for PMT electron transit time?", bool(true));
    AddParameter("CorrectBaselineAndInvertRawWaveforms", "Correct the baseline of raw waveforms and invert them?", bool(true));
    AddParameter("AllowMissingInformationPerPMT", "Allow information to be missing from the frame?", bool(false));
    AddParameter("SkipMissingInformation", "Skip frames that are missing information?", bool(false));
    AddParameter("OutputFrameType", "The type of frame to use in the ouptut. Default: DAQ", I3Frame::DAQ);
}

void AccumulateIndividualChannelWaveforms::Configure() {
    GetParameter("CCMGeometryName", geometry_key_);
    GetParameter("CCMDAQConfigName", daq_config_key_);
    GetParameter("CCMCalibrationName", calibration_key_);
    GetParameter("NIMPulsesName", nim_pulses_key_);
    GetParameter("BCMSummaryName", bcm_summary_key_);
    GetParameter("BaselineEstimatesName", baseline_estimates_key_);
    GetParameter("WaveformsKey", waveforms_key_);
    GetParameter("OutputPrefix", output_prefix_);
    GetParameter("ConsumeFrames", consume_frames_);
    GetParameter("Channels", allowed_channels_);
    GetParameter("TriggerReferenceTime", trigger_reference_time_);
    GetParameter("BCMReferenceTime", bcm_reference_time_);
    GetParameter("ReferenceTimeKey", reference_time_key_);
    GetParameter("CorrectNIMPulseTime", correct_nim_pulse_time_);
    GetParameter("CorrectElectronTransitTime", correct_electron_transit_time_);
    GetParameter("CorrectBaselineAndInvertRawWaveforms", correct_baseline_and_invert_raw_waveforms_);
    GetParameter("AllowMissingInformationPerPMT", allow_missing_);
    GetParameter("SkipMissingInformation", skip_missing_);
    GetParameter("OutputFrameType", output_frame_type_);

    // Preemtively sort the allowed_channels_ so they're ready for the Geometry function
    if(allowed_channels_.size() > 0)
        std::sort(allowed_channels_.begin(), allowed_channels_.end());

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

void AccumulateIndividualChannelWaveforms::Geometry(I3FramePtr frame) {
    // Dump the accumulated waveforms if we get a new geometry
    if(geo_seen_) {
        DumpWaveforms();
    }

    // Assumes allowed_channels_ is already sorted
    CCMGeometryConstPtr geo = frame->Get<CCMGeometryConstPtr>(geometry_key_);
    if (!geo)
        log_fatal("Couldn't find '%s' in the frame!",
                geometry_key_.c_str());
    geo_ = *geo;
    geo_seen_ = true;

    channel_pmt_map_.clear();
    std::map<uint32_t, CCMPMTKey> channel_pmt_map;
    for(std::pair<CCMPMTKey const, uint32_t> const & p : geo->pmt_channel_map) {
        channel_pmt_map[p.second] = p.first;
    }

    CCMAnalysis::Binary::CCMDAQConfigConstPtr daq_config = frame->Get<CCMAnalysis::Binary::CCMDAQConfigConstPtr>(daq_config_key_);
    if (!daq_config)
        log_fatal("Couldn't find '%s' in the frame!",
                daq_config_key_.c_str());

    // Copy and sort the channels from the daq config
    std::vector<uint32_t> geo_channels;
    uint32_t channel_index = 0;
    for(CCMAnalysis::Binary::DigitizerBoard const & board : daq_config->digitizer_boards) {
        std::string board_prefix = "physical_board_";
        size_t char_pos = board.physical_board_id.find(board_prefix);
        char_pos += board_prefix.size();
        int trigger_copy_number = std::atoi(board.physical_board_id.substr(char_pos, std::string::npos).c_str());
        CCMTriggerKey board_trigger_copy_key(CCMTriggerKey::TriggerType::BoardTriggerCopy, size_t(trigger_copy_number));
        for(CCMAnalysis::Binary::ChannelHeader const & channel : board.channels) {
            geo_channels.push_back(channel_index);
            channel_trigger_map_[channel_index] = board_trigger_copy_key;
            channel_index += 1;
        }
    }

    // Assume we're doing all PMTs if none are specified
    if(allowed_channels_.size() == 0) {
        channels_ = geo_channels;
        for(uint32_t const & channel : channels_) {
            if(channel_pmt_map.count(channel)) {
                channel_pmt_map_[channel] = channel_pmt_map[channel];
            }
        }
    } else {
        // Clear the final channel list
        channels_.clear();
        // Fill the final channel list with the intersection of specified channels and those available in the geometry
        // If allowed_channels_ is not sorted then this will fail horribly
        std::set_intersection(geo_channels.begin(), geo_channels.end(), allowed_channels_.begin(), allowed_channels_.end(), std::back_inserter(channels_));

        if(channels_.size() < allowed_channels_.size()) {
            std::vector<uint32_t> missing_channels;
            std::set_difference(allowed_channels_.begin(), allowed_channels_.end(), geo_channels.begin(), geo_channels.end(), std::back_inserter(missing_channels));
            std::stringstream ss;
            ss << "Some specified channels are not present in the geometry: ";
            for(uint32_t const & channel : allowed_channels_)
                ss << " " << channel;
            log_warn(ss.str().c_str());
        }
    }

    // Initialize the accumulated waveforms
    // Fill the channel_pmt_map_ with the PMT keys for each channel
    accumulated_waveforms_.clear();
    for(uint32_t const & channel : channels_) {
        if(accumulated_waveforms_.count(channel) == 0) {
            accumulated_waveforms_[channel] = WaveformAccumulator();
        }
        if(channel_pmt_map.count(channel)) {
            channel_pmt_map_[channel] = channel_pmt_map[channel];
        }
    }

    PushFrame(frame);
}

void AccumulateIndividualChannelWaveforms::DAQ(I3FramePtr frame) {
    ProcessFrame(frame);
    if(not consume_frames_)
        PushFrame(frame);
}

void AccumulateIndividualChannelWaveforms::DumpWaveforms() {
    boost::shared_ptr<I3Map<uint32_t, std::vector<double>>> p_samples = boost::make_shared<I3Map<uint32_t, std::vector<double>>>();
    boost::shared_ptr<I3Map<uint32_t, std::vector<uint32_t>>> p_counts = boost::make_shared<I3Map<uint32_t, std::vector<uint32_t>>>();
    boost::shared_ptr<I3Map<uint32_t, int32_t>> p_fixed_positions = boost::make_shared<I3Map<uint32_t, int32_t>>();

    for(uint32_t const & channel : channels_) {
        std::deque<double> deque_samples = accumulated_waveforms_[channel].GetSummedWaveform();
        std::deque<unsigned int> deque_counts = accumulated_waveforms_[channel].GetCounts();
        int32_t fixed_position = accumulated_waveforms_[channel].GetFixedPosition();

        p_samples->insert(std::make_pair<uint32_t, std::vector<double>>(uint32_t(channel), std::vector<double>(deque_samples.begin(), deque_samples.end())));
        p_counts->insert(std::make_pair<uint32_t, std::vector<uint32_t>>(uint32_t(channel), std::vector<uint32_t>(deque_counts.begin(), deque_counts.end())));
        p_fixed_positions->insert(std::make_pair<uint32_t, int32_t>(uint32_t(channel), std::move(fixed_position)));
    }

    I3FramePtr frame = boost::make_shared<I3Frame>(output_frame_type_);

    frame->Put((output_prefix_).c_str(), p_samples);
    frame->Put((output_prefix_ + "Counts").c_str(), p_counts);
    frame->Put((output_prefix_ + "FixedPosition").c_str(), p_fixed_positions);
    PushFrame(frame);
}

void AccumulateIndividualChannelWaveforms::Finish() {
    DumpWaveforms();
    Flush();
}
