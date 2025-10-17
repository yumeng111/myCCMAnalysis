#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <filesystem>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusBeamTime.h>
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime.h>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <phys-services/I3GSLRandomService.h>
#include <dataio/I3FrameSequence.h>
#include <CCMAnalysis/CCMBinary/BinaryFormat.h>

class OverlayNoise: public I3Module {
    std::string const default_noise_pulse_base_ = "WavedeformPulses";

    std::string input_reco_pulse_name_;
    std::string output_reco_pulse_name_;
    std::string noise_pulse_name_;

    bool make_noise_into_beam_time_pulses_ = false;
    bool make_noise_into_trigger_time_pulses_ = false;

    std::string input_simulated_board_time_offsets_name_;

    bool match_times_to_data_ = true;
    bool match_times_to_simulation_ = false;

    bool restrict_range_ = false;
    double min_noise_time_ = -10000;
    double max_noise_time_ = 6000;
    bool allow_stitching_ = false;
    double min_padding_time_ = 500;
    double max_padding_time_ = 1000;

    std::string randomServiceName_;
    I3RandomServicePtr randomService_;
    dataio::I3FrameSequencePtr frame_sequences;
    std::deque<I3FramePtr> frame_cache;
    std::vector<std::string> file_lists;

    boost::shared_ptr<CCMRecoPulseSeriesMap const> current_noise_pulses = nullptr;
    std::map<CCMPMTKey, double> current_frame_data_board_time_offsets;
    double current_start_time;
    double current_end_time;
public:
    OverlayNoise(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    bool NextFrame();
    std::tuple<boost::shared_ptr<CCMRecoPulseSeriesMap>, std::map<CCMPMTKey, double>> GetNoisePulses(double duration, double zero_time, bool zero_out_width = false);
};

I3_MODULE(OverlayNoise);

OverlayNoise::OverlayNoise(const I3Context& context) : I3Module(context),
    input_reco_pulse_name_("MCRecoPulses"), output_reco_pulse_name_("MCRecoPulsesPlusNoise"), noise_pulse_name_("TriggerTimePulses"),
    input_simulated_board_time_offsets_name_("SimulatedBoardTimeOffsets") {
    AddParameter("InputRecoPulseName", "Simulated reco pulse series name", input_reco_pulse_name_);
    AddParameter("OutputRecoPulseName", "Output reco pulse series name", output_reco_pulse_name_);
    AddParameter("NoisePulseName", "Data reco pulse series name", noise_pulse_name_);
    AddParameter("MakeNoiseIntoBeamTimePulses", "Convert noise pulses into beam time pulses", false);
    AddParameter("MakeNoiseIntoTriggerTimePulses", "Convert noise pulses into trigger time pulses", false);
    AddParameter("InputSimulatedBoardTimeOffsetsName", "Key for map of simulated board time offsets", input_simulated_board_time_offsets_name_);
    AddParameter("MatchTimesToData", "Match times to data", match_times_to_data_);
    AddParameter("MatchTimesToSimulation", "Match times to simulation", match_times_to_simulation_);

    AddParameter("RestrictRange", "Restrict the range of noise pulses to a specific time range", restrict_range_);
    AddParameter("MinNoiseTime", "Minimum time for noise pulses", min_noise_time_);
    AddParameter("MaxNoiseTime", "Maximum time for noise pulses", max_noise_time_);
    AddParameter("AllowStitching", "Allow noise pulses to be stitched together", allow_stitching_);
    AddParameter("MinPaddingTime", "Minimum padding time for stitching noise pulses", min_padding_time_);
    AddParameter("MaxPaddingTime", "Maximum padding time for stitching noise pulses", max_padding_time_);
    randomService_ = I3RandomServicePtr();
    AddParameter("RandomServiceName", "Name of the random service in the context. If empty default random service will be used.", randomServiceName_);
    AddParameter("FileLists", "Files to use for noise", file_lists);
}

void OverlayNoise::Configure() {
    GetParameter("InputRecoPulseName", input_reco_pulse_name_);
    GetParameter("OutputRecoPulseName", output_reco_pulse_name_);
    GetParameter("NoisePulseName", noise_pulse_name_);
    GetParameter("MakeNoiseIntoBeamTimePulses", make_noise_into_beam_time_pulses_);
    GetParameter("MakeNoiseIntoTriggerTimePulses", make_noise_into_trigger_time_pulses_);


    if(make_noise_into_beam_time_pulses_ and make_noise_into_trigger_time_pulses_)
        log_fatal("Cannot make noise pulses into both beam time and trigger time pulses");

    if(make_noise_into_beam_time_pulses_) {
        if(noise_pulse_name_ == "BeamTimePulses") {
            log_warn("Input noise pulses is set to \"BeamTimePulses\" and MakeNoiseIntoBeamTimePulses is true. The input will default to \"%s\" and be converted to beam time pulses", default_noise_pulse_base_.c_str());
        } else if( noise_pulse_name_ == "TriggerTimePulses") {
            log_fatal("Input noise pulses is set to \"TriggerTimePulses\" and MakeNoiseIntoBeamTimePulses is true. We can't double up on conversions and you probably want to pick the right one.");
        } else {
            log_info("Converting input noise pulses \"%s\" into beam time pulses.", noise_pulse_name_.c_str());
        }
    }

    if(make_noise_into_trigger_time_pulses_) {
        if(noise_pulse_name_ == "TriggerTimePulses") {
            log_warn("Input noise pulses is set to \"TriggerTimePulses\" and MakeNoiseIntoTriggerTimePulses is true. The input will default to \"%s\" and be converted to trigger time pulses", default_noise_pulse_base_.c_str());
        } else if( noise_pulse_name_ == "BeamTimePulses") {
            log_fatal("Input noise pulses is set to \"BeamTimePulses\" and MakeNoiseIntoTriggerTimePulses is true. We can't double up on conversions and you probably want to pick the right one.");
        } else {
            log_info("Converting input noise pulses \"%s\" into trigger time pulses.", noise_pulse_name_.c_str());
        }
    }

    GetParameter("InputSimulatedBoardTimeOffsetsName", input_simulated_board_time_offsets_name_);
    GetParameter("MatchTimesToData", match_times_to_data_);
    GetParameter("MatchTimesToSimulation", match_times_to_simulation_);

    if(match_times_to_data_ and match_times_to_simulation_) {
        log_fatal("Both MatchTimesToData and MatchTimesToSimulation are set to true. Please set only one of them.");
    }

    GetParameter("RestrictRange", restrict_range_);
    GetParameter("MinNoiseTime", min_noise_time_);
    GetParameter("MaxNoiseTime", max_noise_time_);
    GetParameter("AllowStitching", allow_stitching_);
    GetParameter("MinPaddingTime", min_padding_time_);
    GetParameter("MaxPaddingTime", max_padding_time_);

    GetParameter("RandomServiceName", randomServiceName_);
    if(randomServiceName_.empty()) {
        randomService_ = I3RandomServicePtr(new I3GSLRandomService(0));
        log_info("+ Random service: I3GSLRandomService  (default)");
    }
    else {
        randomService_ = GetContext().Get<I3RandomServicePtr>(randomServiceName_);
        if(randomService_) log_info("+ Random service: %s  (EXTERNAL)",  randomServiceName_.c_str());
        else log_fatal("No random service \"%s\" in context!", randomServiceName_.c_str());
    }
    GetParameter("FileLists", file_lists);

    // let's set up our frame sequences
    frame_sequences = boost::make_shared<dataio::I3FrameSequence>(file_lists);
}

bool IsBasicPulseSeries(I3FramePtr frame, std::string pulse_series_name) {
    if(not frame->Has(pulse_series_name))
        log_fatal("Frame does not have %s", pulse_series_name.c_str());

    if(not frame->Has("CCMDAQConfig"))
        log_fatal("Frame does not have CCMDAQConfig");

    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulse_series_name);
    CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr trigger_pulses = frame->Get<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr>(pulse_series_name);
    CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr beam_pulses = frame->Get<CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr>(pulse_series_name);
    CCMRecoPulseSeriesMapMaskConstPtr mask = frame->Get<CCMRecoPulseSeriesMapMaskConstPtr>(pulse_series_name);

    if(pulses == nullptr and trigger_pulses == nullptr and beam_pulses == nullptr)
        log_fatal("Frame key \"%s\" does not correspond to any of: CCMRecoPulseSeriesMap, CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime, CCMRecoPulseSeriesMapApplySPECalPlusBeamTime, CCMRecoPulseSeriesMapMask", pulse_series_name.c_str());

    return pulses != nullptr and trigger_pulses == nullptr and beam_pulses == nullptr and mask == nullptr;
}

std::tuple<bool, std::string> AddDummyPulseSeries(I3FramePtr frame, std::string pulse_series_name, size_t dummy_number) {
    CCMRecoPulseSeriesMapConstPtr pulses = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulse_series_name);
    CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr trigger_pulses = frame->Get<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr>(pulse_series_name);
    CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr beam_pulses = frame->Get<CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr>(pulse_series_name);
    CCMRecoPulseSeriesMapMaskConstPtr mask = frame->Get<CCMRecoPulseSeriesMapMaskConstPtr>(pulse_series_name);

    frame->Delete(pulse_series_name);

    std::stringstream ss;
    ss << "DummyPulses" << dummy_number + 1;
    std::string source_name = ss.str();
    ss.str("");
    ss << "DummyPulses" << dummy_number;
    std::string dest_name = ss.str();
    if(mask != nullptr) {
        return {false, mask->GetSource()};
    } else if(trigger_pulses != nullptr) {
        CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>(source_name, trigger_pulses->GetCalibrationSource(), trigger_pulses->GetNIMPulsesSource(), trigger_pulses->GetGeometrySource());
        frame->Put(dest_name, pulses_view);
        return {true, trigger_pulses->GetPulsesSource()};
    } else if(beam_pulses != nullptr) {
        CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime>(source_name, beam_pulses->GetCalibrationSource(), beam_pulses->GetNIMPulsesSource(), beam_pulses->GetGeometrySource(), beam_pulses->GetBCMSummarySource());
        frame->Put(dest_name, pulses_view);
        return {true, beam_pulses->GetPulsesSource()};
    } else {
        return {false, ""};
    }
}

std::map<CCMPMTKey, double> GetDataBoardOffsets(CCMRecoPulseSeriesMapConstPtr pulses) {
    std::map<CCMPMTKey, double> data_board_time_offsets;
    for(auto const & [key, series] : *pulses) {
        if(series.empty()) {
            continue;
        }
        double offset = -(series.front().GetTime());
        data_board_time_offsets.insert(std::make_pair(key, offset));
    }
    return data_board_time_offsets;
}

std::tuple<double, double> GetPulsesTimeRange(CCMRecoPulseSeriesMapConstPtr pulses, bool expansive = false) {
    double min_time;
    double max_time;
    bool found_valid_pulse = false;

    if(expansive) {
        min_time = std::numeric_limits<double>::max();
        max_time = -std::numeric_limits<double>::max();

        for(auto const & [key, series] : *pulses) {
            if(series.empty()) {
                continue;
            }
            found_valid_pulse = true;
            double const & start_time = series.front().GetTime();
            double const & end_time = series.back().GetTime();
            if(start_time < min_time) {
                min_time = start_time;
            }
            if(end_time > max_time) {
                max_time = end_time;
            }
        }
    } else {
        min_time = -std::numeric_limits<double>::max();
        max_time = std::numeric_limits<double>::max();
        for(auto const & [key, series] : *pulses) {
            if(series.empty()) {
                continue;
            }
            found_valid_pulse = true;
            double const & start_time = series.front().GetTime();
            double const & end_time = series.back().GetTime();
            if(start_time > min_time) {
                min_time = start_time;
            }
            if(end_time < max_time) {
                max_time = end_time;
            }
        }
    }
    if(not found_valid_pulse) {
        min_time = 0;
        max_time = 0;
    }
    return {min_time, max_time};
}

std::tuple<double, double, std::map<CCMPMTKey, double>> GetTriggerExtent(I3FramePtr frame, std::string pulse_series_name, bool expansive = false) {
    if(not frame->Has(pulse_series_name)) {
        log_fatal("Frame does not have %s", pulse_series_name.c_str());
    }
    if(not frame->Has("CCMDAQConfig")) {
        log_fatal("Frame does not have CCMDAQConfig");
    }

    I3FramePtr dummy_frame = boost::make_shared<I3Frame>(*frame);
    size_t dummy_number = 0;
    std::string next_pulse_series_name = pulse_series_name;

    // For each layer of pulse series correction we need to put a corresponding correction in place for our dummy pulse series
    // The final result will be read out from DummyPulses0
    while(not IsBasicPulseSeries(frame, next_pulse_series_name)) {
        std::tuple<bool, std::string> x = AddDummyPulseSeries(dummy_frame, next_pulse_series_name, dummy_number);
        bool deeper = std::get<0>(x);
        next_pulse_series_name = std::get<1>(x);
        if(deeper) {
            ++dummy_number;
        }
    }
    std::string const & source_pulse_series_name = next_pulse_series_name;

    CCMAnalysis::Binary::CCMDAQConfigConstPtr daq_config = frame->Get<CCMAnalysis::Binary::CCMDAQConfigConstPtr>("CCMDAQConfig");

    // Prepare the dummy input pulse series
    CCMRecoPulseSeriesMapPtr dummy_pulses = boost::make_shared<CCMRecoPulseSeriesMap>();
    CCMRecoPulse first_pulse;
    CCMRecoPulse last_pulse;
    first_pulse.SetTime(0);
    last_pulse.SetTime(daq_config->machine_configurations[0].num_samples * 2);
    CCMRecoPulseSeries dummy_series{first_pulse, last_pulse};

    CCMRecoPulseSeriesMapConstPtr pulses = dummy_frame->Get<CCMRecoPulseSeriesMapConstPtr>(source_pulse_series_name);
    for(auto & [key, _] : *pulses) {
        dummy_pulses->insert(std::make_pair(key, dummy_series));
    }

    std::stringstream ss;
    ss << "DummyPulses" << dummy_number;
    std::string dummy_pulse_series_name = ss.str();

    // Put the dummy pulse series in the frame
    dummy_frame->Put(dummy_pulse_series_name, dummy_pulses);

    CCMRecoPulseSeriesMapConstPtr corrected_pulses = dummy_frame->Get<CCMRecoPulseSeriesMapConstPtr>("DummyPulses0");

    std::tuple<double, double> range = GetPulsesTimeRange(corrected_pulses, expansive);
    std::map<CCMPMTKey, double> data_board_time_offsets = GetDataBoardOffsets(corrected_pulses);

    double min_time = std::get<0>(range);
    double max_time = std::get<1>(range);
    double mod = min_time - std::floor(min_time / 2) * 2;
    min_time = min_time - mod;
    for(auto & [key, offset] : data_board_time_offsets) {
        offset = offset - min_time;
    }

    return {std::get<0>(range), std::get<1>(range), data_board_time_offsets};
}

std::map<CCMPMTKey, double> GetBoardOffsetCorrections(std::map<CCMPMTKey, double> const & from, std::map<CCMPMTKey, double> const & to, double mod) {
    std::map<CCMPMTKey, double> board_offset_corrections;
    for(auto const & [key, from_offset] : from) {
        board_offset_corrections.insert({key, 0.0});
    }
    if(to.size() == 0)
        return board_offset_corrections;

    size_t count = 0;
    double to_avg = 0.0;
    double from_avg = 0.0;
    for(auto const & [key, from_offset] : from) {
        std::map<CCMPMTKey, double>::const_iterator to_it = to.find(key);
        if(to_it == to.end()) {
            continue;
        }
        double const & to_offset = to_it->second;
        to_avg += to_offset;
        from_avg += from_offset;
        ++count;
    }
    to_avg /= count;
    from_avg /= count;

    double min_diff = std::numeric_limits<double>::max();
    CCMPMTKey min_diff_key;
    for(auto const & [key, from_offset] : from) {
        std::map<CCMPMTKey, double>::const_iterator to_it = to.find(key);
        if(to_it == to.end()) {
            continue;
        }
        double const & to_offset = to_it->second;
        double to_diff = std::abs(to_offset - to_avg);
        double from_diff = std::abs(from_offset - from_avg);
        double min = std::min(to_diff, from_diff);
        if(min < min_diff) {
            min_diff = min;
            min_diff_key = key;
        }
    }

    double to_reference_time = to.at(min_diff_key);
    double to_mod = std::fmod(to_reference_time, mod);
    if(to_mod < 0)
        to_mod += mod;
    to_reference_time -= to_mod;

    for(auto const & [key, from_offset] : from) {
        std::map<CCMPMTKey, double>::const_iterator to_it = to.find(key);
        if(to_it == to.end()) {
            continue;
        }
        double const & to_offset = to_it->second;
        double offset = from_offset - (to_offset - to_reference_time);
        board_offset_corrections.at(key) = offset;
    }
    return board_offset_corrections;
}

bool OverlayNoise::NextFrame() {
    bool failed = false;
    I3FramePtr frame;
    // Grab frames until we get a DAQ frame
    while(true) {
        // Fail if we do not have any more frames
        if(not frame_sequences->more()) {
            if(failed) {
                log_fatal("No DAQ frames found in the input files");
            }
            frame_sequences->rewind();
            failed = true;
            continue;
        }
        frame = frame_sequences->pop_frame();
        if(frame->GetStop() != I3Frame::DAQ) {
            // Skip non-DAQ frames
            continue;
        }

        if(make_noise_into_beam_time_pulses_) {
            // Check if the beam time pulses are already in the frame
            CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr beam_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime const>>(noise_pulse_name_);
            if(beam_pulses != nullptr) {
                // First check if the BCMSummary is in the frame
                std::string const bcm_summary_key = beam_pulses->GetBCMSummarySource();
                if(not frame->Has(bcm_summary_key)) {
                    log_warn("Frame does not have \"%s\", skipping frame.", bcm_summary_key.c_str());
                    continue;
                }
                // If so, go ahead and use them as is
                current_noise_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);
                log_warn("Input noise pulses \"%s\" are already in the frame as beam time pulses, skipping conversion. Consider setting MakeNoiseIntoBeamTimePulses to false.", noise_pulse_name_.c_str());
            } else {
                // First check if the BCMSummary is in the frame
                if(not frame->Has("BCMSummary")) {
                    log_warn("Frame does not have \"BCMSummary\", skipping frame.");
                    continue;
                }
                // We need to convert the existing pulse series into beam time pulses
                current_noise_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);
                std::string name;
                if(current_noise_pulses == nullptr and noise_pulse_name_ == "BeamTimePulses") {
                    log_warn("Defaulting to \"%s\" for noise pulses and converting to beam time pulses", default_noise_pulse_base_.c_str());
                    name = default_noise_pulse_base_;
                } else {
                    if(noise_pulse_name_ == "BeamTimePulses") {
                        log_fatal("Found a \"BeamTimePulses\" key that is not a CCMRecoPulseSeriesMapApplySPECalPlusBeamTime. This is unexpected and likely bad if you're trying to re-apply the corrections, please check your input.");
                    }
                    name = noise_pulse_name_;
                }
                if(not frame->Has(name)) {
                    log_warn("Frame does not have \"%s\", skipping frame.", name.c_str());
                    continue;
                }
                CCMRecoPulseSeriesMapApplySPECalPlusBeamTimePtr beam_pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime>(name, "CCMCalibration", "NIMPulses", "CCMGeometry", "BCMSummary");
                current_noise_pulses = beam_pulses_view->Apply(*frame);
            }
        } else if(make_noise_into_trigger_time_pulses_) {
            // Check if the trigger time pulses are already in the frame
            CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr trigger_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime const>>(noise_pulse_name_);
            if(trigger_pulses != nullptr) {
                // If so, go ahead and use them as is
                current_noise_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);
                log_warn("Input noise pulses \"%s\" are already in the frame as trigger time pulses, skipping conversion. Consider setting MakeNoiseIntoTriggerTimePulses to false.", noise_pulse_name_.c_str());
            } else {
                // We need to convert the existing pulse series into trigger time pulses
                current_noise_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);
                std::string name;
                if(current_noise_pulses == nullptr and noise_pulse_name_ == "TriggerTimePulses") {
                    log_warn("Defaulting to \"%s\" for noise pulses and converting to trigger time pulses", default_noise_pulse_base_.c_str());
                    name = default_noise_pulse_base_;
                } else {
                    if(noise_pulse_name_ == "TriggerTimePulses") {
                        log_fatal("Found a \"TriggerTimePulses\" key that is not a CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime. This is unexpected and likely bad if you're trying to re-apply the corrections, please check your input.");
                    }
                    name = noise_pulse_name_;
                }
                if(not frame->Has(name)) {
                    log_warn("Frame does not have \"%s\", skipping frame.", name.c_str());
                    continue;
                }
                CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimePtr trigger_pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>(name, "CCMCalibration", "NIMPulses", "CCMGeometry");
                current_noise_pulses = trigger_pulses_view->Apply(*frame);
            }
        } else {
            if(not frame->Has(noise_pulse_name_)) {
                log_warn("Frame does not have \"%s\", skipping frame.", noise_pulse_name_.c_str());
                continue;
            }
            current_noise_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);
        }

        break;
    }

    std::tuple<double, double, std::map<CCMPMTKey, double>> noise_extent = GetTriggerExtent(frame, noise_pulse_name_, false);
    current_start_time = std::get<0>(noise_extent);
    current_end_time = std::get<1>(noise_extent);
    current_frame_data_board_time_offsets = std::get<2>(noise_extent);

    if(restrict_range_) {
        current_start_time = std::max(current_start_time, min_noise_time_);
        current_end_time = std::min(current_end_time, max_noise_time_);
    }
    return failed;
}

std::tuple<boost::shared_ptr<CCMRecoPulseSeriesMap>, std::map<CCMPMTKey, double>> OverlayNoise::GetNoisePulses(double duration, double zero_time, bool zero_out_width) {
    zero_out_width = not zero_out_width;

    if(not allow_stitching_) {
        double remaining_time = current_end_time - current_start_time;
        bool already_looped = false;
        while(remaining_time < duration) {
            bool looped = NextFrame();
            if(already_looped and looped) {
                log_fatal("Not enough noise pulses to cover the duration");
            }
            already_looped |= looped;
            remaining_time = current_end_time - current_start_time;
        }
        boost::shared_ptr<CCMRecoPulseSeriesMap> noise_pulses = boost::make_shared<CCMRecoPulseSeriesMap>();
        double random_start_time = randomService_->Uniform(current_start_time, current_end_time - duration);
        double random_end_time = random_start_time + duration;
        for(auto const & [key, current_noise_series] : *current_noise_pulses) {
            if(current_noise_series.empty()) {
                continue;
            }
            CCMRecoPulseSeries noise_series;
            for(CCMRecoPulse const & pulse : current_noise_series) {
                double const & time = pulse.GetTime();
                if(time >= random_start_time and time < random_end_time) {
                    noise_series.push_back(pulse);
                    CCMRecoPulse & noise = noise_series.back();
                    noise.SetTime(time - random_start_time + zero_time);
                    noise.SetWidth(noise.GetWidth() * zero_out_width);
                }
            }
            noise_pulses->insert(std::make_pair(key, noise_series));
        }
        current_start_time = random_start_time;
        return {noise_pulses, current_frame_data_board_time_offsets};
    }


    double accumulated_time = 0.0;
    boost::shared_ptr<CCMRecoPulseSeriesMap> noise_pulses = boost::make_shared<CCMRecoPulseSeriesMap>();

    std::map<CCMPMTKey, double> target_board_time_offsets;
    if(current_noise_pulses != nullptr) {
        target_board_time_offsets = current_frame_data_board_time_offsets;
    }

    std::map<CCMPMTKey, double> board_offset_corrections;

    while(accumulated_time < duration) {
        if(current_noise_pulses == nullptr) {
            NextFrame();
            target_board_time_offsets = current_frame_data_board_time_offsets;
        }
        double remaining_time = current_end_time - current_start_time;
        double random_start_time = randomService_->Uniform(
                current_start_time,
                std::max(current_start_time, current_end_time - (duration - accumulated_time))
        );
        double random_end_time = std::min(current_end_time, random_start_time + std::min(remaining_time, duration - accumulated_time));
        for(auto const & [key, current_noise_series] : *current_noise_pulses) {
            if(not noise_pulses->count(key)) {
                noise_pulses->insert(std::make_pair(key, CCMRecoPulseSeries()));
            }
            double offset_correction = 0;
            std::map<CCMPMTKey, double>::const_iterator offset_it = board_offset_corrections.find(key);
            if(offset_it != board_offset_corrections.end()) {
                offset_correction = offset_it->second;
            }

            CCMRecoPulseSeries & noise_series = noise_pulses->at(key);
            for(CCMRecoPulse const & pulse : current_noise_series) {
                double const & time = pulse.GetTime();
                if(time >= random_start_time and time < random_end_time) {
                    noise_series.push_back(pulse);
                    noise_series.back().SetTime(time - random_start_time + accumulated_time + zero_time + offset_correction);
                }
            }
            noise_pulses->insert(std::make_pair(key, noise_series));
        }
        accumulated_time += random_end_time - random_start_time;
        current_start_time = random_end_time;
        if(random_end_time >= current_end_time) {
            NextFrame();
            board_offset_corrections.clear();
            board_offset_corrections = GetBoardOffsetCorrections(current_frame_data_board_time_offsets, target_board_time_offsets, 2.0);
        }
    }

    return {noise_pulses, target_board_time_offsets};
}

void OverlayNoise::DAQ(I3FramePtr frame) {

    // read in our reco pulse series
    boost::shared_ptr<CCMRecoPulseSeriesMap const> input_reco_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_reco_pulse_name_);
    I3MapPMTKeyDoubleConstPtr sim_board_time_offsets_ptr = frame->Get<I3MapPMTKeyDoubleConstPtr>(input_simulated_board_time_offsets_name_);

    std::tuple<double, double> event_time_range = GetPulsesTimeRange(input_reco_pulses, true);

    double event_start_time = std::get<0>(event_time_range);
    double event_end_time = std::get<1>(event_time_range);
    double event_duration = event_end_time - event_start_time;

    double noise_pre_padding = randomService_->Uniform(min_padding_time_, max_padding_time_);
    double noise_post_padding = randomService_->Uniform(min_padding_time_, max_padding_time_);
    double noise_duration = event_duration + noise_pre_padding + noise_post_padding;

    // Get noise pulses with the correct time offset
    CCMRecoPulseSeriesMapPtr combined_pulses;
    std::map<CCMPMTKey, double> data_board_time_offsets;
    std::tie(combined_pulses, data_board_time_offsets) = GetNoisePulses(noise_duration, event_start_time, true);

    std::map<CCMPMTKey, double> board_offset_corrections;
    if(match_times_to_simulation_) {
        board_offset_corrections = GetBoardOffsetCorrections(data_board_time_offsets, *sim_board_time_offsets_ptr, 2.0);
    } else if(match_times_to_data_) {
        board_offset_corrections = GetBoardOffsetCorrections(*sim_board_time_offsets_ptr, data_board_time_offsets, 2.0);
    }

    if(match_times_to_simulation_) {
        for(auto & [key, pulses] : *combined_pulses) {
            std::map<CCMPMTKey, double>::const_iterator offset_it = board_offset_corrections.find(key);
            if(offset_it == board_offset_corrections.end()) {
                continue;
            }
            double const & offset = offset_it->second;
            for(CCMRecoPulse & pulse : pulses) {
                pulse.SetTime(pulse.GetTime() + offset);
            }
        }
    }

    // Copy the event pulses into the destination map
    for(auto const & [key, source_series] : *input_reco_pulses) {
        // Find the corresponding PMT in the destination map
        CCMRecoPulseSeriesMap::iterator it_dest = combined_pulses->find(key);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == combined_pulses->end()) {
            combined_pulses->insert(std::make_pair(key, CCMRecoPulseSeries()));
            // Update the iterator so it points to our new entry
            it_dest = combined_pulses->find(key);
        }

        // Reference to the destination
        CCMRecoPulseSeries & dest_series = it_dest->second;

        if(match_times_to_data_) {
            std::map<CCMPMTKey, double>::const_iterator offset_it = board_offset_corrections.find(key);
            if(offset_it != board_offset_corrections.end()) {
                double const & offset = offset_it->second;
                for(CCMRecoPulse pulse : source_series) { // Copy the pulse so we can modify it
                    double time = pulse.GetTime();
                    pulse.SetTime(time + offset);
                    dest_series.push_back(pulse);
                }
            } else {
                for(CCMRecoPulse const & pulse : source_series) {
                    dest_series.push_back(pulse);
                }
            }
        } else {
            // Iterate over the vector of CCMRecoPulse in the source map for this PMT
            for(CCMRecoPulse const & pulse : source_series) {
                dest_series.push_back(pulse);
            }
        }

        // Done combining pulse series! Let's sort according to time
        std::sort(dest_series.begin(), dest_series.end(), [](const CCMRecoPulse& a, const CCMRecoPulse& b) { return a.GetTime() < b.GetTime(); });
    }

    // Save
    //TODO save the noise series
    //TODO Add a charge cut on stitched frames
    frame->Put(output_reco_pulse_name_, combined_pulses);
    PushFrame(frame);

}

