#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <tuple>
#include <cctype>
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusBeamTime.h>
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime.h>
#include <dataclasses/CCMRecoPulseSeriesMapApplyOffsets.h>
#include <dataclasses/I3MapCCMPMTKeyUnion.h>
#include <dataclasses/I3MapCCMPMTKeyMask.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <phys-services/I3GSLRandomService.h>
#include <dataio/I3FrameSequence.h>
#include <CCMAnalysis/CCMBinary/BinaryFormat.h>

struct PulsesSegment {
    CCMRecoPulseSeries original_pulses;
    double original_offset;
    double original_start_time;
    double target_offset;
    double target_start_time;
};

struct PulsesMapSegment {
    boost::shared_ptr<CCMRecoPulseSeriesMap const> pulses = nullptr;
    std::map<CCMPMTKey, double> offsets;
    double start_time = 0;
    double end_time = 0;
    double target_start_time = 0;

    PulsesMapSegment() = default;
    PulsesMapSegment(PulsesMapSegment const & other) = default;
    PulsesMapSegment(PulsesMapSegment const & other, double target_start_time) :
        PulsesMapSegment(other) {
        this->target_start_time = target_start_time;
    }
    double duration() const {
        return end_time - start_time;
    }
};

class OverlayNoise: public I3Module {
    std::string const default_noise_pulse_base_ = "WavedeformPulses";

    std::string input_sim_pulses_name_;
    std::string input_sim_offsets_name_;

    std::string input_noise_pulses_name_;

    std::string output_pulses_name_;

    bool make_noise_into_beam_time_pulses_ = false;
    bool make_noise_into_trigger_time_pulses_ = false;


    bool match_times_to_data_ = true;
    bool match_times_to_simulation_ = false;

    bool restrict_range_ = false;
    double min_noise_time_ = -10000;
    double max_noise_time_ = 6000;
    bool allow_stitching_ = false;
    double min_padding_time_ = 500;
    double max_padding_time_ = 1000;

    static constexpr double reco_bin_width_ = 2.0; // ns

    std::string randomServiceName_;
    I3RandomServicePtr randomService_;
    dataio::I3FrameSequencePtr frame_sequences;
    std::deque<I3FramePtr> frame_cache;
    std::vector<std::string> file_lists;

    PulsesMapSegment current_noise_segment;

    //boost::shared_ptr<CCMRecoPulseSeriesMap const> current_noise_pulses = nullptr;
    //std::map<CCMPMTKey, double> current_frame_data_board_time_offsets;
    //double current_start_time;
    //double current_end_time;
public:
    OverlayNoise(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    bool NextFrame();
    std::tuple<std::vector<PulsesMapSegment>, std::map<CCMPMTKey, double>> GetNoisePulses(double duration, double zero_time, bool zero_out_width = false);
    static void PulsesFromSegments(std::vector<PulsesMapSegment> const & segments, std::map<CCMPMTKey, double> const & target_offsets, std::map<CCMPMTKey, CCMRecoPulseSeries> & output, bool zero_out_width = false);
    static inline double Mod(double x, double mod) {
        double result = std::fmod(x, mod);
        result += (result < 0) ? mod : 0;
        return result;
    }
    static inline double SnapToMod(double sum, double original_offset, double target, double mod) {
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

    static inline void SnapPulsesToMod(std::vector<CCMRecoPulse> & pulses, double shift, double original_offset, double target, double mod, bool zero_out_width=false) {
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
};

I3_MODULE(OverlayNoise);

OverlayNoise::OverlayNoise(const I3Context& context) : I3Module(context),
    input_sim_pulses_name_("MCRecoPulses"), output_pulses_name_("MCRecoPulsesPlusNoise"), input_noise_pulses_name_("TriggerTimePulses"),
    input_sim_offsets_name_("SimulatedBoardTimeOffsets") {
    AddParameter("InputRecoPulseName", "Simulated reco pulse series name", input_sim_pulses_name_);
    AddParameter("OutputRecoPulseName", "Output reco pulse series name", output_pulses_name_);
    AddParameter("NoisePulseName", "Data reco pulse series name", input_noise_pulses_name_);
    AddParameter("MakeNoiseIntoBeamTimePulses", "Convert noise pulses into beam time pulses", false);
    AddParameter("MakeNoiseIntoTriggerTimePulses", "Convert noise pulses into trigger time pulses", false);
    AddParameter("InputSimulatedBoardTimeOffsetsName", "Key for map of simulated board time offsets", input_sim_offsets_name_);
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
    GetParameter("InputRecoPulseName", input_sim_pulses_name_);
    GetParameter("OutputRecoPulseName", output_pulses_name_);
    GetParameter("NoisePulseName", input_noise_pulses_name_);
    GetParameter("MakeNoiseIntoBeamTimePulses", make_noise_into_beam_time_pulses_);
    GetParameter("MakeNoiseIntoTriggerTimePulses", make_noise_into_trigger_time_pulses_);


    if(make_noise_into_beam_time_pulses_ and make_noise_into_trigger_time_pulses_)
        log_fatal("Cannot make noise pulses into both beam time and trigger time pulses");

    if(make_noise_into_beam_time_pulses_) {
        if(input_noise_pulses_name_ == "BeamTimePulses") {
            log_warn("Input noise pulses is set to \"BeamTimePulses\" and MakeNoiseIntoBeamTimePulses is true. The input will default to \"%s\" and be converted to beam time pulses", default_noise_pulse_base_.c_str());
        } else if( input_noise_pulses_name_ == "TriggerTimePulses") {
            log_fatal("Input noise pulses is set to \"TriggerTimePulses\" and MakeNoiseIntoBeamTimePulses is true. We can't double up on conversions and you probably want to pick the right one.");
        } else {
            log_info("Converting input noise pulses \"%s\" into beam time pulses.", input_noise_pulses_name_.c_str());
        }
    }

    if(make_noise_into_trigger_time_pulses_) {
        if(input_noise_pulses_name_ == "TriggerTimePulses") {
            log_warn("Input noise pulses is set to \"TriggerTimePulses\" and MakeNoiseIntoTriggerTimePulses is true. The input will default to \"%s\" and be converted to trigger time pulses", default_noise_pulse_base_.c_str());
        } else if( input_noise_pulses_name_ == "BeamTimePulses") {
            log_fatal("Input noise pulses is set to \"BeamTimePulses\" and MakeNoiseIntoTriggerTimePulses is true. We can't double up on conversions and you probably want to pick the right one.");
        } else {
            log_info("Converting input noise pulses \"%s\" into trigger time pulses.", input_noise_pulses_name_.c_str());
        }
    }

    GetParameter("InputSimulatedBoardTimeOffsetsName", input_sim_offsets_name_);
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

std::tuple<double, double, std::map<CCMPMTKey, double>> GetTriggerExtent(I3FramePtr frame, std::string pulse_series_name, double bin_width, bool expansive = false) {
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
    last_pulse.SetTime(daq_config->machine_configurations[0].num_samples * bin_width);
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

    return {std::get<0>(range), std::get<1>(range), data_board_time_offsets};
}

std::map<CCMPMTKey, double> GetBoardOffsetCorrections(std::map<CCMPMTKey, double> const & from, std::map<CCMPMTKey, double> const & to, double mod) {
    std::map<CCMPMTKey, double> board_offset_corrections;
    for(auto const & [key, from_offset] : from) {
        board_offset_corrections.insert({key, 0.0});
    }
    if(to.size() == 0)
        return board_offset_corrections;
    if(from.size() == 0)
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
    size_t n_skipped = 0;
    bool do_warn = false;
    I3FramePtr frame;
    I3FramePtr extent_frame;
    std::string extent_pulse_series_name;
    // Grab frames until we get a DAQ frame
    while(true) {
        if(n_skipped == 10) {
            log_warn("Skipped 10 frames while searching for DAQ frames, enabling warnings...");
            do_warn = true;
        }
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
            ++n_skipped;
            continue;
        }
        extent_frame = boost::make_shared<I3Frame>(*frame);

        if(make_noise_into_beam_time_pulses_) {
            // Check if the beam time pulses are already in the frame
            CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr beam_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime const>>(input_noise_pulses_name_);
            if(beam_pulses != nullptr) {
                // First check if the BCMSummary is in the frame
                std::string const bcm_summary_key = beam_pulses->GetBCMSummarySource();
                if(not frame->Has(bcm_summary_key)) {
                    if(do_warn)
                        log_warn("Frame does not have \"%s\", skipping frame.", bcm_summary_key.c_str());
                    ++n_skipped;
                    continue;
                }
                // If so, go ahead and use them as is
                current_noise_segment.pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_noise_pulses_name_);
                extent_pulse_series_name = input_noise_pulses_name_;
                log_warn("Input noise pulses \"%s\" are already in the frame as beam time pulses, skipping conversion. Consider setting MakeNoiseIntoBeamTimePulses to false.", input_noise_pulses_name_.c_str());
            } else {
                // First check if the BCMSummary is in the frame
                if(not frame->Has("BCMSummary")) {
                    if(do_warn)
                        log_warn("Frame does not have \"BCMSummary\", skipping frame.");
                    ++n_skipped;
                    continue;
                }
                // We need to convert the existing pulse series into beam time pulses
                current_noise_segment.pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_noise_pulses_name_);
                std::string name;
                if(current_noise_segment.pulses == nullptr and input_noise_pulses_name_ == "BeamTimePulses") {
                    log_warn("Defaulting to \"%s\" for noise pulses and converting to beam time pulses", default_noise_pulse_base_.c_str());
                    name = default_noise_pulse_base_;
                } else {
                    if(input_noise_pulses_name_ == "BeamTimePulses") {
                        log_fatal("Found a \"BeamTimePulses\" key that is not a CCMRecoPulseSeriesMapApplySPECalPlusBeamTime. This is unexpected and likely bad if you're trying to re-apply the corrections, please check your input.");
                    }
                    name = input_noise_pulses_name_;
                }
                if(not frame->Has(name)) {
                    if(do_warn)
                        log_warn("Frame does not have \"%s\", skipping frame.", name.c_str());
                    ++n_skipped;
                    continue;
                }
                CCMRecoPulseSeriesMapApplySPECalPlusBeamTimePtr beam_pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime>(name, "CCMCalibration", "NIMPulses", "CCMGeometry", "BCMSummary");
                current_noise_segment.pulses = beam_pulses_view->Apply(*frame);
                extent_pulse_series_name = "DummyBeamTimePulsesView";
                extent_frame->Put(extent_pulse_series_name, beam_pulses_view);
            }
        } else if(make_noise_into_trigger_time_pulses_) {
            // Check if the trigger time pulses are already in the frame
            CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr trigger_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime const>>(input_noise_pulses_name_);
            if(trigger_pulses != nullptr) {
                // If so, go ahead and use them as is
                current_noise_segment.pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_noise_pulses_name_);
                log_warn("Input noise pulses \"%s\" are already in the frame as trigger time pulses, skipping conversion. Consider setting MakeNoiseIntoTriggerTimePulses to false.", input_noise_pulses_name_.c_str());
            } else {
                // We need to convert the existing pulse series into trigger time pulses
                current_noise_segment.pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_noise_pulses_name_);
                std::string name;
                if(current_noise_segment.pulses == nullptr and input_noise_pulses_name_ == "TriggerTimePulses") {
                    log_warn("Defaulting to \"%s\" for noise pulses and converting to trigger time pulses", default_noise_pulse_base_.c_str());
                    name = default_noise_pulse_base_;
                } else {
                    if(input_noise_pulses_name_ == "TriggerTimePulses") {
                        log_fatal("Found a \"TriggerTimePulses\" key that is not a CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime. This is unexpected and likely bad if you're trying to re-apply the corrections, please check your input.");
                    }
                    name = input_noise_pulses_name_;
                }
                if(not frame->Has(name)) {
                    if(do_warn)
                        log_warn("Frame does not have \"%s\", skipping frame.", name.c_str());
                    ++n_skipped;
                    continue;
                }
                CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimePtr trigger_pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>(name, "CCMCalibration", "NIMPulses", "CCMGeometry");
                current_noise_segment.pulses = trigger_pulses_view->Apply(*frame);
                extent_pulse_series_name = "DummyTriggerTimePulsesView";
                extent_frame->Put(extent_pulse_series_name, trigger_pulses_view);
            }
        } else {
            if(not frame->Has(input_noise_pulses_name_)) {
                if(do_warn)
                    log_warn("Frame does not have \"%s\", skipping frame.", input_noise_pulses_name_.c_str());
                ++n_skipped;
                continue;
            }
            current_noise_segment.pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_noise_pulses_name_);
            extent_pulse_series_name = input_noise_pulses_name_;
        }

        break;
    }

    std::tuple<double, double, std::map<CCMPMTKey, double>> noise_extent = GetTriggerExtent(extent_frame, extent_pulse_series_name, reco_bin_width_, false);
    current_noise_segment.start_time = std::get<0>(noise_extent);
    current_noise_segment.end_time = std::get<1>(noise_extent);
    current_noise_segment.offsets = std::get<2>(noise_extent);

    if(restrict_range_) {
        current_noise_segment.start_time = std::max(current_noise_segment.start_time, min_noise_time_);
        current_noise_segment.end_time = std::min(current_noise_segment.end_time, max_noise_time_);
        if(current_noise_segment.end_time <= current_noise_segment.start_time) {
            log_warn("Restricted noise pulse time range [%.2f, %.2f] is invalid. This occurred because the allowed time range is [%.2f, %.2f] while the time range for this pulse series is [%.2f, %.2f]", current_noise_segment.start_time, current_noise_segment.end_time, min_noise_time_, max_noise_time_, std::get<0>(noise_extent), std::get<1>(noise_extent));
        }
    }
    return failed;
}

void OverlayNoise::PulsesFromSegments(std::vector<PulsesMapSegment> const & segments, std::map<CCMPMTKey, double> const & target_offsets, std::map<CCMPMTKey, CCMRecoPulseSeries> & output, bool zero_out_width) {
    for(PulsesMapSegment const & segment : segments) {
        if(segment.pulses == nullptr) {
            continue;
        }
        double shift = segment.target_start_time - segment.start_time;
        for(auto const & [key, series] : *segment.pulses) {
            if(series.empty())
                continue;
            std::map<CCMPMTKey, double>::const_iterator target_it = target_offsets.find(key);
            std::map<CCMPMTKey, double>::const_iterator segment_it = segment.offsets.find(key);
            CCMRecoPulseSeries & output_series = output[key];
            for(CCMRecoPulse const & pulse : series) {
                if(pulse.GetTime() < segment.start_time)
                    continue;
                if(pulse.GetTime() >= segment.end_time)
                    break;
                output_series.push_back(pulse);
            }
            SnapPulsesToMod(output_series, shift, segment_it->second, target_it->second, reco_bin_width_);
        }
    }
}

std::tuple<std::vector<PulsesMapSegment>, std::map<CCMPMTKey, double>> OverlayNoise::GetNoisePulses(double duration, double zero_time, bool zero_out_width) {
    zero_out_width = not zero_out_width;

    if(not allow_stitching_) {
        double remaining_time = current_noise_segment.duration();
        bool already_looped = false;
        size_t n_attempts = 0;
        while(remaining_time < duration) {
            if(max_noise_time_ - min_noise_time_ < duration) {
                log_fatal("Restricted noise pulse time range [%.2f, %.2f] is smaller than the requested duration of %.2f ns.", min_noise_time_, max_noise_time_, duration);
            }
            if(n_attempts >= 100) {
                log_warn("Could not find enough noise pulses to cover the duration of %.2f ns without stitching after %lu attempts, trying again anyway.", duration, n_attempts);
                log_warn("Current noise pulse time range is [%.2f, %.2f] but we need %.2f ns", current_noise_segment.start_time, current_noise_segment.end_time, duration);
            }
            bool looped = NextFrame();
            if(already_looped and looped) {
                log_fatal("Not enough noise pulses to cover the duration of %.2f ns [%.2f, %.2f] without stitching.", duration, current_noise_segment.start_time, current_noise_segment.end_time);
            }
            already_looped |= looped;
            remaining_time = current_noise_segment.duration();
            ++n_attempts;
        }

        // Make the resulting pulses
        std::vector<PulsesMapSegment> segments; segments.reserve(1);
        segments.emplace_back(current_noise_segment);
        PulsesMapSegment & segment = segments.back();
        double random_start_time = randomService_->Uniform(
                current_noise_segment.start_time,
                std::max(current_noise_segment.start_time, current_noise_segment.end_time - duration)
        );
        segment.start_time = random_start_time;
        segment.end_time = segment.start_time + duration;
        segment.target_start_time = zero_time;

        // Update the current segment start time
        current_noise_segment.start_time = random_start_time;

        return {segments, current_noise_segment.offsets};
    }

    // Handle the stitching case

    double accumulated_time = 0.0;

    std::map<CCMPMTKey, double> target_board_time_offsets;
    if(current_noise_segment.pulses != nullptr) {
        target_board_time_offsets = current_noise_segment.offsets;
    }

    std::map<CCMPMTKey, double> board_offset_corrections;

    std::vector<PulsesMapSegment> segments;

    while(accumulated_time < duration) {
        if(current_noise_segment.pulses == nullptr) {
            NextFrame();
            target_board_time_offsets = current_noise_segment.offsets;
        }
        double remaining_time = current_noise_segment.end_time - current_noise_segment.start_time;
        double random_start_time = randomService_->Uniform(
                current_noise_segment.start_time,
                std::max(current_noise_segment.start_time, current_noise_segment.end_time - (duration - accumulated_time))
        );
        double random_end_time = std::min(current_noise_segment.end_time, random_start_time + std::min(remaining_time, duration - accumulated_time));
        segments.emplace_back(current_noise_segment);
        PulsesMapSegment & segment = segments.back();
        segment.start_time = random_start_time;
        segment.end_time = random_end_time;
        segment.target_start_time = accumulated_time + zero_time;

        accumulated_time += random_end_time - random_start_time;
        current_noise_segment.start_time = random_end_time;
        if(random_end_time >= current_noise_segment.end_time) {
            NextFrame();
        }
    }

    return {segments, target_board_time_offsets};
}

void OverlayNoise::DAQ(I3FramePtr frame) {

    // read in our reco pulse series
    boost::shared_ptr<CCMRecoPulseSeriesMap const> input_reco_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_sim_pulses_name_);
    I3MapPMTKeyDoubleConstPtr sim_board_time_offsets_ptr = frame->Get<I3MapPMTKeyDoubleConstPtr>(input_sim_offsets_name_);

    std::tuple<double, double> event_time_range = GetPulsesTimeRange(input_reco_pulses, true);

    double event_start_time = std::get<0>(event_time_range);
    double event_end_time = std::get<1>(event_time_range);
    double event_duration = event_end_time - event_start_time;

    double noise_pre_padding = randomService_->Uniform(min_padding_time_, max_padding_time_);
    double noise_post_padding = randomService_->Uniform(min_padding_time_, max_padding_time_);
    double noise_duration = event_duration + noise_pre_padding + noise_post_padding;

    // Get noise pulses with the correct time offset
    std::vector<PulsesMapSegment> noise_segments;
    std::map<CCMPMTKey, double> data_board_time_offsets;
    std::tie(noise_segments, data_board_time_offsets) = GetNoisePulses(noise_duration, event_start_time, true);

    boost::shared_ptr<CCMRecoPulseSeriesMap> noise_pulses = boost::make_shared<CCMRecoPulseSeriesMap>();
    if(match_times_to_simulation_) {
        // Just shift the noise pulse series before storing it
        PulsesFromSegments(noise_segments, *sim_board_time_offsets_ptr, *noise_pulses, true);
    } else {
        PulsesFromSegments(noise_segments, data_board_time_offsets, *noise_pulses, true);
    }

    // Store noise pulses in the frame
    frame->Put(input_noise_pulses_name_, noise_pulses);

    std::string sim_pulses_key_for_union = input_sim_pulses_name_;

    if(match_times_to_data_) {
        // Store the offsets for the simulation so they can be applied later
        std::string const target_offsets_key = input_noise_pulses_name_ + "_BoardTimeOffsets";

        // Calculate offsets and store them in the frame
        I3MapPMTKeyDoublePtr board_offset_corrections = boost::make_shared<I3MapPMTKeyDouble>();
        *(static_cast<std::map<CCMPMTKey, double>*>(board_offset_corrections.get())) =
            GetBoardOffsetCorrections(*sim_board_time_offsets_ptr, data_board_time_offsets, reco_bin_width_);
        frame->Put(target_offsets_key, board_offset_corrections);

        sim_pulses_key_for_union = input_sim_pulses_name_ + "_OffsetView";

        // Create offset view of the simulation pulses and store it in the frame
        CCMRecoPulseSeriesMapApplyOffsetsPtr offset_view = boost::make_shared<CCMRecoPulseSeriesMapApplyOffsets>(input_sim_pulses_name_, input_sim_offsets_name_, target_offsets_key, reco_bin_width_);
        frame->Put(sim_pulses_key_for_union, offset_view);
    }

    // Create the union of noise+simulation pulses and store it in the frame
    std::vector<std::string> union_inputs = {sim_pulses_key_for_union, input_noise_pulses_name_};
    CCMRecoPulseSeriesMapUnionPtr union_view = boost::make_shared<CCMRecoPulseSeriesMapUnion>(*frame, union_inputs);
    frame->Put(output_pulses_name_, union_view);

    PushFrame(frame);
}

