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
    std::string input_reco_pulse_name_;
    std::string output_reco_pulse_name_;
    std::string noise_pulse_name_;

    bool restrict_range_ = false;
    double min_noise_time_ = -10000;
    double max_noise_time_ = 0;
    bool allow_stitching_ = true;
    double min_padding_time_ = 0;
    double max_padding_time_ = 500;

    std::string randomServiceName_;
    I3RandomServicePtr randomService_;
    dataio::I3FrameSequencePtr frame_sequences;
    std::deque<I3FramePtr> frame_cache;
    std::vector<std::string> file_lists;

    boost::shared_ptr<CCMRecoPulseSeriesMap const> current_noise_pulses = nullptr;
    double current_start_time;
    double current_end_time;
public:
    OverlayNoise(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    bool NextFrame();
    boost::shared_ptr<CCMRecoPulseSeriesMap> GetNoisePulses(double duration, double zero_time, bool zero_out_width = false);
};

I3_MODULE(OverlayNoise);

OverlayNoise::OverlayNoise(const I3Context& context) : I3Module(context),
    input_reco_pulse_name_("MCRecoPulses"), output_reco_pulse_name_("MCRecoPulsesPlusNoise"), noise_pulse_name_("TriggerTimePulses") {
    AddParameter("InputRecoPulseName", "", input_reco_pulse_name_);
    AddParameter("OutputRecoPulseName", "", output_reco_pulse_name_);
    AddParameter("NoisePulseName", "", noise_pulse_name_);

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

std::tuple<double, double> GetPulsesTimeRange(CCMRecoPulseSeriesMapConstPtr pulses, bool expansive = false) {
    double min_time;
    double max_time;

    if(expansive) {
        min_time = std::numeric_limits<double>::max();
        max_time = std::numeric_limits<double>::min();

        for(CCMRecoPulseSeriesMap::const_iterator it = pulses->begin(); it != pulses->end(); ++it) {
            if(it->second.empty()) {
                continue;
            }
            double const & start_time = it->second.front().GetTime();
            double const & end_time = it->second.back().GetTime();
            if(start_time < min_time) {
                min_time = start_time;
            }
            if(end_time > max_time) {
                max_time = end_time;
            }
        }
    } else {
        min_time = std::numeric_limits<double>::min();
        max_time = std::numeric_limits<double>::max();
        for(CCMRecoPulseSeriesMap::const_iterator it = pulses->begin(); it != pulses->end(); ++it) {
            if(it->second.empty()) {
                continue;
            }
            double const & start_time = it->second.front().GetTime();
            double const & end_time = it->second.back().GetTime();
            if(start_time > min_time) {
                min_time = start_time;
            }
            if(end_time < max_time) {
                max_time = end_time;
            }
        }
    }
    return {min_time, max_time};
}

std::tuple<double, double> GetTriggerExtent(I3FramePtr frame, std::string pulse_series_name, bool expansive = false) {
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
    for(CCMRecoPulseSeriesMap::const_iterator it = pulses->begin(); it != pulses->end(); ++it) {
        dummy_pulses->insert(std::make_pair(it->first, dummy_series));
    }

    std::stringstream ss;
    ss << "DummyPulses" << dummy_number;
    std::string dummy_pulse_series_name = ss.str();

    // Put the dummy pulse series in the frame
    dummy_frame->Put(dummy_pulse_series_name, dummy_pulses);

    CCMRecoPulseSeriesMapConstPtr corrected_pulses = dummy_frame->Get<CCMRecoPulseSeriesMapConstPtr>("DummyPulses0");

    return GetPulsesTimeRange(corrected_pulses, expansive);
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
        }
        frame = frame_sequences->pop_frame();
        if(frame->GetStop() != I3Frame::DAQ) {
            // Skip non-DAQ frames
            continue;
        }
        break;
    }
    current_noise_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);
    std::tuple<double, double> noise_extent = GetTriggerExtent(frame, noise_pulse_name_, false);
    current_start_time = std::get<0>(noise_extent);
    current_end_time = std::get<1>(noise_extent);

    if(restrict_range_) {
        current_start_time = std::max(current_start_time, min_noise_time_);
        current_end_time = std::min(current_end_time, max_noise_time_);
    }
    return failed;
}

boost::shared_ptr<CCMRecoPulseSeriesMap> OverlayNoise::GetNoisePulses(double duration, double zero_time, bool zero_out_width) {
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
        for(CCMRecoPulseSeriesMap::const_iterator it = current_noise_pulses->begin(); it != current_noise_pulses->end(); ++it) {
            CCMRecoPulseSeries noise_series;
            for(CCMRecoPulseSeries::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                double time = it2->GetTime();
                if(time >= random_start_time and time < random_end_time + duration) {
                    noise_series.push_back(*it2);
                    CCMRecoPulse & noise = noise_series.back();
                    noise.SetTime(time - random_start_time - zero_time);
                    noise.SetWidth(noise.GetWidth() * zero_out_width);
                }
            }
            noise_pulses->insert(std::make_pair(it->first, noise_series));
        }
        current_start_time = random_start_time;
        return noise_pulses;
    }

    double accumulated_time = 0.0;
    boost::shared_ptr<CCMRecoPulseSeriesMap> noise_pulses = boost::make_shared<CCMRecoPulseSeriesMap>();
    while(accumulated_time < duration) {
        if(current_noise_pulses == nullptr) {
            NextFrame();
        }
        double remaining_time = current_end_time - current_start_time;
        double random_start_time = randomService_->Uniform(
                current_start_time,
                std::max(current_start_time, current_end_time - (duration - accumulated_time))
        );
        double random_end_time = random_start_time + std::min(remaining_time, duration - accumulated_time);
        for(CCMRecoPulseSeriesMap::const_iterator it = current_noise_pulses->begin(); it != current_noise_pulses->end(); ++it) {
            if(not noise_pulses->count(it->first)) {
                noise_pulses->insert(std::make_pair(it->first, CCMRecoPulseSeries()));
            }
            CCMRecoPulseSeries & noise_series = noise_pulses->at(it->first);
            for(CCMRecoPulseSeries::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                double time = it2->GetTime();
                if(time >= random_start_time and time < random_end_time) {
                    noise_series.push_back(*it2);
                    noise_series.back().SetTime(time - random_start_time + accumulated_time - zero_time);
                }
            }
            noise_pulses->insert(std::make_pair(it->first, noise_series));
        }
        accumulated_time += random_end_time - random_start_time;
        current_start_time = random_end_time;
        if(random_end_time >= current_end_time) {
            NextFrame();
        }
    }
    return noise_pulses;
}

void OverlayNoise::DAQ(I3FramePtr frame) {

    // read in our reco pulse series
    boost::shared_ptr<CCMRecoPulseSeriesMap const> input_reco_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_reco_pulse_name_);

    std::tuple<double, double> event_time_range = GetPulsesTimeRange(input_reco_pulses, true);

    double event_start_time = std::get<0>(event_time_range);
    double event_end_time = std::get<1>(event_time_range);
    double event_duration = event_end_time - event_start_time;

    double noise_pre_padding = randomService_->Uniform(min_padding_time_, max_padding_time_);
    double noise_post_padding = randomService_->Uniform(min_padding_time_, max_padding_time_);
    double noise_duration = event_duration + noise_pre_padding + noise_post_padding;

    // Get noise pulses with the correct time offset
    CCMRecoPulseSeriesMapPtr combined_pulses = GetNoisePulses(noise_duration, event_start_time, true);

    // Copy the event pulses into the destination map
    for(CCMRecoPulseSeriesMap::const_iterator it = input_reco_pulses->begin(); it != input_reco_pulses->end(); ++it) {

        // Find the corresponding PMT in the destination map
        CCMRecoPulseSeriesMap::iterator it_dest = combined_pulses->find(it->first);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == combined_pulses->end()) {
            combined_pulses->insert(std::make_pair(it->first, CCMRecoPulseSeries()));
            // Update the iterator so it points to our new entry
            it_dest = combined_pulses->find(it->first);
        }

        // Reference to the destination
        CCMRecoPulseSeries & dest_series = it_dest->second;

        // Iterate over the vector of CCMRecoPulse in the source map for this PMT
        for(CCMRecoPulseSeries::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            dest_series.push_back(*it2);

        // Done combining pulse series! Let's sort according to time
        std::sort(dest_series.begin(), dest_series.end(), [](const CCMRecoPulse& a, const CCMRecoPulse& b) { return a.GetTime() < b.GetTime(); });
    }

    // Save
    //TODO save the noise series
    //TODO Add a charge cut on stitched frames
    frame->Put(output_reco_pulse_name_, combined_pulses);
    PushFrame(frame);

}

