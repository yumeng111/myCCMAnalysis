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
#include <dataclasses/CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <phys-services/I3GSLRandomService.h>
#include <dataio/I3FrameSequence.h>

class OverlayNoise: public I3Module {
    std::string input_reco_pulse_name_;
    std::string output_reco_pulse_name_;
    std::string noise_pulse_name_;

    bool restrict_range_ = false;
    double min_noise_time_ = -10000;
    double max_noise_time_ = 0;
    bool allow_stitching_ = true;
    double padding_time_ = 500;

    std::string randomServiceName_;
    I3RandomServicePtr randomService_;
    dataio::I3FrameSequencePtr frame_sequences;
    std::deque<I3FramePtr> frame_cache;
    std::vector<std::string> file_lists;
public:
    OverlayNoise(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    bool PopFrame();
    void AddViews(I3FramePtr frame);
    void Finish();
};

I3_MODULE(OverlayNoise);

OverlayNoise::OverlayNoise(const I3Context& context) : I3Module(context),
    input_reco_pulse_name_("MCRecoPulses"), output_reco_pulse_name_("MCRecoPulsesPlusNoise"), noise_pulse_name_("TriggerTimePulses"),
    AddParameter("InputRecoPulseName", "", input_reco_pulse_name_);
    AddParameter("OutputRecoPulseName", "", output_reco_pulse_name_);
    AddParameter("NoisePulseName", "", noise_pulse_name_);

    AddParameter("RestrictRange", "Restrict the range of noise pulses to a specific time range", restrict_range_);
    AddParameter("MinNoiseTime", "Minimum time for noise pulses", min_noise_time_);
    AddParameter("MaxNoiseTime", "Maximum time for noise pulses", max_noise_time_);
    AddParameter("AllowStitching", "Allow noise pulses to be stitched together", allow_stitching_);
    AddParameter("PaddingTime", "Padding time for stitching noise pulses", padding_time_);
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
        return {false, mask->GetPulsesSource()};
    } else if(trigger_pulses != nullptr) {
        CCMRecoPulseSeriesMapApplySPECalPlusTriggerTimeConstPtr pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>(source_name, trigger_pulses->GetCalibrationName(), trigger_pulses->GetNIMPulsesName(), trigger_pulses->GetGeometryName());
        frame->Put(dest_name, pulses_view);
        return {true, trigger_pulses->GetPulsesSource()};
    } else if(beam_pulses != nullptr) {
        CCMRecoPulseSeriesMapApplySPECalPlusBeamTimeConstPtr pulses_view = boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusBeamTime>(source_name, beam_pulses->GetCalibrationName(), beam_pulses->GetNIMPulsesName(), beam_pulses->GetGeometryName());
        frame->Put(dest_name, pulses_view);
        return {true, beam_pulses->GetPulsesSource()};
    } else {
        return {false, ""};
    }
}

std::tuple<double, double> GetTimeRange(I3FramePtr frame, std::string pulse_series_name, bool expansive = false) {
    if(not frame->Has(pulse_name)) {
        log_fatal("Frame does not have %s", pulse_name.c_str());
    }
    if(not frame->Has("CCMDAQConfig")) {
        log_fatal("Frame does not have CCMDAQConfig");
    }

    I3FramePtr dummy_frame = boost::make_shared<I3Frame>(*frame);
    size_t dummy_number = 0;
    std::string next_pulse_series_name = pulse_series_name;

    // For each layer of pulse series correction we need to put a corresponding correction in place for our dummy pulse series
    // The final result will be read out from DummyPulses0
    while(not IsBasicPulseSeries(next_pulse_series_name)) {
        std::tuple<bool, std::string> x = AddDummyPulseSeries(dummy_frame, next_pulse_series_name, dummy_number);
        bool deeper = std::get<0>(x);
        next_pulse_series_name = std::get<1>(x);
        if(deeper) {
            ++dummy_number;
        }
    }
    std::string const & source_pulse_series_name = next_pulse_series_name;

    CCMDAQConfigConstPtr daq_config = frame->Get<CCMDAQConfigConstPtr>("CCMDAQConfig");

    // Prepare the dummy input pulse series
    CCMRecoPulseSeriesMapPtr dummy_pulses = boost::make_shared<CCMRecoPulseSeriesMap>();
    CCMRecoPulse first_pulse;
    CCMRecoPulse last_pulse;
    first_pulse.time = 0
    last_pulse.time = daq_config->machine_configurations[0].num_samples * 2
    CCMRecoPulseSeries dummy_series{first_pulse, last_pulse};

    CCMRecopulseSeriesMapConstPtr pulses = dummy_frame->Get<CCMRecoPulseSeriesMapConstPtr>(source_pulse_series_name);
    for(CCMRecoPulseSeriesMap::const_iterator it = pulses->begin(); it != pulses->end(); ++it) {
        dummy_pulses->insert(std::make_pair(it->first, dummy_series));
    }

    std::stringstream ss;
    ss << "DummyPulses" << dummy_number;
    std::string dummy_pulse_series_name = ss.str();

    // Put the dummy pulse series in the frame
    dummy_frame->Put(dummy_pulse_series_name, dummy_pulses);

    CCMRecopulseSeriesMapConstPtr corrected_pulses = dummy_frame->Get<CCMRecoPulseSeriesMapConstPtr>("DummyPulses0");

    double min_time;
    double max_time;

    if(expansive) {
        min_time = std::numeric_limits<double>::max();
        max_time = std::numeric_limits<double>::min();

        for(CCMRecoPulseSeriesMap::const_iterator it = corrected_pulses->begin(); it != corrected_pulses->end(); ++it) {
            double const & start_time = it->second[0].time;
            double const & end_time = it->second[1].time;
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
        for(CCMRecoPulseSeriesMap::const_iterator it = corrected_pulses->begin(); it != corrected_pulses->end(); ++it) {
            double const & start_time = it->second[0].time;
            double const & end_time = it->second[1].time;
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

bool OverlayNoise::PopFrame(){
    I3FramePtr frame;
    // Grab frames until we get a DAQ frame
    while(true) {
        // Fail if we do not have any more frames
        if(not frame_sequences->more())
            return false;
        frame = frame_sequences->pop_frame();
        if(frame->GetStop() != I3Frame::DAQ) {
            // Skip non-DAQ frames
            continue;
        }
        break;
    }

    // Store the frame
    frame_cache.push_back(frame);
    return true;
}

void OverlayNoise::AddViews(I3FramePtr frame){
    // Keys should be an argument
    // Also should have an option for the BeamTimePulses
    CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime corrected_pulses("WavedeformPulses", "CCMCalibration", "NIMPulsesMode", "CCMGeometry");
    frame->Put(noise_pulse_name_, boost::make_shared<CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime>(corrected_pulses));
}

void OverlayNoise::DAQ(I3FramePtr frame){

    // read in our reco pulse series
    boost::shared_ptr<CCMRecoPulseSeriesMap const> input_reco_pulses = frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(input_reco_pulse_name_);

    // set up reco pulses series to hold MC merged with noise
    CCMRecoPulseSeriesMapPtr mc_reco_overlay_noise_dest = boost::make_shared<CCMRecoPulseSeriesMap>();

    // grab a noise frame
    bool pop_frame_success = PopFrame();
    while (not pop_frame_success){
        // oops weve run out of frames! let's reset frame sequences and run over files again
        frame_sequences = boost::make_shared<dataio::I3FrameSequence>(file_lists);
        pop_frame_success = PopFrame();
    }
    I3FramePtr noise_frame = frame_cache[0];

    // check for our noise pulses in the frame
    if (!noise_frame->Has(noise_pulse_name_)){
        // Make this optional or add some logic to do this conditionally
        AddViews(noise_frame);
    }
    boost::shared_ptr<CCMRecoPulseSeriesMap const> noise_reco_pulses = noise_frame->Get<boost::shared_ptr<CCMRecoPulseSeriesMap const>>(noise_pulse_name_);

    // great, now let's overlay MC and noise
    // need random number to get the start time
    double noise_start_time = randomService_->Uniform(noise_region_start_time_, noise_region_end_time_);
    double noise_end_time = noise_start_time + noise_duration_;

    // now let's merge pulse series
    for (CCMRecoPulseSeriesMap::const_iterator it = input_reco_pulses->begin(); it != input_reco_pulses->end(); ++it) {

        // Find the corresponding PMT in the destination map
        CCMRecoPulseSeriesMap::iterator it_dest = mc_reco_overlay_noise_dest->find(it->first);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == mc_reco_overlay_noise_dest->end()) {
            mc_reco_overlay_noise_dest->insert(std::make_pair(it->first, CCMRecoPulseSeries()));
            // Update the iterator so it points to our new entry
            it_dest = mc_reco_overlay_noise_dest->find(it->first);
        }

        // Reference to the destination
        CCMRecoPulseSeries & dest_series = it_dest->second;

        // Iterate over the vector of CCMRecoPulse in the source map for this PMT
        for (CCMRecoPulseSeries::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            // save to dest series
            double new_time = it2->GetTime();
            double new_charge = it2->GetCharge();
            double new_width = it2->GetWidth();
            CCMRecoPulse new_pulse;
            new_pulse.SetCharge(new_charge);
            new_pulse.SetTime(new_time);
            new_pulse.SetWidth(new_width);
            dest_series.push_back(new_pulse);
        }

        // Check if the same PMT key exists in noise pulses
        CCMRecoPulseSeriesMap::const_iterator it_noise = noise_reco_pulses->find(it->first);
        if (it_noise != noise_reco_pulses->end()) {
            // Copy noise pulses
            for (CCMRecoPulseSeries::const_iterator it3 = it_noise->second.begin(); it3 != it_noise->second.end(); ++it3) {
                // check time time of this pulse
                if (it3->GetTime() >= noise_start_time and it3->GetTime() < noise_end_time){
                    // great we want to save this noise pulse!
                    double noise_time = it3->GetTime();
                    double noise_charge = it3->GetCharge();
                    noise_time -= noise_start_time;
                    CCMRecoPulse noise_pulse;
                    noise_pulse.SetCharge(noise_charge);
                    noise_pulse.SetTime(noise_time);
                    noise_pulse.SetWidth(0.0); // using width to track cherenkov light
                    dest_series.push_back(noise_pulse);
                }
            }
        }

        // Done combining pulse series! Let's sort according to time
        std::sort(dest_series.begin(), dest_series.end(), [](const CCMRecoPulse& a, const CCMRecoPulse& b) { return a.GetTime() < b.GetTime(); });
    }

    // remove frame from cache
    frame_cache.pop_front();

    // Save
    frame->Put(output_reco_pulse_name_, mc_reco_overlay_noise_dest);
    PushFrame(frame);

}

void OverlayNoise::Finish() {}

