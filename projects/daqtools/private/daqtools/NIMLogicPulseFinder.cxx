#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <sys/time.h>
#include <sys/resource.h>
#include <string>
#include <iostream>
#include <algorithm>

#include <icetray/I3Frame.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

class NIMLogicPulseFinder : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;
    std::string nim_name_;

    size_t n_frames_for_baseline;
    size_t sample_step;
    double constant_fraction;
    double minimum_nim_pulse_height;

    // Internal state
    bool geo_seen;
    I3Map<CCMTriggerKey, uint32_t> trigger_channel_map_;
    bool estimate_baselines;
    std::deque<I3FramePtr> cached_frames;
    std::map<CCMTriggerKey, std::vector<uint16_t>> baseline_samples;
    std::map<CCMTriggerKey, double> baselines;
    std::map<CCMTriggerKey, double> baseline_stddevs;

    // Constants
    constexpr static double ns_per_sample = 2.0;
    double exp_smoothing_tau_ = 10.0;

public:
    NIMLogicPulseFinder(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    std::vector<double> SmoothWaveform(std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end);
    void AddBaselineSamples(I3FramePtr frame);
    void EstimateBaselines();
    double QuickBaselineEstimate(CCMWaveformUInt16 const & waveform, double baseline, double baseline_stddev);
    bool CheckBaselines(I3FramePtr frame);
    NIMLogicPulseSeries GetNIMPulses(CCMWaveformUInt16 const & waveform, double baseline, double baseline_stddev);
    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(NIMLogicPulseFinder);

NIMLogicPulseFinder::NIMLogicPulseFinder(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_(""), nim_name_(""), geo_seen(false), estimate_baselines(true) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("SampleStep", "The number of steps between samples used for the initial baseline estimate.", size_t(50));
    AddParameter("NFramesForBaseline", "The number of frames to use for a baseline estimate", size_t(10));
    AddParameter("ConstantFraction", "The fraction of the pulse height to use for its start time", double(0.01));
    AddParameter("MinimumPulseHeight", "The minimum pulse height to consider a NIM pulse to be present.", double(1000));
    AddParameter("NIMLogicPulseSeriesMapName", "Name for the output nim pulses map", std::string("NIMPulses"));
}

void NIMLogicPulseFinder::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("SampleStep", sample_step);
    GetParameter("NFramesForBaseline", n_frames_for_baseline);
    GetParameter("ConstantFraction", constant_fraction);
    GetParameter("MinimumPulseHeight", minimum_nim_pulse_height);
    GetParameter("NIMLogicPulseSeriesMapName", nim_name_);
}

void NIMLogicPulseFinder::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    // Cache the trigger channel map
    trigger_channel_map_ = geo.trigger_channel_map;

    // Clear out the stored baseline values
    baseline_samples.clear();
    baselines.clear();

    // Initialize the baseline values
    for(std::pair<CCMTriggerKey, uint32_t> const & key : trigger_channel_map_) {
        baseline_samples.insert({key.first, std::vector<uint16_t>()});
        baselines.insert({key.first, 0});
        baseline_stddevs.insert({key.first, 0});
    }
    geo_seen = true;

    PushFrame(frame);
}

std::vector<double> NIMLogicPulseFinder::SmoothWaveform(std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end) {
    // Invert waveform
    // Smooth with exponential filter
    // Smooth with box filter

    size_t N = std::distance(begin, end);

    std::vector<double> smoothed_wf(N);

    // Average the first N samples to get the starting value for the exponential smoothing
    int init_avg_end_idx = std::min(3 * std::max(int(exp_smoothing_tau_ / ns_per_sample), 1), int(N));
    double exp_start = 0.0;
    std::vector<uint16_t>::const_iterator x_it = begin;
    for(size_t i=0; i < init_avg_end_idx and x_it != end; ++i, ++x_it) {
        exp_start += -double(*x_it);
    }
    exp_start /= init_avg_end_idx;

    double alpha = exp(-ns_per_sample / exp_smoothing_tau_);
    double y_i = exp_start;
    x_it = begin;
    for(size_t i=0; i < smoothed_wf.size() and x_it != end; ++i, ++x_it) {
        double x = -double(*x_it);
        y_i += (1.0 - alpha) * (x - y_i);
        smoothed_wf[i] = y_i;
    }

    // Box smoothing for the first index
    double smoothed_wf_first_value = (smoothed_wf[0] + smoothed_wf[1]) / 2.0;

    // Box smoothing for the last index
    double smoothed_wf_last_value = (smoothed_wf[smoothed_wf.size() - 2] + smoothed_wf[smoothed_wf.size() - 1]) / 2.0;

    // Box smoothing for everything in between
    for(size_t i=1; i < smoothed_wf.size() - 1; ++i) {
        double prev = smoothed_wf[i - 1];
        double current = smoothed_wf[i];
        double next = smoothed_wf[i + 1];
        smoothed_wf[i] = (prev + current + next) / 3.0;
    }
    smoothed_wf[0] = smoothed_wf_first_value;
    smoothed_wf[smoothed_wf.size() - 1] = smoothed_wf_last_value;

    return smoothed_wf;
}

void NIMLogicPulseFinder::AddBaselineSamples(I3FramePtr frame) {
    I3Vector<CCMWaveformUInt16> const & waveforms = frame->Get<I3Vector<CCMWaveformUInt16> const>(waveforms_name_);
    // Iterate over trigger channels
    for(std::pair<CCMTriggerKey, uint32_t> const key : trigger_channel_map_) {
        CCMTriggerKey trigger_key = key.first;
        uint32_t trigger_channel = key.second;
        std::vector<uint16_t> & samples = baseline_samples[trigger_key];
        // Concatenate waveforms from multiple frames
        samples.reserve(samples.size() + waveforms[trigger_channel].GetWaveform().size() / sample_step);
        // samples.insert(samples.end(), waveforms[trigger_channel].GetWaveform().begin(), waveforms[trigger_channel].GetWaveform().end());
        for(size_t i=0; i<waveforms[trigger_channel].GetWaveform().size(); i += sample_step) {
            samples.push_back(waveforms[trigger_channel].GetWaveform()[i]);
        }
    }
}

void NIMLogicPulseFinder::EstimateBaselines() {
    for(std::pair<CCMTriggerKey const, std::vector<uint16_t>> & key_val : baseline_samples) {
        CCMTriggerKey key = key_val.first;
        std::vector<uint16_t> & samples = key_val.second;
        if(samples.size() == 0)
            continue;
        std::vector<double> smoothed = SmoothWaveform(samples.begin(), samples.end());
        std::sort(smoothed.begin(), smoothed.end()); // Must be sorted for Mode function
        // Use the mode of N smoothed waveforms as a baseline estimate
        double baseline = robust_stats::Mode(smoothed.data(), smoothed.size());
        baselines[key] = baseline;
        // Estimate the stddev of the baseline
        double baseline_stddev = robust_stats::MedianAbsoluteDeviation(smoothed.begin(), smoothed.end(), baseline);
        baseline_stddevs[key] = baseline_stddev;
    }

    // Clear out the concatenated waveforms now that we don't need them
    baseline_samples.clear();
    // Indicate that the baselines have already been estimated
    estimate_baselines = false;
}

double NIMLogicPulseFinder::QuickBaselineEstimate(CCMWaveformUInt16 const & waveform, double baseline, double baseline_stddev) {
    // Make a quick estimate of the baseline using a subset of the waveform
    std::vector<double> samples;
    std::vector<uint16_t> const & wf = waveform.GetWaveform();
    samples.reserve(wf.size() / sample_step);
    for(size_t i=0; i<wf.size(); i += sample_step) {
        if(std::abs(-wf[i] - baseline) >= 100 + 3 * baseline_stddev)
            continue;
        samples.push_back(-wf[i]);
    }
    std::sort(samples.begin(), samples.end());
    return robust_stats::Mode(samples.data(), samples.size());
}

NIMLogicPulseSeries NIMLogicPulseFinder::GetNIMPulses(CCMWaveformUInt16 const & waveform, double baseline, double baseline_stddev) {
    // Search for NIM pulses in the waveform
    //
    // First look for a sample that goes above the threshold
    //      Mark the point just before this as the beginning of the search region for the pulse edge
    //      Where "before" is defined as the last point where the waveform was below baseline+1sigma
    //      Start tracking samples to look for the maxiumum value of the pulse
    // Then look for a sample that goes below the threshold
    //      Mark the point just after this i.e. +5 time bins
    // From the peak, search backwards for the point where the waveform drops below the constant fraction
    // From the peak, search forwards for the point where the waveform drops below the constant fraction
    // These two points mark the beginning and end of the pulse
    // Now go back to step one with the next sample

    std::vector<uint16_t> const & wf = waveform.GetWaveform();

    NIMLogicPulseSeries pulse_series;

    bool in_pulse = false;
    size_t start_search_begin_idx = 0;
    double max_sample = 0;
    size_t max_sample_pos = 0;
    double begin_threshold = minimum_nim_pulse_height;
    //double end_threshold = std::min(baseline_stddev, minimum_nim_pulse_height);
    double end_threshold = minimum_nim_pulse_height;
    std::vector<double> nim_samples_for_mode;

    // Iterate over samples
    for(size_t i=0; i<wf.size(); ++i) {
        // Invert and subtract baseline
        double sample = -wf[i] - baseline;
        if(in_pulse) {
            // Update maximum value in the pulse
            if(sample > max_sample) {
                max_sample = sample;
                max_sample_pos = i;
            }
            // Save sample for mode
            nim_samples_for_mode.push_back(sample);
            // Passing below threshold means we reached the end of the pulse
            if(sample <= end_threshold) {
                // Let's calculate the mode
                std::sort(nim_samples_for_mode.begin(), nim_samples_for_mode.end());
                double this_mode = robust_stats::Mode(nim_samples_for_mode.data(), nim_samples_for_mode.size());
                size_t start_search_end_idx = max_sample_pos;
                size_t end_search_begin_idx = max_sample_pos;
                size_t end_search_end_idx = std::min(i + 5, wf.size());
                size_t pulse_start_pos = start_search_begin_idx;
                size_t pulse_end_pos = start_search_begin_idx;
                //double threshold = max_sample * constant_fraction;
                double threshold = this_mode * constant_fraction;
                double pulse_start_time = 0.0;
                // Search for the beginning of the pulse
                for(size_t j=start_search_end_idx-1; j>=start_search_begin_idx; --j) {
                    sample = -wf[j] - baseline;
                    if(sample <= threshold) {
                        pulse_start_pos = j;
                        // let's also interpolate to get our start time
                        double time_below_threshold = j * ns_per_sample;
                        double sample_below_threshold = -wf[j] - baseline;
                        double time_above_threshold = (j+1) * ns_per_sample;
                        double sample_above_threshold = -wf[j+1] - baseline;
                        pulse_start_time = time_below_threshold + ((threshold - sample_below_threshold) * ((time_above_threshold - time_below_threshold) / (sample_above_threshold - sample_below_threshold)));
                        break;
                    }
                }
                double pulse_end_time = 0.0;
                // Search for the end of the pulse
                for(size_t j=end_search_begin_idx; j<end_search_end_idx; ++j) {
                    sample = -wf[j] - baseline;
                    if(sample <= threshold) {
                        pulse_end_pos = j;
                        // let's also interpolate to get our end time
                        double time_below_threshold = j * ns_per_sample;
                        double sample_below_threshold = -wf[j] - baseline;
                        double time_above_threshold = (j-1) * ns_per_sample;
                        double sample_above_threshold = -wf[j-1] - baseline;
                        pulse_end_time = time_below_threshold + ((threshold - sample_below_threshold) * ((time_above_threshold - time_below_threshold) / (sample_above_threshold - sample_below_threshold)));
                        break;
                    }
                }
                // Store the pulse
                NIMLogicPulse pulse;
                //pulse.SetNIMPulseTime(pulse_start_pos * ns_per_sample);
                pulse.SetNIMPulseTime(pulse_start_time);
                pulse.SetNIMPulseLength(pulse_end_time - pulse_start_time);
                pulse_series.push_back(pulse);

                // We're done with the current pulse
                in_pulse = false;
            }
        } else {
            // Check for the beginning of the pulse
            if(sample >= minimum_nim_pulse_height) {
                max_sample = sample;
                max_sample_pos = i;
                for(size_t j=i; true; --j) {
                    double sample = -wf[j] - baseline;
                    if(sample <= baseline_stddev) {
                        start_search_begin_idx = j;
                        break;
                    }
                    if(j == 0)
                        break;
                }

                // We're now dealing with a pulse
                in_pulse = true;
                // Save samples to calculate mode
                nim_samples_for_mode.clear();
                nim_samples_for_mode.push_back(sample);
            }
        }
    }
    return pulse_series;
}

void NIMLogicPulseFinder::ProcessFrame(I3FramePtr frame) {
    // Find the NIM pulses for each trigger channel and store them in the frame
    I3Vector<CCMWaveformUInt16> const & waveforms = frame->Get<I3Vector<CCMWaveformUInt16> const>(waveforms_name_);
    NIMLogicPulseSeriesMapPtr pulse_series_map = boost::make_shared<NIMLogicPulseSeriesMap>();
    for(std::pair<CCMTriggerKey, uint32_t> const key : trigger_channel_map_) {
        CCMTriggerKey const & trigger_key = key.first;
        uint32_t const & trigger_channel = key.second;
        double baseline = baselines[trigger_key];
        double baseline_stddev = baseline_stddevs[trigger_key];
        NIMLogicPulseSeries nim_pulses = GetNIMPulses(waveforms.at(trigger_channel), baseline, baseline_stddev);
        pulse_series_map->insert({trigger_key, nim_pulses});
    }

    frame->Put(nim_name_, pulse_series_map);
}

bool NIMLogicPulseFinder::CheckBaselines(I3FramePtr frame) {
    // Check if a quick estimate of the baselines significantly deviates from the existing estimates
    I3Vector<CCMWaveformUInt16> const & waveforms = frame->Get<I3Vector<CCMWaveformUInt16> const>(waveforms_name_);
    for(std::pair<CCMTriggerKey, uint32_t> const key : trigger_channel_map_) {
        CCMTriggerKey const & trigger_key = key.first;
        uint32_t const & trigger_channel = key.second;
        double baseline = baselines[trigger_key];
        double baseline_stddev = baseline_stddevs[trigger_key];
        if(waveforms[trigger_channel].GetWaveform().size() > 0) {
            double baseline_estimate = QuickBaselineEstimate(waveforms[trigger_channel], baseline, baseline_stddev);
            if(std::abs(baseline - baseline_estimate) > 10.0 * baseline_stddev) {
                // Report that we need a new baseline estimate
                return true;
            }
        }
    }
    // Report that the current baseline estimates are fine
    return false;
}

void NIMLogicPulseFinder::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    // Check if we still need a baseline estimate
    if(estimate_baselines) {
        if(cached_frames.size() == n_frames_for_baseline) {
            EstimateBaselines();
            for(size_t i=0; i<cached_frames.size(); ++i) {
                ProcessFrame(cached_frames[i]);
                PushFrame(cached_frames[i]);
            }
            cached_frames.clear();
        } else {
            AddBaselineSamples(frame);
            cached_frames.push_back(frame);
            return;
        }
    }
    // Check if our baseline estimate is still valid
    if(CheckBaselines(frame)) {
        AddBaselineSamples(frame);
        cached_frames.push_back(frame);
        estimate_baselines = true;
        return;
    }
    ProcessFrame(frame);
    PushFrame(frame);
}

void NIMLogicPulseFinder::Finish() {
    // If we are currently accumulating frames for a new baseline estimate
    // then we need to finish the estimate and clear out all the frames before
    // the tray can finish processing
    if(cached_frames.size()) {
        EstimateBaselines();
        for(size_t i=0; i<cached_frames.size(); ++i) {
            ProcessFrame(cached_frames[i]);
            PushFrame(cached_frames[i]);
        }
        cached_frames.clear();
    }
    Flush();
}
