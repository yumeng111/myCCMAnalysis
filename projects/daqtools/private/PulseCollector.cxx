#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <tuple>
#include <string>
#include <random>
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

#include "daqtools/OnlineRobustStats.h"
#include "daqtools/WaveformSmoother.h"

namespace {

template <class RandomIt, class RandomFunc, class U = typename std::iterator_traits<RandomIt>::value_type>
std::vector<U> ChooseNRandom(size_t N, RandomIt begin, RandomIt end, RandomFunc&& g) {
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    typedef std::uniform_int_distribution<diff_t> distr_t;
    typedef typename distr_t::param_type param_t;
    distr_t D;
    diff_t n = std::distance(begin, end);
    diff_t max_n = std::min(diff_t(N), n - 1);
    std::vector<U> result(begin, end);
    for(diff_t i=0; i<max_n; ++i) {
        std::swap(result[i], result[D(g, param_t(i, n))]);
    }
    result.resize(N);
    return result;
}

size_t CheckForPulse(WaveformSmoother & smoother, size_t start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold, double baseline) {
    smoother.Reset(start_idx);
    // 0 before start of pulse checking [checking for positive derivative]
    // 1 derivative is positive (in rising edge of pulse) [checking for negative derivative]
    // 2 derivative is negative (in falling edge of pulse) [checking for positive derivative]
    // 3 derivative is positive (recovering from droop) [done]
    int state = 1;
    double max_value = smoother.Value() - baseline;
    double max_abs_derivative = std::abs(smoother.Derivative());
    double integral = 0;
    bool found_pulse = false;
    size_t final_idx = start_idx;
    size_t N = std::min(smoother.Size(), start_idx + max_samples);
    for(size_t i=start_idx; i<N; ++i) {
        max_value = std::max(max_value, smoother.Value() - baseline);
        max_abs_derivative = std::max(max_abs_derivative, std::abs(smoother.Derivative()));
        integral += (smoother.Value() - baseline);
        if(state % 2) {
            if(smoother.Derivative() < 0)
                state += 1;
        } else
            if(smoother.Derivative() > 0) {
                state += 1;
        }
        if(state >= 3) {
            found_pulse = true;
            final_idx = i;
            break;
        }
        smoother.Next();
    }
    smoother.Reset(start_idx);
    if(found_pulse
            and (max_value >= value_threshold)
            and (max_abs_derivative >= derivative_threshold)
            and (integral >= integral_threshold)
            and (final_idx - start_idx >= min_length)
            ) {
        return final_idx;
    } else {
        return start_idx;
    }
}

double EstimateBaseline(WaveformSmoother & smoother, size_t samples_for_baseline) {
    smoother.Reset();
    size_t N = smoother.Size();
    N = std::min(N, samples_for_baseline);
    std::vector<double> baseline_samples;
    baseline_samples.reserve(N);
    smoother.Reset();
    for(size_t i=0; i<N; ++i) {
        baseline_samples.push_back(smoother.Value());
        smoother.Next();
    }
    std::sort(baseline_samples.begin(), baseline_samples.end());
    double baseline = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());
    return baseline;
}

std::vector<std::pair<size_t, size_t>> FindPulses(WaveformSmoother & smoother, double baseline, size_t pulse_max_start_sample, double deriv_threshold, size_t max_pulse_width, size_t min_pulse_width, double min_pulse_height, double min_deriv_magnitude, double min_integral) {
    // Determine the maximum starting sample for a pulse
    // We need at least 5 samples after the starting pulse
    size_t N = smoother.Size();
    // Define the maximum sample to consider for the start of a pulse
    // A pulse needs at least 5 samples, so subtract off 5 from the max sample we can look at
    N = std::min(N, pulse_max_start_sample);
    if(N > 5)
        N -= 5;
    size_t pulse_first_index = 0;
    size_t pulse_last_index = 0;
    size_t nearby_charge_window = 300;
    std::vector<std::pair<size_t, size_t>> pulses;

    // Reset the smoother position to the start of the waveform
    smoother.Reset();
    for(size_t i=0; i<N; ++i) {
        // Check the condition for the beginning of a pulse
        if(smoother.Derivative() > deriv_threshold) {
            // Get the last index of the found pulse, 0 if no pulse is found
            pulse_last_index = CheckForPulse(smoother, i, max_pulse_width, min_pulse_width, min_pulse_height, min_deriv_magnitude, min_integral, baseline);
            if(pulse_last_index > i) {
                pulse_first_index = size_t(std::max(ptrdiff_t(0), ptrdiff_t(i) - 5));
                // we found a pulse! now let's do a quick check on charge within the nearby_charge_window
                if (pulse_first_index > nearby_charge_window and pulse_last_index < (N - nearby_charge_window)) {
                    // we are within 300 bins of the edges of our wf
                    double max_pre_pulse = 1;
                    double max_post_pulse = 1;
                    double baseline_subtracted_off;

                    smoother.Reset(pulse_first_index - nearby_charge_window);
                    for (size_t pre_pulse_it = (pulse_first_index - nearby_charge_window); pre_pulse_it < pulse_first_index; ++pre_pulse_it) {
                        baseline_subtracted_off = smoother.RawValue() + baseline;
                        max_pre_pulse = std::min(max_pre_pulse, baseline_subtracted_off); // baseline is negative and RawValue are positive!
                        smoother.Next();
                    }

                    smoother.Reset(pulse_last_index);
                    for (size_t post_pulse_it = pulse_last_index; post_pulse_it < (pulse_last_index + nearby_charge_window); ++post_pulse_it) {
                        baseline_subtracted_off = smoother.RawValue() + baseline;
                        max_post_pulse = std::min(max_post_pulse, baseline_subtracted_off);
                        smoother.Next();
                    }
                    if (max_pre_pulse < 10 and max_pre_pulse > -10 and max_post_pulse < 10 and max_post_pulse > -10) {
                        // ok! no more than 7 counts above the baseline within our nearby_charge_window
                        // now let's save this pulse!
                        //std::cout << "saving pulse from " << pulse_first_index << " until " << pulse_last_index << std::endl;
                        pulses.push_back({pulse_first_index , pulse_last_index});
                    }

                    smoother.Reset(pulse_last_index + nearby_charge_window);
                    i = pulse_last_index;
                } else {
                    smoother.Reset(pulse_last_index);
                    i = pulse_last_index;
                }
            }
        }
        smoother.Next();
    }
    return pulses;
}

} // namespace

class PulseCollector : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;
    std::string nim_name_;

    // Internal state
    bool geo_seen;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    std::map<CCMPMTKey, std::deque<WaveformSmoother>> smooth_waveform_cache;

    size_t samples_for_baseline;

    double smoother_tau;
    double smoother_delta_t;

    double pulse_initial_deriv_threshold;
    size_t pulse_samples_before;
    size_t pulse_max_start_sample;
    size_t pulse_min_width;
    size_t pulse_max_width;
    double pulse_min_height;
    double pulse_min_deriv_magnitude;
    double pulse_min_integral;

    size_t samples_after_pulse;

    double pulse_activity_threshold;
    size_t pulse_max_rising_edges;
    size_t pulse_max_falling_edges;

    size_t max_samples = 100;

    std::mt19937 rand;

public:
    PulseCollector(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();
    void Flush();

    std::vector<std::tuple<size_t, size_t, double>> ProcessWaveform(CCMPMTKey key, CCMWaveformUInt16 const & waveform);

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(PulseCollector);

std::vector<std::tuple<size_t, size_t, double>> PulseCollector::ProcessWaveform(CCMPMTKey key, CCMWaveformUInt16 const & waveform) {
    // vector we will be returning
    std::vector<std::tuple<size_t, size_t, double>> pulse_results;
    size_t nearby_charge_window = 300;
    // Initialize the waveform smoother
    smooth_waveform_cache[key].emplace_back(waveform.GetWaveform().cbegin(), waveform.GetWaveform().cend(), smoother_delta_t, smoother_tau);
    WaveformSmoother & smoother = smooth_waveform_cache[key].back();

    // Get an initial baseline estimate
    double baseline = EstimateBaseline(smoother, samples_for_baseline);

    // Compute the positions of the pulses
    std::vector<std::pair<size_t, size_t>> pulse_positions =
        FindPulses(smoother, baseline, pulse_max_start_sample, pulse_initial_deriv_threshold, pulse_max_width, pulse_min_width, pulse_min_height, pulse_min_deriv_magnitude, pulse_min_integral);

    if (pulse_positions.size() == 0){
        // oops! no pulses
        pulse_results.emplace_back(0, 0, baseline);
    }

    for (size_t pulse_it = 0; pulse_it < pulse_positions.size(); ++pulse_it){
        //bool found_pulse = (pulse_positions[pulse_it].second > pulse_positions[pulse_it].first) and
        //    (pulse_positions[pulse_it].first < (pulse_samples_before - 1));
        //if(not found_pulse)
        //    continue;

        std::vector<double> baseline_samples;
        baseline_samples.reserve(pulse_positions[pulse_it].second);
        smoother.Reset();
        for(size_t i=0; i<pulse_positions[pulse_it].second; ++i) {
            baseline_samples.push_back(smoother.Value());
            smoother.Next();
        }
        std::sort(baseline_samples.begin(), baseline_samples.end());
        baseline = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());
        size_t pulse_start = pulse_positions[pulse_it].first - nearby_charge_window;
        size_t pulse_end = pulse_positions[pulse_it].second + nearby_charge_window;
        pulse_results.emplace_back(pulse_start, pulse_end, baseline);
    }
    //double activity_threshold = 0;
    //smoother.Reset(pulse_positions.first);
    //for(size_t i=0; i<(pulse_positions.second - pulse_positions.first); ++i) {
    //    double value = smoother.Value() - baseline;
    //    activity_threshold = std::max(activity_threshold, value);
    //    smoother.Next();
    //}
    //activity_threshold = std::min(0.75 * activity_threshold, pulse_activity_threshold);

    //smoother.Reset(pulse_positions.second);
    //double prev_value = smoother.Value() - baseline;
    //double current_value;
    //size_t n_rising = 0;
    //size_t n_falling = 0;
    //for(size_t i=0; i<samples_after_pulse; ++i) {
    //    smoother.Next();
    //    current_value = smoother.Value() - baseline;

    //    if(prev_value < activity_threshold
    //            and current_value >= activity_threshold) {
    //        n_rising += 1;
    //        if(n_rising >= pulse_max_rising_edges)
    //            return {0, baseline};
    //    } else if(current_value < activity_threshold
    //            and prev_value >= activity_threshold) {
    //        n_falling += 1;
    //        if(n_falling >= pulse_max_falling_edges)
    //            return {0, baseline};
    //    }
    //    prev_value = current_value;
    //}
    return pulse_results;
}

void PulseCollector::ProcessFrame(I3FramePtr frame) {
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series const>(waveforms_name_);
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<std::vector<double>>>> pulse_samples = boost::make_shared<I3Map<CCMPMTKey, std::vector<std::vector<double>>>>();
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
        CCMWaveformUInt16 const & waveform = waveforms[channel];
        if(waveform.GetWaveform().size() > 0) {
            std::vector<std::tuple<size_t, size_t, double>> wf_result = ProcessWaveform(key, waveform);
            // std::cout << "length of pulses = " << wf_result.size() << std::endl;
            // now looping over this vector of pulse information
            std::vector<std::vector<double>> vector_of_results;
            WaveformSmoother const & smoother = smooth_waveform_cache[key].back();
            std::pair<std::vector<uint16_t>::const_iterator, std::vector<uint16_t>::const_iterator> iterators = smoother.GetRawWaveform();
            for (size_t wf_it = 0; wf_it < wf_result.size(); ++wf_it){
                size_t first_pulse_sample = std::get<0>(wf_result[wf_it]);
                size_t last_pulse_sample = std::get<1>(wf_result[wf_it]);
                double baseline = std::get<2>(wf_result[wf_it]);
                if(last_pulse_sample > 0) {
                    // Loop over the pulses in the waveform
                    std::vector<double> result;
                    result.reserve(last_pulse_sample - first_pulse_sample + 1);
                    std::vector<uint16_t>::const_iterator wf_it = iterators.first + first_pulse_sample;
                    size_t pulse_idx = first_pulse_sample;
                    while(wf_it != iterators.second and pulse_idx <= last_pulse_sample) {
                        double value = (-double(*wf_it)) - baseline;
                        result.push_back(value);
                        ++wf_it;
                        ++pulse_idx;
                    }
                    vector_of_results.push_back(result);
                }
            }
            pulse_samples->insert(std::make_pair(key, vector_of_results));

        }
        smooth_waveform_cache[key].clear();
    }
    frame->Put("PulseSamples", pulse_samples);
}

PulseCollector::PulseCollector(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_(""), nim_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("NumSamplesForBaseline", "Number of samples to use for initial baseline estimate", size_t(400));
    AddParameter("SmoothingTau", "Time constant for waveform smoothing in ns", double(10.0));
    AddParameter("SmoothingDeltaT", "Bin width in ns for smoothing", double(2.0));
    AddParameter("InitialDerivativeThreshold", "Initial positive derivative threshold for a pulse", double(0.3));
    AddParameter("NumSamplesBeforePulse", "Number of samples required before a pulse", size_t(300));
    AddParameter("MaxPulseStartSample", "Maximum sample from which to start a pulse search", size_t(4400));
    AddParameter("MinPulseWidth", "Minimum width for defining a pulse", size_t(5));
    AddParameter("MaxPulseWidth", "Maxiumum width for defining a pulse", size_t(100));
    AddParameter("MinPulseHeight", "Minimum height for defining a pulse", double(5.0));
    AddParameter("MinPulseDerivativeMagnitude", "Minimum derivative magnitude for defining a pulse", double(0.65));
    AddParameter("MinPulseIntegral", "Minimum integral for defining a pulse", double(25.0));
    AddParameter("NumSamplesAfterPulse", "Number of samples to save after pulse", size_t(300));
    AddParameter("PulseActivityThreshold", "Threshold in ADC counts to indicate downstream", double(20.0));
    AddParameter("MaxRisingEdges", "Maximum number of allowed pulse rising edges within window", size_t(1));
    AddParameter("MaxFallingEdges", "Maximum number of allowed pulse falling edges within window", size_t(1));
}

void PulseCollector::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("NumSamplesForBaseline", samples_for_baseline);
    GetParameter("SmoothingTau", smoother_tau);
    GetParameter("SmoothingDeltaT", smoother_delta_t);
    GetParameter("InitialDerivativeThreshold", pulse_initial_deriv_threshold);
    GetParameter("NumSamplesBeforePulse", pulse_samples_before);
    GetParameter("MaxPulseStartSample", pulse_max_start_sample);
    GetParameter("MinPulseWidth", pulse_min_width);
    GetParameter("MaxPulseWidth", pulse_max_width);
    GetParameter("MinPulseHeight", pulse_min_height);
    GetParameter("MinPulseDerivativeMagnitude", pulse_min_deriv_magnitude);
    GetParameter("MinPulseIntegral", pulse_min_integral);
    GetParameter("NumSamplesAfterPulse", samples_after_pulse);
    GetParameter("PulseActivityThreshold", pulse_activity_threshold);
    GetParameter("MaxRisingEdges", pulse_max_rising_edges);
    GetParameter("MaxFallingEdges", pulse_max_falling_edges);
}

void PulseCollector::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);


    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;

    geo_seen = true;
    PushFrame(frame);
}

void PulseCollector::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    ProcessFrame(frame);
    PushFrame(frame);
}

void PulseCollector::Flush() {
}

void PulseCollector::Finish() {
    Flush();
}
