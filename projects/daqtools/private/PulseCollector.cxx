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
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>

namespace {

class OnlineRobustStats {
    std::deque<double> buffer;
    std::multiset<double> sorted_samples;
public:
    OnlineRobustStats() {}

    void AddValue(double x) {
        buffer.push_back(x);
        sorted_samples.insert(x);
    }

    void RemoveValue() {
        if(buffer.size() == 0)
            return;
        {
            double const & to_remove = buffer.front();
            sorted_samples.erase(to_remove);
        }
        buffer.pop_front();
    }

    double Median() {
        const size_t N = sorted_samples.size();
        const size_t half = N / 2;
        std::multiset<double>::iterator it = sorted_samples.begin();
        for(size_t i=0; i<(half-1); ++i)
            ++it;

        // Odd count: return middle
        if(N % 2) {
            ++it;
            return *it;
        }

        // Even count: return average of middle two.
        double ret = *it;
        ++it;
        ret += *it;
        ret /= 2;
        return ret;
    }

    double Mode() {
        return robust_stats::Mode(sorted_samples.begin(), sorted_samples.end());
    }

    double Stddev(double median) {
        return robust_stats::MedianAbsoluteDeviation(
                sorted_samples.begin(),
                sorted_samples.end(),
                median);
    }
};

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

class WaveformSmoother {
    std::vector<uint16_t>::const_iterator begin;
    std::vector<uint16_t>::const_iterator end;
    size_t N;
    std::vector<double> smoothed_wf;
    std::vector<double> derivative;
    size_t index;
    size_t max_computed;
    double delta_t;
    double tau;
    double current_value;
    double current_derivative;
    double alpha;
    double y_i;
    double exp_prevprev, exp_prev, exp_current, exp_next, exp_nextnext;
    double x;
public:
    WaveformSmoother(std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end, double delta_t, double tau) :
        begin(begin), end(end), N(std::distance(begin, end)), smoothed_wf(N), derivative(N), index(0), max_computed(0), delta_t(delta_t), tau(tau) {
        int init_avg_end_idx = std::min(3 * std::max(int(tau / delta_t), 1), int(N));
        double exp_start = 0.0;
        std::vector<uint16_t>::const_iterator x_it = begin;
        for(size_t i=0; i < init_avg_end_idx and x_it != end; ++i, ++x_it) {
            exp_start += -double(*x_it);
        }
        exp_start /= init_avg_end_idx;
        alpha = exp(-delta_t / tau);
        y_i = exp_start;

        // Exponentially smooth the first element
        x = -double(begin[0]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_prevprev = y_i;

        // Exponentially smooth the second element
        x = -double(begin[1]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_prev = y_i;

        // Exponentially smooth the third element
        x = -double(begin[2]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_current = y_i;

        // Exponentially smooth the fourth element
        x = -double(begin[3]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_next = y_i;

        // Box smoothing for the first element
        smoothed_wf[0] = (exp_prevprev + exp_prev) / 2.0;

        // Box smoothing for the second element
        smoothed_wf[1] = (exp_prevprev + exp_prev + exp_current) / 3.0;

        // Box smoothing for the third element
        smoothed_wf[2] = (exp_prev + exp_current + exp_next) / 3.0;

        // Derivative for the first element
        derivative[0] = (-3.0 * smoothed_wf[0] + 4.0 * smoothed_wf[1] - smoothed_wf[2]) / (2.0 * delta_t);

        index = 0;

        current_value = smoothed_wf[index];
        current_derivative = derivative[index];
    }

    void Reset() {
        index = 0;
    }

    void Reset(size_t reset_index) {
        index = std::min(index, reset_index);
    }

    void Reset(std::vector<double>::const_iterator end) {
        index = std::min(ptrdiff_t(index), std::distance(smoothed_wf.cbegin(), end) - 1);
    }

    std::pair<std::vector<uint16_t>::const_iterator, std::vector<uint16_t>::const_iterator> GetRawWaveform() const {
        return {begin, begin + index+1};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetSmoothedWaveform() const {
        return {smoothed_wf.cbegin(), smoothed_wf.cbegin() + index+1};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetFullSmoothedWaveform() const {
        return {smoothed_wf.cbegin(), smoothed_wf.cend()};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetDerivative() const {
        return {derivative.cbegin(), derivative.cbegin() + index+1};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetFullDerivative() const {
        return {derivative.cbegin(), derivative.cend()};
    }

    int CurrentIndex() const {
        return index;
    }

    uint16_t RawValue() const {
        return begin[index];
    }

    double Value() const {
        return smoothed_wf[index];
    }

    double Derivative() const {
        return derivative[index];
    }

    void Next() {
        if(index < N-1)
            index += 1;
        if(index <= max_computed)
            return;
        size_t idx = index;
        // Exponential smoothing for the nextnext element
        // Box smoothing for the next element
        // Derivative for the current element
        if(idx < N-2) {
            x = -double(begin[idx+2]);
            y_i += (1.0 - alpha) * (x - y_i);
            exp_nextnext = y_i;

            // Box smoothing for everything in between
            smoothed_wf[idx+1] = (exp_current + exp_next + exp_nextnext) / 3.0;
            exp_current = exp_next;
            exp_next = exp_nextnext;

            // Finite difference derivative
            derivative[idx] = (smoothed_wf[idx+1] - smoothed_wf[idx-1]) / (2.0 * delta_t);
        } else if(idx == N-2) {
            // Box smoothing for the last element
            smoothed_wf[N-1] = (exp_current + exp_next) / 2.0;

            // Derivative for the second to last element
            derivative[N-2] = (smoothed_wf[N-1] - smoothed_wf[N-3]) / (2.0 * delta_t);

            // Derivative for the last element
            derivative[N-1] = (3.0 * smoothed_wf[N-1] - 4.0 * smoothed_wf[N-2] + smoothed_wf[N-3]) / (2.0 * delta_t);
        }
        current_value = smoothed_wf[idx];
        current_derivative = derivative[idx];
        max_computed = idx;
    }

    size_t Size() {
        return N;
    }
};

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

std::pair<size_t, size_t> FindFirstPulse(WaveformSmoother & smoother, double baseline, size_t pulse_max_start_sample, double deriv_threshold, size_t max_pulse_width, size_t min_pulse_width, double min_pulse_height, double min_deriv_magnitude, double min_integral) {
    size_t N = smoother.Size();
    N = std::min(N, pulse_max_start_sample);
    if(N > 5)
        N -= 5;
    size_t pulse_first_index = 0;
    size_t pulse_last_index = 0;

    smoother.Reset();
    for(size_t i=0; i<N; ++i) {
        if(smoother.Derivative() > deriv_threshold) {
            pulse_last_index = CheckForPulse(smoother, i, max_pulse_width, min_pulse_width, min_pulse_height, min_deriv_magnitude, min_integral, baseline);
            if(pulse_last_index > i) {
                pulse_first_index = size_t(std::max(ptrdiff_t(0), ptrdiff_t(i) - 5));
                break;
            }
        }
        smoother.Next();
    }
    return {pulse_first_index, pulse_last_index};
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

    std::pair<size_t, double> ProcessWaveform(CCMPMTKey key, CCMWaveformUInt16 const & waveform);

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(PulseCollector);

std::pair<size_t, double> PulseCollector::ProcessWaveform(CCMPMTKey key, CCMWaveformUInt16 const & waveform) {
    smooth_waveform_cache[key].emplace_back(waveform.GetWaveform().cbegin(), waveform.GetWaveform().cend(), smoother_delta_t, smoother_tau);
    WaveformSmoother & smoother = smooth_waveform_cache[key].back();
    double baseline = EstimateBaseline(smoother, samples_for_baseline);
    std::pair<size_t, size_t> pulse_positions =
        FindFirstPulse(smoother, baseline, pulse_max_start_sample, pulse_initial_deriv_threshold, pulse_max_width, pulse_min_width, pulse_min_height, pulse_min_deriv_magnitude, pulse_min_integral);
    bool found_pulse = pulse_positions.second > pulse_positions.first;
    if(not found_pulse)
        return {0, baseline};
    if(pulse_positions.first < pulse_samples_before - 1)
        return {0, baseline};

    std::vector<double> baseline_samples;
    baseline_samples.reserve(pulse_positions.second);
    smoother.Reset();
    for(size_t i=0; i<pulse_positions.second; ++i) {
        baseline_samples.push_back(smoother.Value());
        smoother.Next();
    }
    std::sort(baseline_samples.begin(), baseline_samples.end());
    baseline = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());

    double activity_threshold = 0;
    smoother.Reset(pulse_positions.first);
    for(size_t i=0; i<(pulse_positions.second - pulse_positions.first); ++i) {
        double value = smoother.Value() - baseline;
        activity_threshold = std::max(activity_threshold, value);
        smoother.Next();
    }
    activity_threshold = std::min(0.75 * activity_threshold, pulse_activity_threshold);

    smoother.Reset(pulse_positions.second);
    double prev_value = smoother.Value() - baseline;
    double current_value;
    size_t n_rising = 0;
    size_t n_falling = 0;
    for(size_t i=0; i<samples_after_pulse; ++i) {
        smoother.Next();
        current_value = smoother.Value() - baseline;

        if(prev_value < activity_threshold
                and current_value >= activity_threshold) {
            n_rising += 1;
            if(n_rising >= pulse_max_rising_edges)
                return {0, baseline};
        } else if(current_value < activity_threshold
                and prev_value >= activity_threshold) {
            n_falling += 1;
            if(n_falling >= pulse_max_falling_edges)
                return {0, baseline};
        }
        prev_value = current_value;
    }
    return {pulse_positions.second, baseline};
}

void PulseCollector::ProcessFrame(I3FramePtr frame) {
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series const>(waveforms_name_);
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> pulse_samples = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
        CCMWaveformUInt16 const & waveform = waveforms[channel];
        if(waveform.GetWaveform().size() > 0) {
            std::pair<size_t, double> wf_result = ProcessWaveform(key, waveform);
            size_t last_pulse_sample = wf_result.first;
            double baseline = wf_result.second;
            if(last_pulse_sample > 0) {
                WaveformSmoother const & smoother = smooth_waveform_cache[key].back();
                std::pair<std::vector<uint16_t>::const_iterator, std::vector<uint16_t>::const_iterator> iterators = smoother.GetRawWaveform();
                std::vector<double> result;
                result.reserve(smoother.CurrentIndex());
                std::vector<uint16_t>::const_iterator it = iterators.first;
                while(it != iterators.second) {
                    double value = (-double(*it)) - baseline;
                    result.push_back(value);
                    ++it;
                }
                pulse_samples->insert(std::make_pair(key, result));
            }
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
    AddParameter("NumSamplesBeforePulse", "Number of samples required before a pulse", size_t(200));
    AddParameter("MaxPulseStartSample", "Maximum sample from which to start a pulse search", size_t(4400));
    AddParameter("MinPulseWidth", "Minimum width for defining a pulse", size_t(5));
    AddParameter("MaxPulseWidth", "Maxiumum width for defining a pulse", size_t(100));
    AddParameter("MinPulseHeight", "Minimum height for defining a pulse", double(5.0));
    AddParameter("MinPulseDerivativeMagnitude", "Minimum derivative magnitude for defining a pulse", double(0.65));
    AddParameter("MinPulseIntegral", "Minimum integral for defining a pulse", double(25.0));
    AddParameter("NumSamplesAfterPulse", "Number of samples to save after pulse", size_t(3000));
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
