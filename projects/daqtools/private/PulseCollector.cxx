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

    int CurrentIndex() {
        return index;
    }

    uint16_t RawValue() {
        return begin[index];
    }

    double Value() {
        return smoothed_wf[index];
    }

    double Derivative() {
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

std::pair<size_t, size_t> GetSamplesBeforeThreshold(WaveformSmoother & smoother, double deriv_threshold, size_t max_samples) {
    size_t N = smoother.Size();
    N = std::min(max_samples + 5, N);
    size_t last_index = 0;

    std::vector<double> baseline_samples;
    baseline_samples.reserve(400);
    smoother.Reset();
    for(size_t i=0; i<N; ++i) {
        baseline_samples.push_back(smoother.Value());
        smoother.Next();
    }
    std::sort(baseline_samples.begin(), baseline_samples.end());
    double baseline = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());

    smoother.Reset();
    for(size_t i=0; i<N; ++i) {
        if(smoother.Derivative() > deriv_threshold) {
            size_t pulse_last_idx = CheckForPulse(smoother, i, 100, 5, 5.0, 0.65, 25.0, baseline);
            if(pulse_last_idx > i) {
                last_index = size_t(std::max(ptrdiff_t(0), ptrdiff_t(i) - 5));
                break;
            }
        }
        smoother.Next();
    }
    return {0, last_index};
}

class OnlineStats {
    double w_sum = 0;
    double w_sum2 = 0;
    double mean = 0;
    double S = 0;
public:
    void AddWeightedValue(double x, double w) {
        w_sum = w_sum + w;
        w_sum2 = w_sum2 + std::copysign(w * w, w);
        double mean_old = mean;
        mean = mean_old + (w / w_sum) * (x - mean_old);
        S = S + w * (x - mean_old) * (x - mean);
    }
    void AddValue(double x) {
        AddWeightedValue(x, 1.0);
    }
    void RemoveValue(double x) {
        AddWeightedValue(x, -1.0);
    }
    double PopulationVariance() {
        return S / w_sum;
    }
    double SampleFrequencyVariance() {
        return S / (w_sum - 1);
    }
    double SampleReliabilityVariance() {
        return S / (w_sum - w_sum2 / w_sum);
    }
    double Mean() {
        return mean;
    }
};

} // namespace


class BaselineEstimator {
    double deriv_threshold;
    double tau = 10.0;
    double delta_t = 2.0;
    size_t max_samples = 100;
    size_t min_samples = 10;
    size_t target_samples = 50;

    std::deque<std::vector<double>> ordered_baseline_samples;
    std::multiset<double> sorted_baseline_samples;
public:
    BaselineEstimator(double deriv_threshold, double tau, double delta_t, size_t max_samples, size_t min_samples, size_t target_samples) :
        deriv_threshold(deriv_threshold),
        tau(tau),
        delta_t(delta_t),
        max_samples(max_samples),
        min_samples(min_samples),
        target_samples(target_samples) {
    }

    void AddBaselineSamples(std::vector<double> const & samples) {
        ordered_baseline_samples.push_back(samples);
        sorted_baseline_samples.insert(samples.begin(), samples.end());
    }

    void RemoveBaselineSamples() {
        if(ordered_baseline_samples.size() == 0)
            return;
        {
            std::vector<double> const & to_remove = ordered_baseline_samples.front();
            for(double const & d : to_remove)
                sorted_baseline_samples.erase(d);
        }
        ordered_baseline_samples.pop_front();
    }

    double EstimateBaseline() {
        return robust_stats::Mode(sorted_baseline_samples.begin(), sorted_baseline_samples.end());
    }

    double EstimateBaselineStddev(double median) {
        return robust_stats::MedianAbsoluteDeviation(
                sorted_baseline_samples.begin(),
                sorted_baseline_samples.end(),
                median);
    }
};


class PulseCollector : public I3Module {
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
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    std::map<CCMPMTKey, std::deque<WaveformSmoother>> smooth_waveform_cache;

    double deriv_threshold = 0.3;
    double tau = 10.0;
    double delta_t = 2.0;
    double max_samples = 100;
    size_t min_samples = 10;
    size_t target_samples = 50;
    std::mt19937 rand;

public:
    PulseCollector(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();
    void Flush();

    size_t ProcessWaveform(CCMPMTKey key, size_t frame_number, CCMWaveformUInt16 const & waveform);

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(PulseCollector);

size_t PulseCollector::ProcessWaveform(CCMPMTKey key, size_t frame_number, CCMWaveformUInt16 const & waveform) {
    smooth_waveform_cache[key].emplace_back(waveform.GetWaveform().cbegin(), waveform.GetWaveform().cend(), delta_t, tau);
    WaveformSmoother & smoother = smooth_waveform_cache[key].back();
    std::pair<size_t, size_t> samples_before_threshold =
        GetSamplesBeforeThreshold(smoother, deriv_threshold, max_samples);
    smoother.Reset(samples_before_threshold.second);
    if(samples_before_threshold.second > samples_before_threshold.first) {
        for(size_t i=0; i<3000; ++i) {
            smoother.Next();
        }
    }
    return samples_before_threshold.second;
}

void PulseCollector::ProcessFrame(I3FramePtr frame) {
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series const>(waveforms_name_);
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> pulse_samples = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
        CCMWaveformUInt16 const & waveform = waveforms[channel];
        if(waveform.GetWaveform().size() > 0) {
            size_t last_pulse_sample = ProcessWaveform(key, 0, waveform);
            if(last_pulse_sample > 0) {
                WaveformSmoother const & smoother = smooth_waveform_cache[key].back();
                std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoothed_iterators = smoother.GetSmoothedWaveform();
                pulse_samples->emplace(std::make_pair(key, std::vector<double>(smoothed_iterators.first, smoothed_iterators.second)));
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
}

void PulseCollector::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
}

void PulseCollector::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);


    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;

    for(std::pair<CCMPMTKey const, uint32_t> p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        smooth_waveform_cache.insert({key, std::deque<WaveformSmoother>()});
    //    baseline_estimators.insert({key, std::deque<std::pair<size_t, BaselineEstimator>>()});
    }

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
    //// If we are currently accumulating frames for a new baseline estimate
    //// then we need to finish the estimate and clear out all the frames before
    //// the tray can finish processing
    //if(cached_frames.size()) {
    //    EstimateBaselines();
    //    for(size_t i=0; i<cached_frames.size(); ++i) {
    //        ProcessFrame(cached_frames[i]);
    //        PushFrame(cached_frames[i]);
    //    }
    //    cached_frames.clear();
    //}
    Flush();
}
