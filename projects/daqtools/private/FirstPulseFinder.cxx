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

size_t CheckForPulse(WaveformSmoother smoother, size_t start_idx, size_t max_samples, size_t min_length, double value_threshold, double derivative_threshold, double integral_threshold, double baseline) {
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

double EstimateBaseline(WaveformSmoother smoother, size_t samples_for_baseline) {
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

class FirstPulseFinder : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;
    std::string nim_name_;

    // Internal state
    bool geo_seen;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;

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

    std::mt19937 rand;

public:
    FirstPulseFinder(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();
    void Flush();

    void ProcessFrame(I3FramePtr frame);
};

I3_MODULE(FirstPulseFinder);

void FirstPulseFinder::ProcessFrame(I3FramePtr frame) {
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series const>(waveforms_name_);
    boost::shared_ptr<I3Map<CCMPMTKey, std::pair<size_t, size_t>>> pulse_positions = boost::make_shared<I3Map<CCMPMTKey, std::pair<size_t, size_t>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, WaveformSmoother>> smooth_waveforms = boost::make_shared<I3Map<CCMPMTKey, WaveformSmoother>>();
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey const & key = p.first;
        uint32_t const & channel = p.second;
        CCMWaveformUInt16 const & waveform = waveforms[channel];
        if(waveform.GetWaveform().size() > 0) {
            smooth_waveforms->insert({key, WaveformSmoother(waveform.GetWaveform().cbegin(), waveform.GetWaveform().cend(), smoother_delta_t, smoother_tau)});
            WaveformSmoother & smoother = smooth_waveforms->at(key);
            double baseline = EstimateBaseline(smoother, samples_for_baseline);
            pulse_positions->operator[](key) =
                FindFirstPulse(smoother,
                        baseline,
                        pulse_max_start_sample,
                        pulse_initial_deriv_threshold,
                        pulse_max_width,
                        pulse_min_width,
                        pulse_min_height,
                        pulse_min_deriv_magnitude,
                        pulse_min_integral);
        }
    }
    frame->Put("PulsePositions", pulse_positions);
    frame->Put("WaveformSmoothers", smooth_waveforms);
}

FirstPulseFinder::FirstPulseFinder(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_(""), nim_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("NumSamplesForBaseline", "Number of samples to use for initial baseline estimate", size_t(400));
    AddParameter("SmoothingTau", "Time constant for waveform smoothing in ns", double(10.0));
    AddParameter("SmoothingDeltaT", "Bin width in ns for smoothing", double(2.0));
    AddParameter("InitialDerivativeThreshold", "Initial positive derivative threshold for a pulse", double(0.3));
    AddParameter("NumSamplesBeforePulse", "Number of samples required before a pulse", size_t(400));
    AddParameter("MaxPulseStartSample", "Maximum sample from which to start a pulse search", size_t(4400));
    AddParameter("MinPulseWidth", "Minimum width for defining a pulse", size_t(5));
    AddParameter("MaxPulseWidth", "Maxiumum width for defining a pulse", size_t(100));
    AddParameter("MinPulseHeight", "Minimum height for defining a pulse", double(5.0));
    AddParameter("MinPulseDerivativeMagnitude", "Minimum derivative magnitude for defining a pulse", double(0.65));
    AddParameter("MinPulseIntegral", "Minimum integral for defining a pulse", double(25.0)); // Usual integral is between 30 and 40 with an error 10, probably good to go lower than 25, possibly 10
}

void FirstPulseFinder::Configure() {
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
}

void FirstPulseFinder::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);


    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;

    geo_seen = true;
    PushFrame(frame);
}

void FirstPulseFinder::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    ProcessFrame(frame);
    PushFrame(frame);
}

void FirstPulseFinder::Flush() {
}

void FirstPulseFinder::Finish() {
    Flush();
}
