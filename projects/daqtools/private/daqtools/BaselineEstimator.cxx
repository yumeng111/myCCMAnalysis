#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <nlopt.h>
#include <nlopt.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <thread>
#include <cmath>
#include <functional>

#include <icetray/ctpl.h>
#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WindowedStats.h"
#include <dataclasses/calibration/BaselineEstimate.h>

struct BaselineEstimatorJob {
    std::atomic<bool> running = false;
    std::thread thread;
    size_t thread_index = 0;
    I3FramePtr frame = nullptr;
    size_t frame_index = 0;
};

struct BaselineEstimatorResult {
    I3FramePtr frame = nullptr;
    bool done = false;
};

void LinearFlatFit(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end, double & linear_b, double & linear_m, double & linear_sigma, double & flat_mean, double & flat_sigma) {
    double N = std::distance(begin, end);
    double X = 0.0;
    double X2 = 0.0;
    double T = 0.0;
    double T2 = 0.0;
    double XT = 0.0;

    size_t t = 0;
    std::vector<double>::const_iterator it = begin;
    for(; it != end; ++it) {
        double x = *it;
        X += x;
        X2 += x*x;
        T += t;
        T2 += t*t;
        XT += t*x;

        ++t;
    }
    double b = X/N - (N*XT*T-X*T*T) / (N*N*T2 - N*T*T);
    double m = (N*XT - X*T) / (N*T2 - T*T);
    double sigma2 = (X2 + N*b*b + m*m*T2 - 2*b*X - 2*m*XT) / N;
    linear_b = b;
    linear_m = m;
    linear_sigma = sqrt(sigma2);
    flat_mean = X / N;

    flat_sigma = 0;
    it = begin;
    for(; it != end; ++it) {
        flat_sigma += ((*it) - flat_mean) * ((*it) - flat_mean);
    }
    flat_sigma /= N;
}

double FlatChi2PerDOF(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end, double flat_mean, double flat_sigma) {
    size_t N = std::distance(begin, end);
    size_t DOF = std::max(N, size_t(3)) - 2;
    double stddev = flat_sigma;
    double pred = flat_mean;
    double chi2 = 0.0;
    for(; begin != end; ++begin) {
        double z = ((*begin) - pred) / stddev;
        chi2 += z * z;
    }
    chi2 /= 2.0;
    chi2 += log(stddev) + log(sqrt(2.0 * M_PI));
    return chi2 /= DOF;
}

double LinearChi2PerDOF(
        std::vector<double>::const_iterator begin,
        std::vector<double>::const_iterator end,
        double linear_intercept,
        double linear_slope,
        double linear_sigma) {
    size_t N = std::distance(begin, end);
    size_t DOF = std::max(N, size_t(4)) - 3;
    double stddev = linear_sigma;
    double chi2 = 0.0;
    size_t t = 0;
    for(; begin != end; ++begin, ++t) {
        double pred = linear_intercept + linear_slope * t;
        double z = ((*begin) - pred) / stddev;
        chi2 += z * z;
    }
    chi2 /= 2.0;
    chi2 += log(stddev) + log(sqrt(2.0 * M_PI));
    return chi2 /= DOF;
}

double ExponentialChi2PerDOF(
        std::vector<double>::const_iterator begin,
        std::vector<double>::const_iterator end,
        double a,
        double b,
        double c) {
    size_t N = std::distance(begin, end);
    size_t DOF = std::max(N, size_t(4)) - 3;
    double stddev = b;
    double chi2 = 0.0;
    size_t t = 0;
    for(; begin != end; ++begin, t += 2) {
        double pred = a + b * std::exp(c * t);
        double z = ((*begin) - pred) / stddev;
        chi2 += z * z;
    }
    chi2 /= 2.0;
    chi2 += log(stddev) + log(sqrt(2.0 * M_PI));
    return chi2 /= DOF;
}

class BaselineEstimator: public I3Module {
    std::exception_ptr teptr = nullptr;
    std::string geometry_key_;
    std::string daq_config_key_;
    std::string ccm_waveforms_key_;
    std::string output_key_;

    bool geo_seen_;
    CCMGeometry geo_;

    I3Vector<CCMPMTKey> allowed_pmt_keys_;
    I3Vector<uint32_t> allowed_channels_;
    std::vector<CCMPMTKey> pmt_keys_;
    std::vector<uint32_t> channels_;

    bool output_channels_;

    size_t num_threads;
    size_t max_cached_frames;

    size_t num_samples;
    size_t num_exp_samples;

    size_t frame_index = 0;
    size_t min_frame_idx = 0;
    std::deque<BaselineEstimatorJob *> free_jobs;
    std::deque<BaselineEstimatorJob *> running_jobs;
    std::deque<BaselineEstimatorResult> results;
public:
    BaselineEstimator(const I3Context&);
    void Configure();
    void Finish();
    void DAQ(I3FramePtr frame);
    void Geometry(I3FramePtr frame);
};

template <
    class Iterator,
    class U = typename std::iterator_traits<Iterator>::value_type>
double Average(Iterator begin, Iterator end) {
    double total = 0;
    double count = std::distance(begin, end);;
    for(; begin != end; ++begin) {
        total += *begin;
    }
    return total / count;
}

template<typename InputIterator, typename OutputIterator>
void OutlierFilter(double starting_estimate, InputIterator begin, InputIterator end, OutputIterator result_begin) {
    double delta_tau = 20;
    double prev_tau = 2.0;
    double next_tau = 2.0;
    double value = starting_estimate;

    // now let's loop over the waveform
    size_t N = std::distance(begin, end);
    for(size_t wf_it = 0; wf_it < N; ++wf_it) {
        double delta = (*begin) - value;
        double e = std::fabs(delta) / delta_tau;

        if(wf_it > 0) {
            e += std::fabs((*begin) - (*(begin-1))) / prev_tau;
        }

        if(wf_it + 1 < N) {
            e += std::fabs((*begin) - (*(begin+1))) / next_tau;
        }
        delta *= std::exp(-e);
        value += delta;
        (*result_begin) = value;
        ++result_begin;
        ++begin;
    }
}

void FitExponential(std::vector<double> const & y, double & a, double & b, double & c) {
    std::vector<double> x (y.size());
    // let's fill our x vals (aka time)
    for (size_t time_it = 0 ; time_it < y.size(); ++time_it){
        x[time_it] = time_it * 2;
    }

    // now fitting
    std::vector<double> S(y.size());
    double S2_sum = 0;
    S[0] = 0;
    double sum_S = 0;

    for (size_t result_it = 1; result_it < S.size(); ++result_it){
        sum_S += 0.5 * (y[result_it] + y[result_it-1]) * (x[result_it] - x[result_it-1]);
        S[result_it] = sum_S;
        S2_sum += sum_S * sum_S;
    }

    double m00 = 0;
    double m01 = 0;
    double m11 = S2_sum;
    double v0 = 0;
    double v1 = 0;

    for (size_t x_it = 0; x_it < x.size(); ++x_it){
        m00 += std::pow(x[x_it] - x[0] , 2);
        m01 += (x[x_it] - x[0]) * S[x_it];
        v0 += (y[x_it] - y[0]) * (x[x_it] - x[0]);
        v1 += (y[x_it] - y[0]) * S[x_it];
    }

    double m10 = m01;

    // now let's find the inverse of m
    double m_inv00= 0;
    double m_inv01 = 0;
    double m_inv10 = 0;
    double m_inv11 = 0;
    double prefactor = 1/(m00 * m11 - m01 * m10);

    m_inv00 = prefactor * m11;
    m_inv01 = prefactor * -1 * m01;
    m_inv10 = prefactor * -1 * m10;
    m_inv11 = prefactor * m00;

    // let's now multiply m_inv.v
    double mult_0 = 0;
    double mult_1 = 0;

    mult_0 = m_inv00 * v0 + m_inv01 * v1;
    mult_1 = m_inv01 * v0 + m_inv11 * v1;

    c = double(mult_1); // one of the parameters we're solving for!!!

    // let's now do it again!
    m00 = 0;
    m01 = 0;
    m10 = 0;
    m11 = 0;
    v0 = 0;
    v1 = 0;

    m00 = x.size();

    for (size_t i = 0; i < x.size(); ++i){
        m01 += std::exp(c * x[i]);
        m11 += std::exp(2.0 * c *  x[i]);
        v0 += y[i];
        v1 += y[i] * std::exp(c * x[i]);
    }
    m10 = m01;

    // now let's invert m again
    prefactor = 1/(m00 * m11 - m01 * m10);

    m_inv00 = prefactor * m11;
    m_inv01 = prefactor * -1 * m01;
    m_inv10 = prefactor * -1 * m10;
    m_inv11 = prefactor * m00;

    a = double(m_inv00 * v0 + m_inv01 * v1);
    b = double(m_inv01 * v0 + m_inv11 * v1);
}

void ProcessWaveform(BaselineEstimate & baseline, std::vector<short unsigned int> const & samples, size_t target_num_samples, size_t target_exp_samples) {
    baseline.target_num_frames = 1;
    if (samples.size() == 0) {
        baseline.baseline = std::numeric_limits<double>::quiet_NaN();
        baseline.stddev = std::numeric_limits<double>::quiet_NaN();
        baseline.num_frames = 1;
        baseline.num_samples = 0;
        return;
    }

    // vector to store results of outlier filter
    size_t N = std::min(samples.size(), target_num_samples);
    std::vector<double> outlier_filter_results;
    outlier_filter_results.reserve(N);
    std::vector<uint16_t> starting_samples(samples.begin(), samples.begin() + std::min(size_t(100), N));
    std::sort(starting_samples.begin(), starting_samples.end());
    double starting_value = robust_stats::Mode(starting_samples.begin(), starting_samples.end());
    OutlierFilter(starting_value, samples.begin(), samples.begin() + N, std::back_inserter(outlier_filter_results));
    starting_value = outlier_filter_results.back();

    double flat_mean = 0;
    double flat_sigma = 0;
    double linear_intercept = 0;
    double linear_slope = 0;
    double linear_sigma = 0;
    LinearFlatFit(outlier_filter_results.begin(), outlier_filter_results.end(), linear_intercept, linear_slope, linear_sigma, flat_mean, flat_sigma);

    double flat_chi2_per_dof_threshold = 1.0;
    double flat_sigma_threshold = 35.0;
    double flat_chi2_per_dof = FlatChi2PerDOF(outlier_filter_results.begin(), outlier_filter_results.end(), flat_mean, flat_sigma);

    bool bad_flat_fit = isnan(flat_mean) or isnan(flat_sigma);
    bool good_flat_fit = flat_chi2_per_dof < flat_chi2_per_dof_threshold and flat_sigma < flat_sigma_threshold;

    if(bad_flat_fit) {
        // Extend the outlier filter to the whole waveform
        OutlierFilter(starting_value, samples.begin() + N, samples.end(), std::back_inserter(outlier_filter_results));
    }

    if(bad_flat_fit or good_flat_fit) {
        // Do nothing
    } else {
        // Flat fit is not good enough, let's try an exponential
        // initializing the exponential fit params
        double a;
        double b;
        double c;

        // Run the outlier filter on the rest of the waveform
        OutlierFilter(starting_value, samples.begin() + N, samples.begin() + std::min(std::max(N, target_exp_samples - target_num_samples), samples.size()), std::back_inserter(outlier_filter_results));
        starting_value = outlier_filter_results.back();
        FitExponential(outlier_filter_results, a, b, c);

        bool bad_exp_fit = isnan(a) or isnan(b) or isnan(c);
        if(bad_exp_fit) {
            OutlierFilter(starting_value, samples.begin() + outlier_filter_results.size(), samples.end(), std::back_inserter(outlier_filter_results));
            // Use the outlier filter result of the whole waveform as is
            // i.e. do nothing
        } else {
            // Subtract off the exponential component from the outlier filter
            for (size_t exp_it = 0; exp_it < outlier_filter_results.size(); ++exp_it){
                outlier_filter_results[exp_it] -= b * std::exp(c * (exp_it * 2.0));
            }
        }
    }

    std::sort(outlier_filter_results.begin(), outlier_filter_results.end());
    double baseline_mode_val = robust_stats::Mode(outlier_filter_results.begin(), outlier_filter_results.end());
    double baseline_std = robust_stats::MedianAbsoluteDeviation(outlier_filter_results.begin(), outlier_filter_results.end(), baseline_mode_val);
    baseline.baseline = -1 * baseline_mode_val;
    baseline.stddev = baseline_std;
    baseline.num_frames = 1;
    baseline.num_samples = outlier_filter_results.size();
    outlier_filter_results.clear();
    outlier_filter_results.shrink_to_fit();
    return;
}

void ProcessWaveformVariance(BaselineEstimate & baseline, std::vector<short unsigned int> const & samples, size_t window_size, double variance_percentile) {
    baseline.target_num_frames = 1;
    if (samples.size() == 0) {
        baseline.baseline = std::numeric_limits<double>::quiet_NaN();
        baseline.stddev = std::numeric_limits<double>::quiet_NaN();
        baseline.num_frames = 1;
        baseline.num_samples = 0;
        return;
    }

    std::vector<double> smoothed_samples(samples.size());
    WaveformSmoother smoother(samples.begin(), samples.end(), 2, 10);
    for(size_t i=0; i<samples.size(); ++i) {
        smoothed_samples[i] = smoother.Value();
        smoother.Next();
    }

    // vector to store results of outlier filter
    size_t N = smoothed_samples.size();
    std::vector<double> outlier_filter_results;
    outlier_filter_results.reserve(N);
    std::vector<uint16_t> starting_samples(smoothed_samples.begin(), smoothed_samples.begin() + std::min(size_t(100), N));
    std::sort(starting_samples.begin(), starting_samples.end());
    double starting_value = robust_stats::Mode(starting_samples.begin(), starting_samples.end());
    OutlierFilter(starting_value, smoothed_samples.begin(), smoothed_samples.begin() + N, std::back_inserter(outlier_filter_results));
    starting_value = outlier_filter_results.back();

    double flat_mean = 0;
    double flat_sigma = 0;
    double linear_intercept = 0;
    double linear_slope = 0;
    double linear_sigma = 0;
    LinearFlatFit(outlier_filter_results.begin(), outlier_filter_results.begin() + 1000, linear_intercept, linear_slope, linear_sigma, flat_mean, flat_sigma);

    double flat_chi2_per_dof_threshold = 1.0;
    double flat_sigma_threshold = 35.0;
    double flat_chi2_per_dof = FlatChi2PerDOF(outlier_filter_results.begin(), outlier_filter_results.begin() + 1000, flat_mean, flat_sigma);

    bool bad_flat_fit = isnan(flat_mean) or isnan(flat_sigma);
    bool good_flat_fit = flat_chi2_per_dof < flat_chi2_per_dof_threshold and flat_sigma < flat_sigma_threshold;

    if(bad_flat_fit or good_flat_fit) {
        // Do nothing
    } else {
        // initializing the exponential fit params
        double a;
        double b;
        double c;

        FitExponential(outlier_filter_results, a, b, c);
        double exp_chi2_per_dof = ExponentialChi2PerDOF(outlier_filter_results.begin(), outlier_filter_results.begin() + std::min(size_t(1000), N), a, b, c);

        bool bad_exp_fit = isnan(a) or isnan(b) or isnan(c);
        if(not bad_exp_fit and c < 0 and b > 0 and c > -1.0/1000.0 and exp_chi2_per_dof < 35.0) {
            // Subtract off the exponential component from the outlier filter
            for (size_t exp_it = 0; exp_it < outlier_filter_results.size(); ++exp_it) {
                outlier_filter_results[exp_it] -= b * std::exp(c * (exp_it * 2.0));
            }
        }
    }

    // vector to store results of outlier filter
    N = std::min(N, window_size);
    if(N == 0) {
        baseline.baseline = std::numeric_limits<double>::quiet_NaN();
        baseline.stddev = 0;
        baseline.num_frames = 1;
        baseline.num_samples = 0;
        return;
    }

    size_t N_for_variance = std::min(smoothed_samples.size(), std::max(size_t(smoothed_samples.size() * variance_percentile), size_t(3)));
    size_t N_half = N_for_variance / 2;
    size_t N_opposite_half = N - std::min(N, N_half);
    std::vector<std::tuple<double, double>> variance_and_values;
    variance_and_values.reserve(N+1);
    WindowedStats stats;
    for(size_t i=0; i<N; ++i) {
        stats.AddValue(smoothed_samples.at(i));
    }
    std::function<bool(std::tuple<double, double>, std::tuple<double, double>)> comp = [](std::tuple<double, double> const & a, std::tuple<double, double> const & b) {
        return std::get<0>(a) < std::get<0>(b);
    };
    for(size_t i=0; i<N_half; ++i) {
        double const var = stats.Variance();
        double const val = smoothed_samples.at(i);

        if(variance_and_values.size() < N_for_variance or var < std::get<0>(variance_and_values.back())) {
            std::vector<std::tuple<double, double>>::iterator it = std::lower_bound(variance_and_values.begin(), variance_and_values.end(), std::make_tuple(var, val), comp);
            variance_and_values.insert(it, std::make_tuple(var, val)); // insert before iterator it
            if(variance_and_values.size() > N_for_variance) {
                variance_and_values.pop_back();
            }
        }
    }
    for(size_t i=N_half; i + N_opposite_half<smoothed_samples.size(); ++i) {
        stats.AddValue(smoothed_samples[i+N_opposite_half]);
        stats.RemoveValue();
        double const var = stats.Variance();
        double const val = smoothed_samples.at(i);
        if(variance_and_values.size() < N_for_variance or var < std::get<0>(variance_and_values.back())) {
            std::vector<std::tuple<double, double>>::iterator it = std::lower_bound(variance_and_values.begin(), variance_and_values.end(), std::make_tuple(var, val), comp);
            variance_and_values.insert(it, std::make_tuple(var, val)); // insert before iterator it
            if(variance_and_values.size() > N_for_variance) {
                variance_and_values.pop_back();
            }
        }
    }
    for(size_t i=smoothed_samples.size()-N_opposite_half; i<smoothed_samples.size(); ++i) {
        double const var = stats.Variance();
        double const val = smoothed_samples.at(i);
        if(variance_and_values.size() < N_for_variance or var < std::get<0>(variance_and_values.back())) {
            std::vector<std::tuple<double, double>>::iterator it = std::lower_bound(variance_and_values.begin(), variance_and_values.end(), std::make_tuple(var, val), comp);
            variance_and_values.insert(it, std::make_tuple(var, val)); // insert before iterator it
            if(variance_and_values.size() > N_for_variance) {
                variance_and_values.pop_back();
            }
        }
    }

    std::vector<double> selected_samples;
    selected_samples.reserve(variance_and_values.size());
    for(std::tuple<double, double> const & var_and_val : variance_and_values) {
        selected_samples.push_back(std::get<1>(var_and_val));
    }

    std::sort(selected_samples.begin(), selected_samples.end());
    double baseline_mode_val = robust_stats::Mode(selected_samples.begin(), selected_samples.end());
    double baseline_std = robust_stats::MedianAbsoluteDeviation(selected_samples.begin(), selected_samples.end(), baseline_mode_val);
    baseline.baseline = baseline_mode_val;
    baseline.stddev = baseline_std;
    baseline.num_frames = 1;
    baseline.num_samples = selected_samples.size();
    return;
}

I3_MODULE(BaselineEstimator);

BaselineEstimator::BaselineEstimator(const I3Context& context) : I3Module(context) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMAnalysis::Binary::CCMDAQConfig>::value()));
    AddParameter("CCMWaveformsName", "Key for CCMWaveforms object", std::string("CCMWaveforms"));
    AddParameter("PMTKeys", "PMTKeys to run over", I3Vector<CCMPMTKey>());
    AddParameter("Channels", "PMTKeys to run over", I3Vector<uint32_t>());
    AddParameter("OutputChannels", "Output based on channel number rather than PMTKey", bool(false));
    AddParameter("NumThreads", "Number of worker threads to use for baseline estimation", (size_t)(0));
    AddParameter("NumSamples", "Number of samples to use for the baseline estimation", (size_t)(800));
    AddParameter("NumExpSamples", "Number of samples to use for the exponential fit", (size_t)(4000));
    AddParameter("OutputName", "Key to save output I3Vector<BaselineEstimate> to", std::string("BaselineEstimates"));
    AddParameter("MaxCachedFrames", "The maximum number of frames this module is allowed to have cached", (size_t)(1000));
}

void BaselineEstimator::Configure() {
    GetParameter("CCMGeometryName", geometry_key_);
    GetParameter("CCMDAQConfigName", daq_config_key_);
    GetParameter("CCMWaveformsName", ccm_waveforms_key_);
    GetParameter("PMTKeys", allowed_pmt_keys_);
    GetParameter("Channels", allowed_channels_);
    GetParameter("OutputChannels", output_channels_);
    GetParameter("NumThreads", num_threads);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
    GetParameter("NumSamples", num_samples);
    GetParameter("NumExpSamples", num_exp_samples);
    GetParameter("OutputName", output_key_);
    GetParameter("MaxCachedFrames", max_cached_frames);
}

void BaselineEstimator::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_key_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_key_.c_str());
    }

    CCMGeometryConstPtr geo = frame->Get<CCMGeometryConstPtr>(geometry_key_);
    if (!geo)
        log_fatal("Couldn't find '%s' in the frame!",
                geometry_key_.c_str());
    geo_ = *geo;
    geo_seen_ = true;

    channels_.clear();
    pmt_keys_.clear();

    CCMAnalysis::Binary::CCMDAQConfigConstPtr daq_config = frame->Get<CCMAnalysis::Binary::CCMDAQConfigConstPtr>(daq_config_key_);
    if (!daq_config)
        log_fatal("Couldn't find '%s' in the frame!",
                daq_config_key_.c_str());

    if(output_channels_) {
        // Copy and sort the channels from the daq config
        std::vector<uint32_t> geo_channels;
        uint32_t channel_index = 0;
        for(CCMAnalysis::Binary::DigitizerBoard const & board : daq_config->digitizer_boards) {
            for(CCMAnalysis::Binary::ChannelHeader const & channel : board.channels) {
                geo_channels.push_back(channel_index);
                channel_index += 1;
            }
        }

        // Assume we're doing all PMTs if none are specified
        if(allowed_channels_.size() == 0) {
            channels_ = geo_channels;
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
                log_warn("%s", ss.str().c_str());
            }
            log_info("Using %zu channels", channels_.size());
        }
    } else {
        // Copy and sort the pmt keys from the geometry
        std::vector<CCMPMTKey> geo_keys; geo_keys.reserve(geo_.pmt_channel_map.size());
        for(std::pair<CCMPMTKey const, uint32_t> const & p : geo_.pmt_channel_map)
            geo_keys.push_back(p.first);
        std::sort(geo_keys.begin(), geo_keys.end());

        // Assume we're doing all PMTs if none are specified
        if(allowed_pmt_keys_.size() == 0) {
            pmt_keys_ = geo_keys;
        } else {
            // Clear the final pmt key list
            pmt_keys_.clear();
            // Fill the final pmt key list with the intersection of specified pmt keys and those available in the geometry
            // If allowed_pmt_keys_ is not sorted then this will fail horribly
            std::set_intersection(geo_keys.begin(), geo_keys.end(), allowed_pmt_keys_.begin(), allowed_pmt_keys_.end(), std::back_inserter(pmt_keys_));

            if(pmt_keys_.size() < allowed_pmt_keys_.size()) {
                std::vector<CCMPMTKey> missing_keys;
                std::set_difference(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end(), pmt_keys_.begin(), pmt_keys_.end(), std::back_inserter(missing_keys));
                std::stringstream ss;
                ss << "Some specified CCMPMTKeys are not present in the geometry:";
                for(CCMPMTKey const & key : missing_keys)
                    ss << " " << key;
                log_warn("%s", ss.str().c_str());
            }
        }
    }

    PushFrame(frame);
}

void FrameThread(std::atomic<bool> & running, I3Frame * frame, std::string const & ccm_waveforms_key_, std::string const & output_key_, CCMGeometry const & geo_, std::vector<CCMPMTKey> const & pmt_keys_, std::vector<uint32_t> const & channels_, bool output_channels_, size_t num_samples, size_t num_exp_samples) {
    // let's read in our waveform
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>(ccm_waveforms_key_);
    if(waveforms == nullptr) {
        log_fatal("No CCMWaveformUInt16Series under key name \"%s\"", ccm_waveforms_key_.c_str());
    }

    size_t num_baselines;
    if(output_channels_) {
        num_baselines = channels_.size();
    } else {
        num_baselines = pmt_keys_.size();
    }
    std::vector<BaselineEstimate> baseline_estimates;
    baseline_estimates.resize(num_baselines);
    std::vector<BaselineEstimate> baseline_estimates_variance;
    baseline_estimates_variance.resize(num_baselines);

    // loop over each channel in waveforms
    for(size_t i=0; i<num_baselines; ++i) {
        // let's get our baseline
        uint32_t channel;
        if(output_channels_) {
            channel = channels_[i];
        } else {
            channel = geo_.pmt_channel_map.at(pmt_keys_[i]);
        }
        ProcessWaveform(
            baseline_estimates[i],
            std::cref(waveforms->at(channel).GetWaveform()),
            num_samples,
            num_exp_samples
        );
        ProcessWaveformVariance(
            baseline_estimates_variance[i],
            std::cref(waveforms->at(channel).GetWaveform()),
            5,
            0.1
        );
    }

    if(output_channels_) {
        boost::shared_ptr<I3Map<uint32_t, BaselineEstimate>> baselines = boost::make_shared<I3Map<uint32_t, BaselineEstimate>>();
        for(size_t i = 0; i < num_baselines; ++i) {
            baselines->emplace(channels_[i], baseline_estimates[i]);
        }
        frame->Put(output_key_, baselines);
        boost::shared_ptr<I3Map<uint32_t, BaselineEstimate>> baselines_variance = boost::make_shared<I3Map<uint32_t, BaselineEstimate>>();
        for(size_t i = 0; i < num_baselines; ++i) {
            baselines_variance->emplace(channels_[i], baseline_estimates_variance[i]);
        }
        frame->Put(output_key_ + "Variance", baselines_variance);
    } else {
        // I3Map to store pmt key and baselines
        boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> baselines = boost::make_shared<I3Map<CCMPMTKey, BaselineEstimate>>();
        for(size_t i = 0; i < num_baselines; ++i) {
            baselines->emplace(pmt_keys_[i], baseline_estimates[i]);
        }
        frame->Put(output_key_, baselines);
        boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> baselines_variance = boost::make_shared<I3Map<CCMPMTKey, BaselineEstimate>>();
        for(size_t i = 0; i < num_baselines; ++i) {
            baselines_variance->emplace(pmt_keys_[i], baseline_estimates_variance[i]);
        }
        frame->Put(output_key_ + "Variance", baselines_variance);
    }

    running.store(false);
}

void RunFrameThread(BaselineEstimatorJob * job,
    std::string const & ccm_waveforms_key_,
    std::string const & output_key_,
    CCMGeometry & geo_,
    std::vector<CCMPMTKey> const & pmt_keys_,
    std::vector<uint32_t> const & channels_,
    bool output_channels_,
    size_t num_samples,
    size_t num_exp_samples) {

    job->running.store(true);

    job->thread = std::thread(FrameThread,
        std::ref(job->running),
        job->frame.get(),
        std::cref(ccm_waveforms_key_),
        std::cref(output_key_),
        std::cref(geo_),
        std::cref(pmt_keys_),
        std::cref(channels_),
        output_channels_,
        num_samples,
        num_exp_samples
    );
}

void BaselineEstimator::DAQ(I3FramePtr frame) {
    if(not geo_seen_) {
        log_fatal("No geometry seen before first DAQ frame");
        PushFrame(frame);
        return;
    }
    while(true) {
        // Check if any jobs have finished
        for(int i=int(running_jobs.size())-1; i>=0; --i) {
            if (teptr) {
                try{
                    std::rethrow_exception(teptr);
                }
                catch(const std::exception &ex)
                {
                    std::cerr << "Thread exited with exception: " << ex.what() << "\n";
                }
            }
            if(not running_jobs[i]->running.load()) {
                BaselineEstimatorJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                results[job->frame_index - min_frame_idx].done = true;
            } else {
                BaselineEstimatorJob * job = running_jobs[i];
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                PushFrame(results[i].frame);
                results[i].frame = nullptr;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_frame_idx += results_done;
        }

        if(not frame)
            break;

        // Attempt to queue up a new job for the frame
        BaselineEstimatorJob * job = nullptr;

        if(free_jobs.size() > 0) {
            job = free_jobs.front();
            job->running.store(false);
            free_jobs.pop_front();
        } else if(running_jobs.size() < num_threads) {
            job = new BaselineEstimatorJob();
            job->running.store(false);
            job->thread_index = running_jobs.size();
        }

        if(job != nullptr and results.size() < max_cached_frames) {
            job->running.store(true);
            running_jobs.push_back(job);
            job->frame = frame;
            job->frame_index = frame_index;
            results.emplace_back();
            results.back().frame = frame;
            results.back().done = false;
            frame_index += 1;
            RunFrameThread(job,
                    ccm_waveforms_key_,
                    output_key_,
                    geo_,
                    pmt_keys_,
                    channels_,
                    output_channels_,
                    num_samples,
                    num_exp_samples);
            break;
        } else if(job != nullptr) {
            free_jobs.push_back(job);
        }
    }
}

void BaselineEstimator::Finish() {
    while(running_jobs.size() > 0) {
        // Check if any jobs have finished
        for(int i=int(running_jobs.size())-1; i>=0; --i) {
            if(not running_jobs[i]->running.load()) {
                BaselineEstimatorJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                results[job->frame_index - min_frame_idx].done = true;
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                PushFrame(results[i].frame);
                results[i].frame = nullptr;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_frame_idx += results_done;
        }
    }
}
