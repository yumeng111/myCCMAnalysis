#include <icetray/IcetrayFwd.h>

#include <iostream>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iterator>
#include <utility>
#include <numeric>
#include <functional>

#include <boost/scoped_ptr.hpp>
#include <boost/make_shared.hpp>

#include <icetray/I3Frame.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <phys-services/I3RandomService.h>

#include "daqtools/OnlineRobustStats.h"
#include "daqtools/WaveformSmoother.h"
#include "dataclasses/calibration/BaselineEstimate.h"

struct WFInfo {
    WaveformSmoother & smoother;
    std::vector<double> values;
    double linear_mean;
    double linear_intercept;
    double linear_slope;
    double linear_error;
    double flat_mean;
    double flat_error;
    double exp_mean;
    double exp_input;
    double exp_sigma;
    bool suitable_for_estimate;
    WFInfo(WFInfo && wfinfo) = default;
    WFInfo(WFInfo & wfinfo) = default;
    WFInfo(WFInfo const & wfinfo) = default;
    WFInfo(WaveformSmoother & s) : smoother(s) {}
    WFInfo & operator=(WFInfo && o) {
        smoother = o.smoother;
        values = std::move(o.values);
        linear_mean = o.linear_mean;
        linear_intercept = o.linear_intercept;
        linear_slope = o.linear_slope;
        linear_error = o.linear_error;
        flat_mean = o.flat_mean;
        flat_error = o.flat_error;
        exp_mean = o.exp_mean;
        exp_input = o.exp_input;
        exp_sigma = o.exp_sigma;
        suitable_for_estimate = o.suitable_for_estimate;
    }
};

struct EstimateInfo {
    I3FramePtr cached_frame;
    WFInfo waveform_properties;
    size_t frames_passed;
    bool contributed_to_estimator;
    size_t contributed_estimator_index;
    EstimateInfo(I3FramePtr frame, WFInfo && info) :
        cached_frame(frame),
        waveform_properties(std::move(info)),
        frames_passed(0),
        contributed_to_estimator(false),
        contributed_estimator_index(0) {}
    EstimateInfo(EstimateInfo && o) :
        cached_frame(std::move(o.cached_frame)),
        waveform_properties(o.waveform_properties),
        frames_passed(std::move(o.frames_passed)),
        contributed_to_estimator(std::move(o.contributed_to_estimator)),
        contributed_estimator_index(std::move(o.contributed_estimator_index)) {};
    EstimateInfo& operator=(EstimateInfo && o) {
        cached_frame = std::move(o.cached_frame);
        waveform_properties = std::move(o.waveform_properties);
        frames_passed = o.frames_passed;
        contributed_to_estimator = o.contributed_to_estimator;
        contributed_estimator_index = o.contributed_estimator_index;
        return *this;
    }
};

struct BaselinePMTKeyInfo {
    std::map<size_t, OnlineRobustStatsBatched> baseline_estimators;
    std::map<size_t, size_t> n_frames_since_last_estimator_update;
    std::deque<EstimateInfo> estimate_info;
    size_t num_estimators = 0;
};

struct BaselineFrameInfo {
    std::set<CCMPMTKey> channels_pending;
    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> baseline_estimates;
    BaselineFrameInfo() :
        channels_pending(),
        baseline_estimates(boost::make_shared<I3Map<CCMPMTKey, BaselineEstimate>>())
    {}
};


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

template <class RandomIt, class U = typename std::iterator_traits<RandomIt>::value_type>
std::vector<U> ChooseNRandom(size_t N, RandomIt begin, RandomIt end, I3RandomServicePtr r) {
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    diff_t n = std::distance(begin, end);
    diff_t max_n = std::min(diff_t(N), n - 1);
    std::vector<U> result(begin, end);
    for(diff_t i=0; i<max_n; ++i) {
        std::swap(result[i], result[r->Integer(n-i) + i]);
    }
    result.resize(N);
    return result;
}

} // namespace

void LinearFlatFit(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end, double & linear_b, double & linear_m, double & linear_sigma, double & flat_mean, double & flat_sigma , double & exp_mean, double & exp_input, double & exp_sigma, double tau) {
    double N = std::distance(begin, end);
    double X = 0.0;
    double X2 = 0.0;
    double T = 0.0;
    double T2 = 0.0;
    double XT = 0.0;

    double E = 0.0;
    double E2 = 0.0;
    double XE = 0.0;

    size_t t = 0;
    std::vector<double>::const_iterator it = begin;
    for(; it != end; ++it) {
        double x = *it;
        X += x;
        X2 += x*x;
        T += t;
        T2 += t*t;
        XT += t*x;

        double e = exp(-double(t)/(tau));
        E += e;
        E2 += e*e;
        XE += x * e;

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

    exp_mean = (-E2 * X + E * XE) / (E*E - E2 * N);
    exp_input = (E * X - N * XE) / (E*E - E2 * N);
    exp_sigma = sqrt(
            (-2.0 * E2 * X + E2 * X * X - E * E * X2 + E2 * N * X2 + 2 * E * XE - 4.0 * E * X * XE + 3 * N * XE * XE) /
            (- E * E * N + E2 * N * N)
            );
}

BaselineEstimate SingleWFBaselineEstimate(WFInfo & wf_info, double derivative_threshold, size_t target_size, size_t min_size, double tau, bool use_exp_fit) {
    WaveformSmoother & smoother = wf_info.smoother;
    smoother.Reset();
    size_t largest_region_size = 0;
    int largest_region_idx = -1;
    std::vector<std::pair<size_t, size_t>> region_idxs;
    while(true) {
        bool in_region = false;
        size_t first_idx = 0;
        size_t end_idx = 0;
        for(size_t i=0; i<smoother.Size(); ++i) {
            double deriv_value = smoother.Derivative();
            bool within_thresh = std::abs(deriv_value) <= derivative_threshold;
            if(within_thresh) {
                if(in_region) {
                    end_idx += 1;
                } else {
                    first_idx = i;
                }
                in_region = true;
            } else {
                if(in_region) {
                    end_idx = i;
                    region_idxs.emplace_back(first_idx, end_idx);
                    size_t size = end_idx - first_idx;
                    if(size > largest_region_size) {
                        largest_region_size = size;
                        largest_region_idx = region_idxs.size() - 1;
                    }
                    if(size > target_size)
                        break;
                }
                in_region = false;
            }
            smoother.Next();
        }
        if(largest_region_size < min_size) {
            derivative_threshold *= 1.2;
            largest_region_size = 0;
            largest_region_idx = -1;
            region_idxs.clear();
            continue;
        } else {
            break;
        }
    }
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> wf_its = smoother.GetSmoothedWaveform();
    std::vector<double>::const_iterator baseline_begin = wf_its.first + region_idxs[largest_region_idx].first;
    std::vector<double>::const_iterator baseline_end = wf_its.first + region_idxs[largest_region_idx].second;

    LinearFlatFit(baseline_begin, baseline_end,
            wf_info.linear_intercept,
            wf_info.linear_slope,
            wf_info.linear_error,
            wf_info.flat_mean,
            wf_info.flat_error,
            wf_info.exp_mean,
            wf_info.exp_input,
            wf_info.exp_sigma,
            tau);

    wf_info.values = std::vector<double>(baseline_begin, baseline_end);
    if(use_exp_fit) {
        for(size_t i=0; i<wf_info.values.size(); ++i) {
            wf_info.values[i] -= wf_info.exp_input * exp(-double(i)/tau);
        }
    }
    std::sort(wf_info.values.begin(), wf_info.values.end());
    double mode = robust_stats::Mode(wf_info.values.begin(), wf_info.values.end());

    BaselineEstimate baseline_estimate(
            mode,
            robust_stats::MedianAbsoluteDeviation(wf_info.values.begin(), wf_info.values.end(), mode),
            0,
            1,
            wf_info.values.size());
    return baseline_estimate;
}

double BaselineComparisonScore(OnlineRobustStatsBatched & stats, std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) {
    //double mode = stats.Mode();
    double mode = stats.Mean();
    double stddev = stats.Stddev(mode);
    size_t N = std::distance(begin, end);
    double chi2 = 0.0;
    for(; begin != end; ++begin) {
        double z = ((*begin) - mode) / stddev;
        chi2 += z;
    }
    chi2 /= 2.0;
    chi2 /= N;
    chi2 += log(stddev) + log(sqrt(2.0 * M_PI));
    return chi2 /= N;
}

size_t ChooseEstimator(WFInfo const & wf_info, BaselinePMTKeyInfo & pmt_key_info, bool allow_new_estimator=false, double max_score=1.5) {
    std::map<size_t, OnlineRobustStatsBatched> & baseline_estimators = pmt_key_info.baseline_estimators;
    if(baseline_estimators.size() == 0) {
        return 0;
    }
    std::map<size_t, double> comparison_scores;
    for(std::pair<size_t const, OnlineRobustStatsBatched> & p : baseline_estimators) {
        comparison_scores[p.first] = BaselineComparisonScore(p.second, wf_info.values.cbegin(), wf_info.values.cend());
    }
    std::map<size_t, double>::iterator min_score =
        std::min_element(comparison_scores.begin(), comparison_scores.end(),
                [](decltype(comparison_scores)::value_type& l, decltype(comparison_scores)::value_type& r) -> bool { return l.second < r.second; });
    size_t best_idx = min_score->first;
    if(allow_new_estimator and min_score->second > max_score) {
        return pmt_key_info.num_estimators;
    } else {
        return best_idx;
    }
}

WFInfo ComputeWFInfo(WaveformSmoother & smoother, std::pair<size_t, size_t> & pulse_positions, I3RandomServicePtr random, size_t num_samples, double tau, bool use_exp_fit) {
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> wf_its = smoother.GetSmoothedWaveform();
    std::vector<double>::const_iterator wf_begin = wf_its.first;
    std::vector<double>::const_iterator wf_end = wf_its.second;
    size_t wf_size = std::distance(wf_begin, wf_end);
    std::vector<double>::const_iterator baseline_begin = wf_begin + std::max(ptrdiff_t(std::min(ptrdiff_t(200), ptrdiff_t(wf_size)-1)), ptrdiff_t(0));
    size_t N_baseline_samples = size_t(std::max(ptrdiff_t(0), ptrdiff_t(std::min(ptrdiff_t(std::distance(baseline_begin, wf_end)), ptrdiff_t(pulse_positions.first))) - 5));
    std::vector<double>::const_iterator baseline_end = baseline_begin + N_baseline_samples;

    WFInfo wf_info(smoother);
    if(N_baseline_samples == 0) {
        std::cout << "No baseline samples" << std::endl;
        wf_info.linear_mean = 0;
        wf_info.linear_intercept = 0;
        wf_info.linear_slope = 0;
        wf_info.linear_error = 0;
        wf_info.flat_mean = 0;
        wf_info.flat_error = 0;
        wf_info.suitable_for_estimate = false;
        return wf_info;
    }

    //wf_info.values = ChooseNRandom(num_samples, baseline_begin, baseline_end, random);
    //std::copy(baseline_begin, baseline_end, std::back_inserter(wf_info.values));

    LinearFlatFit(baseline_begin, baseline_end,
            wf_info.linear_intercept,
            wf_info.linear_slope,
            wf_info.linear_error,
            wf_info.flat_mean,
            wf_info.flat_error,
            wf_info.exp_mean,
            wf_info.exp_input,
            wf_info.exp_sigma,
            tau);

    wf_info.values = std::vector<double>(baseline_begin, baseline_end);
    if(use_exp_fit) {
        for(size_t i=0; i<wf_info.values.size(); ++i) {
            wf_info.values[i] -= wf_info.exp_input * exp(-double(i)/tau);
        }
    }

    wf_info.values = ChooseNRandom(num_samples, baseline_begin, baseline_end, random);

    wf_info.linear_mean = wf_info.linear_intercept + wf_info.linear_slope * N_baseline_samples / 2.0;
    wf_info.suitable_for_estimate = true;

    if(std::abs(wf_info.linear_mean - wf_info.flat_mean) / std::min(wf_info.linear_error, wf_info.flat_error) > 1) {
        std::cout << "Linear mean does not match flat mean: " << wf_info.linear_mean << "+/-" << wf_info.linear_error << " vs " <<  wf_info.flat_mean << "+/-" << wf_info.flat_error << std::endl;
        wf_info.suitable_for_estimate = false;
    }
    if(wf_info.flat_error / wf_info.linear_error > 1.5) {
        std::cout << "Flat error much larger than linear error: " << wf_info.flat_error << " vs " << wf_info.linear_error << std::endl;
        wf_info.suitable_for_estimate = false;
    }
    if(std::abs(wf_info.linear_slope * N_baseline_samples) > 2.0) {
        std::cout << "Delta ADC is too large: " << wf_info.linear_slope * N_baseline_samples << std::endl;
        wf_info.suitable_for_estimate = false;
    }
    if(std::abs(wf_info.linear_slope) > 0.2) {
        std::cout << "Linear slope is too large: " << wf_info.linear_slope << std::endl;
        wf_info.suitable_for_estimate = false;
    }
    return wf_info;
}

class BaselineEstimator : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string input_name_;
    std::string output_name_;

    size_t num_frames_for_estimate;
    size_t num_samples;
    size_t max_frames_waiting_for_estimate;
    size_t max_estimator_lifetime;

    double droop_tau;
    bool use_exp_fit;

    size_t frames_seen = 0;

    // Internal state
    bool geo_seen;
    boost::shared_ptr<CCMGeometry const> geo;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;

    std::deque<I3FramePtr> cached_frames;
    std::map<CCMPMTKey, BaselinePMTKeyInfo> baseline_pmt_info;
    std::map<I3FramePtr, BaselineFrameInfo> baseline_frame_info;

    public:
    BaselineEstimator(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    void UpdateEstimators(I3FramePtr frame);
    void UpdateEstimates();
    void UpdateAndPushCompletedFrames();
    void RemoveOldEstimators();
};

I3_MODULE(BaselineEstimator);

BaselineEstimator::BaselineEstimator(const I3Context& context) : I3Module(context),
    geometry_name_(""),  geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("Output", "Ouput key prefix", std::string(""));
    AddParameter("Input", "Input key prefix", std::string(""));
    AddParameter("NumFramesForEstimate", "Number of frames to use for baseline estimate", size_t(100));
    AddParameter("NumSamples", "Number of samples from each frame to use for baseline estimate", size_t(200));
    AddParameter("MaxFramesWaitingForEstimate", "Maximum number of frames to postpone computing an estimate for frame that has an incomplete estimator", size_t(500));
    AddParameter("MaxEstimatorLifetime", "Maximum number of frames since receiving an estimate that we should keep an estimator", size_t(1000));
    AddParameter("DroopTau", "Time constant for droop correction", double(3e3));
    AddParameter("UseExpFit", "Use the exponential fit to correct for droop effects", bool(false));
}

void BaselineEstimator::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("Output", output_name_);
    GetParameter("Input", input_name_);
    GetParameter("NumFramesForEstimate", num_frames_for_estimate);
    GetParameter("NumSamples", num_samples);
    GetParameter("MaxFramesWaitingForEstimate", max_frames_waiting_for_estimate);
    GetParameter("MaxEstimatorLifetime", max_estimator_lifetime);
    GetParameter("DroopTau", droop_tau);
    GetParameter("UseExpFit", use_exp_fit);
}

void BaselineEstimator::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    geo_seen = true;
    geo = frame->Get<boost::shared_ptr<CCMGeometry const>>(geometry_name_);


    // Cache the trigger channel map
    pmt_channel_map_ = geo->pmt_channel_map;
    PushFrame(frame);
}

void BaselineEstimator::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    ++frames_seen;
    UpdateEstimators(frame);
    UpdateEstimates();
    UpdateAndPushCompletedFrames();
    RemoveOldEstimators();
}

void BaselineEstimator::UpdateEstimators(I3FramePtr frame) {
    I3RandomServicePtr random = GetService<I3RandomServicePtr>("I3RandomService");
    I3Map<CCMPMTKey, WaveformSmoother> smoothed_wfs = frame->Get<I3Map<CCMPMTKey, WaveformSmoother>>("WaveformSmoothers");
    I3Map<CCMPMTKey, std::pair<size_t, size_t>> pulse_positions = frame->Get<I3Map<CCMPMTKey, std::pair<size_t, size_t>>>("PulsePositions");

    std::set<CCMPMTKey> & channels_pending = baseline_frame_info[frame].channels_pending;
    cached_frames.push_back(frame);

    for(std::pair<CCMPMTKey const, WaveformSmoother> & p : smoothed_wfs) {
        CCMPMTKey pmt_key = p.first;
        CCMOMGeo::OMType type = geo->pmt_geo.at(pmt_key).omtype;
        if(type == CCMOMGeo::OMType::BeamCurrentMonitor) {
            continue;
        }
        WaveformSmoother & smoother = p.second;
        BaselinePMTKeyInfo & pmt_info = baseline_pmt_info[pmt_key];
        std::deque<EstimateInfo> & estimate_info = baseline_pmt_info[pmt_key].estimate_info;
        std::map<size_t, OnlineRobustStatsBatched> & baseline_estimators = baseline_pmt_info[pmt_key].baseline_estimators;
        std::map<size_t, size_t> & n_frames_since_last_estimator_update = baseline_pmt_info[pmt_key].n_frames_since_last_estimator_update;

        std::pair<size_t, size_t> & pos = pulse_positions.at(pmt_key);
        EstimateInfo pending_estimate(frame, ComputeWFInfo(smoother, pos, random, num_samples, droop_tau / 2.0, use_exp_fit));
        WFInfo & wf_info = pending_estimate.waveform_properties;
        //std::cout << pmt_key << " is suitable for estimate? " << (wf_info.suitable_for_estimate ? "Yes" : "No") << std::endl;
        if(wf_info.suitable_for_estimate) {
            // Choose from existing estimator
            // Return the end index if an appropriate estimator does not exist
            // or if there are no estimators
            size_t estimator_index = ChooseEstimator(wf_info, pmt_info, true);
            //std::cout << "\tChose estimator[" << estimator_index << "]" << std::endl;
            if(estimator_index >= pmt_info.num_estimators) {
                // Instantiate a new estimator
                baseline_estimators[pmt_info.num_estimators];
                n_frames_since_last_estimator_update[pmt_info.num_estimators] = 0;
                pmt_info.num_estimators += 1;
            }
            OnlineRobustStatsBatched & estimator = baseline_estimators.at(estimator_index);
            // Add samples to the estimator
            estimator.AddValues(wf_info.values);
            pending_estimate.contributed_to_estimator = true;
            pending_estimate.contributed_estimator_index = estimator_index;
            if(estimator.NBatches() > num_frames_for_estimate) {
                estimator.RemoveValue();
            }
            // Reset counter for current estimator
            n_frames_since_last_estimator_update[estimator_index] = 0;
            // Update the counter for the number of frames since the last estimator update
            for(std::pair<size_t const, OnlineRobustStatsBatched> & p_est : baseline_estimators) {
                if(p_est.first != estimator_index) {
                    n_frames_since_last_estimator_update[p_est.first] += 1;
                }
            }
        }
        //for(std::pair<size_t const, OnlineRobustStatsBatched> & p_est : baseline_estimators) {
        //    std::cout << "\testimator[" << p_est.first << "]: " << p_est.second.NBatches() << " frames, " << p_est.second.NSamples() << " samples" << std::endl;
        //}
        estimate_info.emplace_back(std::move(pending_estimate));
        channels_pending.insert(pmt_key);
    }
}

void BaselineEstimator::RemoveOldEstimators() {
    for(std::pair<CCMPMTKey const, BaselinePMTKeyInfo> & p : baseline_pmt_info) {
        std::map<size_t, OnlineRobustStatsBatched> & baseline_estimators = p.second.baseline_estimators;
        std::map<size_t, size_t> & n_frames_since_last_estimator_update = p.second.n_frames_since_last_estimator_update;
        std::vector<size_t> estimators_to_remove;
        for(std::pair<size_t const, size_t> & it : n_frames_since_last_estimator_update) {
            if(it.second > max_estimator_lifetime) {
                estimators_to_remove.push_back(it.first);
            }
        }
        for(size_t const & est_idx : estimators_to_remove) {
            baseline_estimators.erase(est_idx);
        }
    }
}

void BaselineEstimator::UpdateEstimates() {
    for(std::pair<CCMPMTKey const, BaselinePMTKeyInfo> & p : baseline_pmt_info) {
        CCMPMTKey pmt_key = p.first;
        BaselinePMTKeyInfo & pmt_info = p.second;
        std::deque<EstimateInfo> & estimate_info = pmt_info.estimate_info;
        //std::cout << pmt_key << " has " << estimate_info.size() << " pending estimates." << std::endl;
        std::map<size_t, OnlineRobustStatsBatched> & baseline_estimators = pmt_info.baseline_estimators;

        std::vector<size_t> estimates_to_remove;
        size_t estimate_index = 0;
        bool skip = false;
        // Check if cached frames can have estimates added
        for(EstimateInfo & estimate : estimate_info) {
            I3FramePtr cached_frame = estimate.cached_frame;
            boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> baseline_estimates = baseline_frame_info[cached_frame].baseline_estimates;
            std::set<CCMPMTKey> & channels_pending = baseline_frame_info[cached_frame].channels_pending;
            WFInfo & wf_info = estimate.waveform_properties;
            if(baseline_estimators.size() == 0) {
                // We can't produce an estimate right now
                //std::cout << "\tNo estimators for " << estimate_index << std::endl;
                estimate.frames_passed += 1;
                estimate_index += 1;
                skip = true;
                continue;
            }
            if(skip and estimate.frames_passed <= max_frames_waiting_for_estimate) {
                // A previous estimate failed and this estimate is not old enough to force an estimate
                // We will assume that this estimate cannot be computed until the previous one succeeds
                estimate.frames_passed += 1;
                estimate_index += 1;
                continue;
            }
            size_t estimator_index;
            bool have_estimator;
            if(estimate.contributed_to_estimator) {
                estimator_index = estimate.contributed_estimator_index;
                have_estimator = true;
            } else {
                // Choose the best estimator, requiring that an existing one be chosen
                estimator_index = ChooseEstimator(wf_info, p.second, false);
                have_estimator = estimator_index < pmt_info.num_estimators;
            }
            if(have_estimator) {
                OnlineRobustStatsBatched & estimator = baseline_estimators.at(estimator_index);
                if(estimator.NBatches() >= num_frames_for_estimate) {
                    //std::cout << "\tAdding estimate for " << estimate_index << std::endl;
                    // If the estimator has enough samples then perform the estimate
                    // Compute the estimator
                    double current_baseline_estimate = estimator.Mode();
                    std::vector<double> baseline_samples;
                    baseline_samples.reserve(wf_info.values.size());
                    for(size_t i=0; i<wf_info.values.size(); ++i) {
                        double delta = (wf_info.values[i] - current_baseline_estimate);
                        double e = std::abs(delta) / 10.0;
                        if(i > 0) {
                            e += std::abs(wf_info.values[i] - wf_info.values[i-1]);
                        }
                        if(i + 1 < wf_info.values.size()) {
                            e += std::abs(wf_info.values[i] - wf_info.values[i+1]);
                        }
                        delta = delta * exp(-e);
                        current_baseline_estimate += delta;
                        baseline_samples.push_back(current_baseline_estimate);
                    }
                    std::sort(baseline_samples.begin(), baseline_samples.end());
                    double mode = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());
                    double stddev = robust_stats::MedianAbsoluteDeviation(baseline_samples.begin(), baseline_samples.end(), mode);
                    BaselineEstimate baseline_estimate(
                            mode,
                            stddev,
                            //estimator.Mode(),
                            //estimator.Stddev(estimator.Mode()),
                            //estimator.Mean(),
                            //estimator.Stddev(estimator.Mean()),
                            num_frames_for_estimate,
                            estimator.NBatches(),
                            estimator.NSamples());
                    // Put the estimate into the map
                    baseline_estimates->insert({pmt_key, baseline_estimate});
                    channels_pending.erase(pmt_key);
                    estimates_to_remove.push_back(estimate_index);
                } else if(estimate.frames_passed > max_frames_waiting_for_estimate) {
                    //std::cout << "\tAdding estimate for " << estimate_index << std::endl;
                    // If too many frames have gone by without an estimate for this frame
                    // then perform an estimate anyway
                    double current_baseline_estimate = estimator.Mode();
                    std::vector<double> baseline_samples;
                    baseline_samples.reserve(wf_info.values.size());
                    for(size_t i=0; i<wf_info.values.size(); ++i) {
                        double delta = (wf_info.values[i] - current_baseline_estimate);
                        double e = std::abs(delta) / 10.0;
                        if(i > 0) {
                            e += std::abs(wf_info.values[i] - wf_info.values[i-1]);
                        }
                        if(i + 1 < wf_info.values.size()) {
                            e += std::abs(wf_info.values[i] - wf_info.values[i+1]);
                        }
                        delta = delta * exp(-e);
                        current_baseline_estimate += delta;
                    }
                    std::sort(baseline_samples.begin(), baseline_samples.end());
                    double mode = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());
                    double stddev = robust_stats::MedianAbsoluteDeviation(baseline_samples.begin(), baseline_samples.end(), mode);
                    BaselineEstimate baseline_estimate(
                            estimator.Mode(),
                            estimator.Stddev(estimator.Mode()),
                            //estimator.Mean(),
                            //estimator.Stddev(estimator.Mean()),
                            num_frames_for_estimate,
                            estimator.NBatches(),
                            estimator.NSamples());
                    baseline_estimates->insert({pmt_key, baseline_estimate});
                    channels_pending.erase(pmt_key);
                    estimates_to_remove.push_back(estimate_index);
                } else {
                    //std::cout << "\tNot enough frames for estimate " << estimate_index << std::endl;
                    // We can't produce an estimate right now
                    estimate.frames_passed += 1;
                    skip = true;
                }
            } else {
                BaselineEstimate baseline_estimate = SingleWFBaselineEstimate(wf_info, 0.3, num_samples * num_frames_for_estimate, 400, droop_tau, use_exp_fit);
                baseline_estimates->insert({pmt_key, baseline_estimate});
                channels_pending.erase(pmt_key);
                estimates_to_remove.push_back(estimate_index);
            }
            estimate_index += 1;
        }
        //std::cout << "Finished " << estimates_to_remove.size()  << " estimates" << std::endl;
        if(estimates_to_remove.size() > 0) {
            // Remove estimates in reverse order
            size_t remove_idx = estimates_to_remove.size();
            size_t idx = estimate_info.size();
            while(remove_idx > 0 and idx > 0) {
                if((idx-1) == estimates_to_remove[remove_idx-1]) {
                    estimate_info.erase(estimate_info.begin() + estimates_to_remove[remove_idx-1]);
                    --remove_idx;
                }
                --idx;
            }
        }
    }
}

void BaselineEstimator::UpdateAndPushCompletedFrames() {
    while(cached_frames.size() > 0) {
        I3FramePtr frame = cached_frames.front();
        bool frame_complete = false;
        {
            BaselineFrameInfo & frame_info = baseline_frame_info[frame];
            frame_complete = (frame_info.channels_pending.size() == 0);
            if(frame_complete) {
                frame->Put((output_name_ + "BaselineEstimates").c_str(), frame_info.baseline_estimates);
                PushFrame(frame);
            }
        }
        if(frame_complete) {
            baseline_frame_info.erase(frame);
            cached_frames.pop_front();
        } else {
            break;
        }
    }
}

void BaselineEstimator::Finish() {
    Flush();
}
