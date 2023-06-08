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

struct WFInfo {
    std::vector<double> values;
    double linear_mean;
    double linear_intercept;
    double linear_slope;
    double linear_error;
    double flat_mean;
    double flat_error;
    bool suitable_for_estimate;
};

struct EstimateInfo {
    I3FramePtr cached_frame;
    WFInfo waveform_properties;
    size_t frames_passed;
    EstimateInfo(I3FramePtr frame, WFInfo info) :
        cached_frame(frame), waveform_properties(info), frames_passed(0) {}
};

struct BaselinePMTKeyInfo {
    std::deque<OnlineRobustStatsBatched> baseline_estimators;
    std::deque<size_t> n_frames_since_last_estimator_update;
    std::deque<EstimateInfo> estimate_info;
};

struct BaselineEstimate {
    double baseline;
    double stddev;
    size_t target_num_frames;
    size_t num_frames;
    size_t num_samples;
    BaselineEstimate() {}
    BaselineEstimate(double baseline, double stddev, size_t target_num_frames, size_t num_frames, size_t num_samples) :
        baseline(baseline), stddev(stddev), target_num_frames(target_num_frames), num_frames(num_frames), num_samples(num_samples) {}
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

} // namespace

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
        flat_sigma += (*it) - flat_mean;
    }
    flat_sigma /= N;
}

double BaselineComparisonScore(OnlineRobustStatsBatched & stats, std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) {
    double mode = stats.Mode();
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

size_t ChooseEstimator(WFInfo const & wf_info, std::deque<OnlineRobustStatsBatched> & baseline_estimators, bool allow_new_estimator=false, double max_score=1.5) {
    std::vector<double> comparison_scores(baseline_estimators.size());
    for(size_t i=0; i<baseline_estimators.size(); ++i) {
        comparison_scores[i] = BaselineComparisonScore(baseline_estimators[i], wf_info.values.cbegin(), wf_info.values.cend());
    }
    std::vector<double>::iterator min_score = std::min_element(comparison_scores.begin(), comparison_scores.end());
    size_t best_idx = std::distance(comparison_scores.begin(), min_score);
    if(allow_new_estimator and (*min_score) > max_score) {
        return baseline_estimators.size() + 1;
    } else {
        return best_idx;
    }
}

WFInfo ComputeWFInfo(WaveformSmoother & smoother, std::pair<size_t, size_t> & pulse_positions) {
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> wf_its = smoother.GetSmoothedWaveform();
    std::vector<double>::const_iterator wf_begin = wf_its.first;
    std::vector<double>::const_iterator wf_end = wf_its.second;
    size_t wf_size = std::distance(wf_begin, wf_end);
    std::vector<double>::const_iterator baseline_begin = wf_begin + std::max(ptrdiff_t(std::min(ptrdiff_t(200), ptrdiff_t(wf_size)-1)), ptrdiff_t(0));
    size_t N_baseline_samples = size_t(std::max(ptrdiff_t(0), ptrdiff_t(std::min(ptrdiff_t(std::distance(baseline_begin, wf_end)), ptrdiff_t(pulse_positions.first))) - 5));
    std::vector<double>::const_iterator baseline_end = baseline_begin + N_baseline_samples;

    WFInfo wf_info;
    if(N_baseline_samples == 0) {
        wf_info.linear_mean = 0;
        wf_info.linear_intercept = 0;
        wf_info.linear_slope = 0;
        wf_info.linear_error = 0;
        wf_info.flat_mean = 0;
        wf_info.flat_error = 0;
        wf_info.suitable_for_estimate = false;
        return wf_info;
    }

    std::copy(baseline_begin, baseline_end, std::back_inserter(wf_info.values));

    LinearFlatFit(baseline_begin, baseline_end,
            wf_info.linear_intercept, 
            wf_info.linear_slope,
            wf_info.linear_error, 
            wf_info.flat_mean,
            wf_info.flat_error);

    wf_info.linear_mean = wf_info.linear_intercept + wf_info.linear_slope * N_baseline_samples / 2.0;

    if(std::abs(wf_info.linear_mean - wf_info.flat_mean) / std::min(wf_info.linear_error, wf_info.flat_error) > 1) {
        wf_info.suitable_for_estimate = false;
    } else if(wf_info.flat_error / wf_info.linear_error > 1.5) {
        wf_info.suitable_for_estimate = false;
    } else if(std::abs(wf_info.linear_slope * N_baseline_samples) > 2.0) {
        wf_info.suitable_for_estimate = false;
    } else if(std::abs(wf_info.linear_slope) > 0.2) {
        wf_info.suitable_for_estimate = false;
    } else {
        wf_info.suitable_for_estimate = true;
    }
    return wf_info;
}

class BaselineEstimator : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string input_name_;

    size_t num_frames_for_estimate;
    size_t num_samples;
    size_t max_frames_waiting_for_estimate;

    // Internal state
    bool geo_seen;
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
    void Flush();

    void UpdateEstimators(I3FramePtr frame);
    void UpdateEstimates();
    void UpdateAndPushCompletedFrames();
    void RemoveOldEstimators();
};

I3_MODULE(BaselineEstimator);

BaselineEstimator::BaselineEstimator(const I3Context& context) : I3Module(context),
    geometry_name_(""),  geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
        AddParameter("Input", "Input key prefix", std::string(""));
        AddParameter("NumFramesForEstimate", "Number of frames to use for baseline estimate", size_t(100));
        AddParameter("NumSamples", "Number of samples from each frame to use for baseline estimate", size_t(200));
        AddParameter("MaxFramesWaitingForEstimate", "Maximum number of frames to postpone computing an estimate for frame that has an incomplete estimator", size_t(500));
    }

void BaselineEstimator::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("Input", input_name_);
    GetParameter("NumFramesForEstimate", num_frames_for_estimate);
    GetParameter("NumSamples", num_samples);
    GetParameter("MaxFramesWaitingForEstimate", max_frames_waiting_for_estimate);
}

void BaselineEstimator::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);


    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;

    geo_seen = true;
    PushFrame(frame);
}

void BaselineEstimator::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    UpdateEstimators(frame);
    UpdateEstimates();
    UpdateAndPushCompletedFrames();
}

void BaselineEstimator::UpdateEstimators(I3FramePtr frame) {
    I3Map<CCMPMTKey, WaveformSmoother> smoothed_wfs = frame->Get<I3Map<CCMPMTKey, WaveformSmoother>>("WaveformSmoothers");
    I3Map<CCMPMTKey, std::pair<size_t, size_t>> pulse_positions = frame->Get<I3Map<CCMPMTKey, std::pair<size_t, size_t>>>("PulsePositions");

    std::set<CCMPMTKey> & channels_pending = baseline_frame_info[frame].channels_pending;
    cached_frames.push_back(frame);

    for(std::pair<CCMPMTKey const, WaveformSmoother> & p : smoothed_wfs) {
        CCMPMTKey pmt_key = p.first;
        WaveformSmoother & smoother = p.second;
        std::deque<EstimateInfo> & estimate_info = baseline_pmt_info[pmt_key].estimate_info;
        std::deque<OnlineRobustStatsBatched> & baseline_estimators = baseline_pmt_info[pmt_key].baseline_estimators;
        std::deque<size_t> & n_frames_since_last_estimator_update = baseline_pmt_info[pmt_key].n_frames_since_last_estimator_update;

        std::pair<size_t, size_t> & pos = pulse_positions.at(pmt_key);
        WFInfo wf_info = ComputeWFInfo(smoother, pos);
        if(wf_info.suitable_for_estimate) {
            // Choose from existing estimator
            // Return the end index if an appropriate estimator does not exist
            // or if there are no estimators
            size_t estimator_index = ChooseEstimator(wf_info, baseline_estimators, true);
            if(estimator_index > baseline_estimators.size()) {
                // Instantiate a new estimator
                baseline_estimators.emplace_back();
                n_frames_since_last_estimator_update.push_back(0);
            }
            OnlineRobustStatsBatched & estimator = baseline_estimators[estimator_index];
            // Add samples to the estimator
            estimator.AddValues(wf_info.values);
            if(estimator.NBatches() > num_frames_for_estimate) {
                estimator.RemoveValue();
            }
            // Reset counter for current estimator
            n_frames_since_last_estimator_update[estimator_index] = 0;
            // Update the counter for the number of frames since the last estimator update
            for(size_t i=0; i<baseline_estimators.size(); ++i) {
                if(i != estimator_index) {
                    n_frames_since_last_estimator_update[i] += 1;
                }
            }
        }
        estimate_info.emplace_back(frame, wf_info);
        channels_pending.insert(pmt_key);
    }
}

void BaselineEstimator::UpdateEstimates() {
    for(std::pair<CCMPMTKey const, BaselinePMTKeyInfo> & p : baseline_pmt_info) {
        CCMPMTKey pmt_key = p.first;
        BaselinePMTKeyInfo & pmt_info = p.second;
        std::deque<EstimateInfo> & estimate_info = pmt_info.estimate_info;
        std::deque<OnlineRobustStatsBatched> & baseline_estimators = pmt_info.baseline_estimators;
        std::deque<size_t> & n_frames_since_last_estimator_update = pmt_info.n_frames_since_last_estimator_update;

        std::vector<size_t> estimates_to_remove;
        size_t estimate_index = 0;
        // Check if cached frames can have estimates added
        for(EstimateInfo & estimate : estimate_info) {
            I3FramePtr cached_frame = estimate.cached_frame;
            boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> baseline_estimates = baseline_frame_info[cached_frame].baseline_estimates;
            std::set<CCMPMTKey> & channels_pending = baseline_frame_info[cached_frame].channels_pending;
            WFInfo & wf_info = estimate.waveform_properties;
            // Choose the best estimator, requiring that an existing one be chosen
            size_t estimator_index = ChooseEstimator(wf_info, baseline_estimators, false);
            OnlineRobustStatsBatched & estimator = baseline_estimators[estimator_index];
            if(estimator.NBatches() >= num_frames_for_estimate) {
                // If the estimator has enough samples then perform the estimate
                // Compute the estimator
                BaselineEstimate baseline_estimate(
                        estimator.Mode(),
                        estimator.Stddev(estimator.Mode()),
                        num_frames_for_estimate,
                        estimator.NBatches(),
                        estimator.NSamples());
                // Put the estimate into the map
                baseline_estimates->insert({pmt_key, baseline_estimate});
                channels_pending.erase(pmt_key);
                estimates_to_remove.push_back(estimate_index);
            } else if(estimate.frames_passed > max_frames_waiting_for_estimate) {
                // If too many frames have gone by without an estimate for this frame
                // then perform an estimate anyway
                BaselineEstimate baseline_estimate(
                        estimator.Mode(),
                        estimator.Stddev(estimator.Mode()),
                        num_frames_for_estimate,
                        estimator.NBatches(),
                        estimator.NSamples());
                baseline_estimates->insert({pmt_key, baseline_estimate});
                channels_pending.erase(pmt_key);
                estimates_to_remove.push_back(estimate_index);
            } else {
                // We can't produce an estimate right now
                estimate.frames_passed += 1;
            }
            estimate_index += 1;
        }
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
                frame->Put("BaselineEstimates", frame_info.baseline_estimates);
                PushFrame(frame);
            }
        }
        if(frame_complete) {
            baseline_frame_info.erase(frame);
            cached_frames.pop_front();
        } else {
            continue;
        }
    }
}

void BaselineEstimator::Flush() {
}

void BaselineEstimator::Finish() {
    Flush();
}
