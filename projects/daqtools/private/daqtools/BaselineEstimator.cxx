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

class BaselineEstimator: public I3Module {
    std::exception_ptr teptr = nullptr;
    std::string geometry_name_;
    bool geo_seen;
    std::string ccm_waveforms_name_;
    std::string output_name_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;

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
    void Process();
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

template<typename Iterator>
void OutlierFilter(double starting_estimate, std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end, Iterator result_begin) {
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
    double S2_sum;
    S[0] = 0;
    double sum_S;

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

I3_MODULE(BaselineEstimator);

BaselineEstimator::BaselineEstimator(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key for CCMWaveforms object", std::string("CCMWaveforms"));
    AddParameter("NumThreads", "Number of worker threads to use for baseline estimation", (size_t)(0));
    AddParameter("NumSamples", "Number of samples to use for the baseline estimation", (size_t)(800));
    AddParameter("NumExpSamples", "Number of samples to use for the exponential fit", (size_t)(4000));
    AddParameter("OutputName", "Key to save output I3Vector<BaselineEstimate> to", std::string("BaselineEstimates"));
    AddParameter("MaxCachedFrames", "The maximum number of frames this module is allowed to have cached", (size_t)(1000));
}

void BaselineEstimator::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", ccm_waveforms_name_);
    GetParameter("NumThreads", num_threads);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
    GetParameter("NumSamples", num_samples);
    GetParameter("NumExpSamples", num_exp_samples);
    GetParameter("OutputName", output_name_);
    GetParameter("MaxCachedFrames", max_cached_frames);
}

void BaselineEstimator::Process() {
  if (inbox_)
    log_trace("%zu frames in inbox", inbox_->size());

  I3FramePtr frame = PopFrame();

  if(!frame or frame->GetStop() == I3Frame::DAQ) {
    DAQ(frame);
    return;
  }

  if(frame->GetStop() == I3Frame::Physics && ShouldDoPhysics(frame))
    {
      Physics(frame);
    }
  else if(frame->GetStop() == I3Frame::Geometry && ShouldDoGeometry(frame))
    Geometry(frame);
  else if(frame->GetStop() == I3Frame::Calibration && ShouldDoCalibration(frame))
    Calibration(frame);
  else if(frame->GetStop() == I3Frame::DetectorStatus && ShouldDoDetectorStatus(frame))
    DetectorStatus(frame);
  else if(frame->GetStop() == I3Frame::Simulation && ShouldDoSimulation(frame))
    Simulation(frame);
  else if(frame->GetStop() == I3Frame::DAQ && ShouldDoDAQ(frame)) {
    DAQ(frame);
  } else if(ShouldDoOtherStops(frame))
    OtherStops(frame);
}


void BaselineEstimator::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void FrameThread(std::atomic<bool> & running, I3Frame * frame, std::string const & ccm_waveforms_name_, std::string const & output_name_, I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map_, size_t num_samples, size_t num_exp_samples) {
    // let's read in our waveform
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>(ccm_waveforms_name_);
    if(waveforms == nullptr) {
        log_fatal("No CCMWaveformUInt16Series under key name \"%s\"", ccm_waveforms_name_.c_str());
    }

    std::vector<CCMPMTKey> pmt_keys;
    pmt_keys.reserve(pmt_channel_map_.size());
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        pmt_keys.push_back(p.first);
    }
    size_t num_baselines = pmt_keys.size();
    std::vector<BaselineEstimate> baseline_estimates;
    baseline_estimates.resize(num_baselines);

    // loop over each channel in waveforms
    for(size_t i=0; i<num_baselines; ++i) {
        // let's get our baseline
        CCMPMTKey key = pmt_keys[i];
        uint32_t channel = pmt_channel_map_.at(key);
        ProcessWaveform(
            baseline_estimates[i],
            std::cref(waveforms->at(channel).GetWaveform()),
            num_samples,
            num_exp_samples
        );
    }

    // I3Map to store pmt key and baselines
    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> baselines = boost::make_shared<I3Map<CCMPMTKey, BaselineEstimate>>();
    for(size_t i = 0; i < num_baselines; ++i) {
        baselines->emplace(pmt_keys[i], baseline_estimates[i]);
    }

    frame->Put(output_name_, baselines);
    running.store(false);
}

void RunFrameThread(BaselineEstimatorJob * job,
    std::string const & ccm_waveforms_name_,
    std::string const & output_name_,
    I3Map<CCMPMTKey, uint32_t> const & pmt_channel_map_,
    size_t num_samples,
    size_t num_exp_samples) {

    job->running.store(true);

    job->thread = std::thread(FrameThread,
        std::ref(job->running),
        job->frame.get(),
        std::cref(ccm_waveforms_name_),
        std::cref(output_name_),
        std::cref(pmt_channel_map_),
        num_samples,
        num_exp_samples
    );
}

void BaselineEstimator::DAQ(I3FramePtr frame) {
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
                ccm_waveforms_name_,
                output_name_,
                pmt_channel_map_,
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
