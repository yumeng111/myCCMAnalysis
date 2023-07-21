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
    size_t DOF = std::max(N, 3) - 2;
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
    size_t DOF = std::max(N, 4) - 3;
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

class darcy_baselines: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    void ProcessWaveform(std::vector<short unsigned int> const & samples, BaselineEstimate & baseline);
    void FitExponential(std::vector<double> const & y, double & a, double & b, double & c);
    void OutlierFilter(std::vector<short unsigned int> const & samples, std::vector<double> & outlier_filter_results);
    public:
    darcy_baselines(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
};

I3_MODULE(darcy_baselines);

darcy_baselines::darcy_baselines(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    }


void darcy_baselines::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
}


void darcy_baselines::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void darcy_baselines::OutlierFilter(std::vector<short unsigned int> const & samples, std::vector<double> & outlier_filter_results){

    // first let's find the average of the first 10 bins of our wf as the starting value
    double delta_tau = 20;
    double prev_tau = 2.0;
    double next_tau = 2.0;
    double starting_val = 0;
    double counter = 0;

    for(size_t it = 0; it < 10; ++it) {
        starting_val += samples[it];
        counter += 1;
    }

    double value = starting_val/counter;

    // now let's loop over the waveform
    for(size_t wf_it = 0; wf_it < samples.size(); ++wf_it) {
        double delta = samples[wf_it] - value;
        double e = std::fabs(delta) / delta_tau;

        if(wf_it > 0) {
            e += std::fabs(samples[wf_it] - samples[wf_it - 1]) / prev_tau;
        }

        if(wf_it + 1 < samples.size()) {
            e += std::fabs(samples[wf_it] - samples[wf_it + 1]) / next_tau;
        }
        delta *= std::exp(-e);
        value += delta;
        outlier_filter_results[wf_it] = value;
    }

}

void darcy_baselines::FitExponential(std::vector<double> const & y, double & a, double & b, double & c){

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

    c = mult_1; // one of the parameters we're solving for!!!

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

    a = m_inv00 * v0 + m_inv01 * v1;
    b = m_inv01 * v0 + m_inv11 * v1;
}


void darcy_baselines::ProcessWaveform(std::vector<short unsigned int> const & samples, BaselineEstimate & baseline){

    if (samples.size() == 0) {
        return;
    }

    // vector to store results of outlier filter
    std::vector<double> outlier_filter_results(samples.size());
    OutlierFilter(samples, outlier_filter_results);

    double flat_mean = 0;
    double flat_sigma = 0;
    double linear_intercept = 0;
    double linear_slope = 0;
    double linear_sigma = 0;
    LinearFlatFit(outlier_filter_results.begin(), outlier_filter_results.end(), linear_intercept, linear_slope, linear_sigma, flat_mean, flat_sigma);

    double flat_chi2_per_dof_threshold = 1.0;
    double flat_sigma_threshold = 30.0;
    double flat_chi2_per_dof = FlatChi2PerDOF(outlier_filter_results.begin(), outlier_filter_results.end(), flat_mean, flat_sigma);
    if(flat_chi2_per_dof < flat_chi2_per_dof_threshold and flat_sigma < flat_sigma_threshold) {
        // a flat fit is good enough
        // Let's use the mode of the outlier filter as the baseline estimate
    
        return;
    }

    double linear_chi2_per_dof_threshold = 1.0;
    double linear_sigma_threshold = 10.0;
    double linear_slope_threshold = 1.0;
    double linear_chi2_per_dof = LinearChi2PerDOF(outlier_filter_results.begin(), outlier_filter_results.end(), linear_intercept, linear_slope, linear_sigma);
    if(linear_chi2_per_dof < linear_chi2_per_dof_threshold and linear_sigma < linear_sigma_threshold and linear_slope < linear_slope_threshold) {
        // a linear fit is good enough
        // Let's subtract off the linear component and use the mode of result as the baseline

        return;
    }

    // Flat and linear are not good enough, let's try an exponential

    // now let's put the result through the exponential fitter
    double a = 0;
    double b = 0;
    double c = 0;
    FitExponential(outlier_filter_results, a, b, c);

    if (a == 0 or b == 0 or c == 0) {
        // oops! exponential fit didnt go right
        // let's just use the mode of the outlier filter for the baselines
        std::sort(outlier_filter_results.begin(), outlier_filter_results.end());
        double baseline_mode_val = robust_stats::Mode(outlier_filter_results.begin(), outlier_filter_results.end());
        double baseline_std = robust_stats::MedianAbsoluteDeviation(outlier_filter_results.begin(), outlier_filter_results.end(), baseline_mode_val);
        // now let's save it to our BaselineEstimate object baseline
        baseline.baseline = baseline_mode_val;
        baseline.stddev = baseline_std;
        baseline.target_num_frames = 0;
        baseline.num_frames = 0;
        baseline.num_samples = 0;
    } else {
        // now let's get the mode of our exponential fit!!
        std::vector<double> exp_result(samples.size());
        double current_time;
        for (size_t exp_it = 0; exp_it < samples.size(); ++exp_it){
            current_time = exp_it * 2;
            exp_result[exp_it] = a + b * std::exp(c * current_time);
        }
        std::sort(exp_result.begin(), exp_result.end());
        double baseline_mode_val = robust_stats::Mode(exp_result.begin(), exp_result.end());
        double baseline_std = robust_stats::MedianAbsoluteDeviation(exp_result.begin(), exp_result.end(), baseline_mode_val);
        // now let's save it to our BaselineEstimate object baseline
        baseline.baseline = baseline_mode_val;
        baseline.stddev = baseline_std;
        baseline.target_num_frames = 0;
        baseline.num_frames = 0;
        baseline.num_samples = 0;
    }
}

void darcy_baselines::DAQ(I3FramePtr frame) {

    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // let's read in our waveform
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");

    // I3Map to store pmt key and baselines
    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate>> Baselines = boost::make_shared<I3Map<CCMPMTKey, BaselineEstimate>>();

    // loop over each channel in waveforms
    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey key = p.first;
        uint32_t channel = p.second;
        CCMWaveformUInt16 const & waveform = waveforms->at(channel);
        std::vector<short unsigned int> const & samples = waveform.GetWaveform();

        // let's initalize our baseline object
        BaselineEstimate baseline;

        // let's get our baseline
        ProcessWaveform(samples, baseline);

        // let's save our baseline
        Baselines->insert({key, baseline});
    }

    frame->Put("BaselineEstimates", Baselines);
    PushFrame(frame);
}


