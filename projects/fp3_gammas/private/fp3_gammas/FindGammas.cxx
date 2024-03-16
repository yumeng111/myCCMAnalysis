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
#include "daqtools/WaveformSmoother.h"

#include "daqtools/OnlineRobustStats.h"

// Start an anonymous namespace to keep the functions local to this file
namespace {

struct Extreme {
    size_t position;
    double value;
    double derivative;
    double second_derivative;
    double local_average;
};

size_t FindNeutronPeakIndex(std::vector<uint16_t> const & raw_waveform) {
    return std::distance(raw_waveform.begin(), std::max_element(raw_waveform.begin(), raw_waveform.end()));
}

double EstimateBaseline(WaveformSmootherDerivative & smoother, size_t samples_for_baseline) {
    smoother.Reset();
    size_t N = smoother.Size();
    N = std::min(N, samples_for_baseline);
    std::vector<double> baseline_samples;
    baseline_samples.reserve(N);
    for(size_t i=0; i<N; ++i) {
        baseline_samples.push_back(smoother.Value());
        smoother.Next();
    }
    std::sort(baseline_samples.begin(), baseline_samples.end());
    double baseline = robust_stats::Mode(baseline_samples.begin(), baseline_samples.end());
    return baseline;
}

double ComputeSecondDerivative(WaveformSmootherDerivative & smoother, size_t position, double bin_width) {
    size_t index = smoother.CurrentIndex();
    double result = 0;
    smoother.Reset(position+1);
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoothed_its = smoother.GetFullSmoothedWaveform();
    if(position == 0) {
        std::vector<double> data(smoothed_its.first, smoothed_its.first + 4);
        result = (2.0*data[0] - 5.0*data[1] + 4.0*data[2] - data[3]) / (bin_width * bin_width);
    } else if(position == smoother.Size() - 1) {
        std::vector<double> data(smoothed_its.second - 4, smoothed_its.second);
        result = (2.0*data[3] - 5.0*data[2] + 4.0*data[1] - data[0]) / (bin_width * bin_width);
    } else {
        std::vector<double> data(smoothed_its.first + position - 1, smoothed_its.first + position + 2);
        result = (data[0] - 2.0*data[1] + data[2]) / (bin_width * bin_width);
    }
    smoother.Reset(index);
    return result;
}

std::vector<Extreme> ComputeExtrema(WaveformSmootherDerivative & smoother, double baseline, size_t neutron_peak_index, double bin_width) {
    size_t index = smoother.CurrentIndex();
    smoother.Reset();
    std::vector<Extreme> extrema;
    size_t N = smoother.Size();
    double derivative = smoother.Derivative();
    double value = smoother.Value() - baseline;
    for(size_t i=0; i<N-1; ++i) {
        smoother.Next();
        double next_derivative = smoother.Derivative();
        if(derivative == 0 or (derivative > 0 and next_derivative < 0) or (derivative < 0 and next_derivative > 0)) {
            double second_derivative = ComputeSecondDerivative(smoother, i, bin_width);
            std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoother_its = smoother.GetFullSmoothedWaveform();
            std::vector<double> local_average_samples(smoother_its.first + (std::max(i, size_t(2)) - 2), smoother_its.first + std::min(i+2, N));
            double local_average = std::accumulate(local_average_samples.begin(), local_average_samples.end(), 0.0) / local_average_samples.size();
            extrema.push_back({i, value, derivative, second_derivative, local_average});
        }
        value = smoother.Value() - baseline;
        derivative = next_derivative;
    }
    smoother.Reset(index);
    return extrema;
}

void SortExtrema(std::vector<Extreme> & extrema) {
    std::sort(extrema.begin(), extrema.end(), [](Extreme const & a, Extreme const & b) {
        return a.value > b.value; // Sort in descending order
    });
}

std::vector<Extreme> SelectNoiseExtrema(std::vector<Extreme> const & extrema, size_t neutron_peak_index, double threshold) {
    std::vector<Extreme> noise_extrema;
    for(Extreme const & e : extrema) {
        if(e.position < neutron_peak_index and e.value < threshold) {
            noise_extrema.push_back(e);
        }
    }
    return noise_extrema;
}

double AverageNoise(std::vector<Extreme> const & noise_extrema) {
    std::vector<double> noise_values;
    noise_values.reserve(noise_extrema.size());
    for(Extreme const & e : noise_extrema) {
        noise_values.push_back(e.value);
    }
    return std::accumulate(noise_values.begin(), noise_values.end(), 0.0) / noise_values.size();
}

double MaxNoise(std::vector<Extreme> const & noise_extrema) {
    double max_noise = 0;
    for(Extreme const & e : noise_extrema) {
        max_noise = std::max(max_noise, e.value);
    }
    return max_noise;
}

size_t FindStartIndex(WaveformSmootherDerivative & smoother, size_t peak_index, double noise_threshold, double baseline) {
    size_t index = smoother.CurrentIndex();
    smoother.Reset(peak_index);
    size_t N = std::max(size_t(1), std::min(smoother.Size(), peak_index + 1)) - 1;
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoothed_its = smoother.GetFullSmoothedWaveform();
    size_t start_index = 0;
    for(std::vector<double>::const_iterator it=smoothed_its.first+N; it!=smoothed_its.first; --it) {
        double value = *it - baseline;
        if(value < noise_threshold) {
            start_index = std::distance(smoothed_its.first, it) + 1;
            break;
        }
    }
    smoother.Reset(index);
    return start_index;
}

size_t FindEndIndex(WaveformSmootherDerivative & smoother, size_t peak_index, double noise_threshold, double baseline) {
    size_t index = smoother.CurrentIndex();
    smoother.Reset(peak_index);
    size_t N = smoother.Size() - std::min(smoother.Size(), peak_index + 1);
    size_t end_index = smoother.Size() - 1;
    for(size_t i=0; i<N; ++i) {
        double value = smoother.Value() - baseline;
        if(value < noise_threshold) {
            end_index = i + peak_index - 1;
            break;
        }
        smoother.Next();
    }
    smoother.Reset(index);
    return end_index;
}

double SumBetweenIndicesInclusive(WaveformSmootherDerivative & smoother, size_t start_index, size_t end_index, double baseline) {
    size_t index = smoother.CurrentIndex();
    smoother.Reset(end_index);
    size_t N = std::min(smoother.Size(), end_index + 1);
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoothed_its = smoother.GetFullSmoothedWaveform();
    double sum = std::accumulate(smoothed_its.first + start_index, smoothed_its.first + N, 0.0);
    smoother.Reset(index);
    return sum;
}

std::vector<Extreme> SelectGammaExtrema(std::vector<Extreme> const & extrema, size_t neutron_start_index, double gamma_threshold, double local_average_threshold) {
    std::vector<Extreme> gamma_extrema;
    for(Extreme const & e : extrema) {
        if(e.position < neutron_start_index and e.value > gamma_threshold and e.local_average > local_average_threshold and e.second_derivative < 0) {
            gamma_extrema.push_back(e);
        }
    }
    return gamma_extrema;
}

} // namespace

class FindGammas : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;
    std::string nim_name_;

    // Internal state
    bool geo_seen;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;

    CCMPMTKey fp3_key_ = CCMPMTKey(9, 1, 0);

    size_t samples_for_baseline;

    double smoother_tau;
    double smoother_delta_t;
public:
    FindGammas(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
};

I3_MODULE(FindGammas);

FindGammas::FindGammas(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_(""), nim_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("NumSamplesForBaseline", "Number of samples to use for initial baseline estimate", size_t(400));
    AddParameter("SmoothingTau", "Time constant for waveform smoothing in ns", double(0.0));
    AddParameter("SmoothingDeltaT", "Bin width in ns for smoothing", double(2.0));
}

void FindGammas::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("NumSamplesForBaseline", samples_for_baseline);
    GetParameter("SmoothingTau", smoother_tau);
    GetParameter("SmoothingDeltaT", smoother_delta_t);
}

void FindGammas::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);

    // Cache the trigger channel map
    pmt_channel_map_ = geo.pmt_channel_map;

    geo_seen = true;
    PushFrame(frame);
}

void FindGammas::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    if(not frame->Has(waveforms_name_)) {
        log_fatal("Could not find CCMWaveforms object with the key named \"%s\" in the DAQ frame.", waveforms_name_.c_str());
    }
    PushFrame(frame);
}
