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
#include <dataclasses/physics/CCMBCMSummary.h>
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

struct GammaEntry {
    double baseline;
    double average_noise_level;
    double max_noise_level;
    size_t num_noise_peaks;
    size_t waveform_length;
    double gamma_start_time;
    double gamma_end_time;
    double gamma_charge_sum;
    double gamma_height;
    double gamma_derivative;
    double gamma_second_derivative;
    double gamma_local_average;
    size_t gamma_peak_time;
    double neutron_height;
    double neutron_derivative;
    double neutron_second_derivative;
    double neutron_local_average;
    size_t neutron_peak_time;
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
    double result = 0;
    smoother.ComputeTo(position+1);
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
    return result;
}

double ComputeLocalAverage(WaveformSmootherDerivative & smoother, size_t position, size_t samples_for_local_average) {
    size_t N = smoother.Size();
    size_t N_before = std::min(position, size_t(samples_for_local_average / 2));
    size_t N_after = std::min(N - position, samples_for_local_average - N_before);
    smoother.ComputeTo(position + N_after - 1);
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoothed_its = smoother.GetFullSmoothedWaveform();
    size_t start = std::max(int(0), int(position) - int(N_before));
    size_t end = std::min(N, position + N_after);
    std::vector<double> local_average_samples(smoothed_its.first + start, smoothed_its.first + end);
    double local_average = std::accumulate(local_average_samples.begin(), local_average_samples.end(), 0.0) / local_average_samples.size();
    return local_average;
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
            double local_average = ComputeLocalAverage(smoother, i, 5);
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
    smoother.ComputeTo(peak_index);
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
    smoother.ComputeTo(end_index);
    size_t N = std::min(smoother.Size(), end_index + 1);
    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> smoothed_its = smoother.GetFullSmoothedWaveform();
    double sum = std::accumulate(smoothed_its.first + start_index, smoothed_its.first + N, 0.0);
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

double ComputeTimeOffset(bool & valid_time_offset, NIMLogicPulseSeriesMap const & nim_pulses, CCMTriggerKey const & bcm_board_key, CCMTriggerKey const & fp3_board_key, CCMBCMSummary const & bcm_summary) {
    if(not nim_pulses.count(bcm_board_key)) {
        valid_time_offset = false;
        return 0;
    }
    if(not nim_pulses.at(bcm_board_key).size()) {
        valid_time_offset = false;
        return 0;
    }
    double bcm_board_time = nim_pulses.at(bcm_board_key).front().GetNIMPulseTime();

    if(not nim_pulses.count(fp3_board_key)) {
        valid_time_offset = false;
        return 0;
    }
    if(not nim_pulses.at(fp3_board_key).size()) {
        valid_time_offset = false;
        return 0;
    }
    double fp3_board_time = nim_pulses.at(fp3_board_key).front().GetNIMPulseTime();

    double time_offset = bcm_summary.bcm_start_time - (bcm_board_time - fp3_board_time);

    valid_time_offset = true;
    return time_offset;
}

std::vector<GammaEntry> ComputeGammaEntries(std::vector<Extreme> const & gamma_extrema, WaveformSmootherDerivative & smoother, double baseline, double noise_average, double noise_max, size_t neutron_peak_index, size_t neutron_start_index, size_t num_noise_peaks, double bin_width) {
    double neutron_height = smoother.Value(neutron_peak_index) - baseline;
    double neutron_derivative = smoother.Derivative(neutron_peak_index);
    double neutron_second_derivative = ComputeSecondDerivative(smoother, neutron_peak_index, bin_width);
    double neutron_local_average = ComputeLocalAverage(smoother, neutron_peak_index, 5);
    double neutron_time = neutron_peak_index * bin_width;

    std::vector<GammaEntry> gamma_entries;
    for(size_t i=0; i<gamma_extrema.size(); ++i) {
        Extreme const & gamma = gamma_extrema[i];
        size_t start_index = FindStartIndex(smoother, gamma.position, noise_max, baseline);
        size_t end_index = FindEndIndex(smoother, gamma.position, noise_max, baseline);
        double start_time = start_index * bin_width;
        double end_time = end_index * bin_width;
        double peak_time = gamma.position * bin_width;
        double integral = SumBetweenIndicesInclusive(smoother, start_index, end_index, baseline) * bin_width;
        double height = gamma.value - baseline;
        double derivative = gamma.derivative;
        double second_derivative = gamma.second_derivative;
        double local_average = gamma.local_average;

        gamma_entries.push_back({baseline, noise_average, noise_max, num_noise_peaks, smoother.Size(), start_time, end_time, integral, height, derivative, second_derivative, local_average, peak_time, neutron_height, neutron_derivative, neutron_second_derivative, neutron_local_average, neutron_time});

    }
    return gamma_entries;
}


} // namespace

class FindGammas : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;
    std::string nim_pulses_name_;
    std::string bcm_summary_name_;

    // Internal state
    bool geo_seen;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;

    CCMPMTKey fp3_key_;
    CCMPMTKey bcm_key_;
    size_t fp3_channel_;

    CCMTriggerKey bcm_board_key;
    CCMTriggerKey fp3_board_key;

    size_t samples_for_baseline_;
    size_t samples_for_local_average_;
    double noise_maximum_;
    double gamma_threshold_;

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
    geometry_name_(""), waveforms_name_(""), nim_pulses_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("NIMPulsesName", "Key to input NIMLogicPulseSeriesMap", std::string("NIMPulses"));
    AddParameter("BCMSummaryName", "Key to input CCMBCMSummary", std::string("CCMBCMSummary"));
    AddParameter("FlightPath3Key", "CCMPMTKey for the flight path 3", CCMPMTKey(9, 1, 0));
    AddParameter("BCMKey", "CCMPMTKey for the BCM", CCMPMTKey(10, 1, 0));
    AddParameter("NumSamplesForBaseline", "Number of samples to use for initial baseline estimate", size_t(1000));
    AddParameter("NumSamplesForLocalAverage", "Number of samples to use for local average", size_t(50));
    AddParameter("NoiseMaximum", "Maximum for noise", double(25.0));
    AddParameter("GammaThreshold", "Threshold for gamma peaks", double(80.0));
    AddParameter("SmoothingTau", "Time constant for waveform smoothing in ns", double(0.0));
    AddParameter("SmoothingDeltaT", "Bin width in ns for smoothing", double(2.0));
}

void FindGammas::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
    GetParameter("BCMSummaryName", bcm_summary_name_);
    GetParameter("FlightPath3Key", fp3_key_);
    GetParameter("BCMKey", bcm_board_key);
    GetParameter("NumSamplesForBaseline", samples_for_baseline_);
    GetParameter("NumSamplesForLocalAverage", samples_for_local_average_);
    GetParameter("NoiseMaximum", noise_maximum_);
    GetParameter("GammaThreshold", gamma_threshold_);
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

    // Find the flight path 3 channel
    I3Map<CCMPMTKey, uint32_t>::const_iterator it = pmt_channel_map_.find(fp3_key_);
    if(it == pmt_channel_map_.end()) {
        log_fatal("Could not find the flight path 3 key in the geometry.");
    }
    fp3_channel_ = it->second;

    bcm_board_key = geo.trigger_copy_map.at(bcm_key_);
    fp3_board_key = geo.trigger_copy_map.at(fp3_key_);

    geo_seen = true;
    PushFrame(frame);
}

void FindGammas::DAQ(I3FramePtr frame) {
    // Check if we have the prerequisite information
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    if(not frame->Has(waveforms_name_)) {
        log_fatal("Could not find CCMWaveforms object with the key named \"%s\" in the DAQ frame.", waveforms_name_.c_str());
    }
    if(not frame->Has(nim_pulses_name_)) {
        log_fatal("Could not find NIMLogicPulseSeriesMap object with the key named \"%s\" in the DAQ frame.", nim_pulses_name_.c_str());
    }
    if(not frame->Has("CCMBCMSummary")) {
        log_fatal("Could not find CCMBCMSummary object with the key named \"CCMBCMSummary\" in the DAQ frame.");
    }

    // Get what we need to compute the time offset
    I3Map<CCMTriggerKey, NIMLogicPulseSeries> const & nim_pulses = frame->Get<I3Map<CCMTriggerKey, NIMLogicPulseSeries> const>(nim_pulses_name_);
    CCMBCMSummary const & bcm_summary = frame->Get<CCMBCMSummary const>(bcm_summary_name_);

    bool valid_time_offset = false;
    double time_offset = ComputeTimeOffset(valid_time_offset, nim_pulses, bcm_board_key, fp3_board_key, bcm_summary);
    if(not valid_time_offset) {
        log_warn("Could not compute a valid time offset. Skipping this frame.");
        PushFrame(frame);
        return;
    }

    // Get the waveforms
    boost::shared_ptr<I3Vector<CCMWaveformUInt16> const> waveforms = frame->Get<boost::shared_ptr<I3Vector<CCMWaveformUInt16>const>>(waveforms_name_);
    if(waveforms->empty()) {
        log_fatal("No waveforms found in the CCMWaveforms object.");
    }

    // Get the waveform for flight path 3
    std::vector<uint16_t> const & waveform = waveforms->at(fp3_channel_).GetWaveform();

    // Create a smoother that we can use for all the calculations
    WaveformSmootherDerivative smoother(waveform.begin(), waveform.end(), smoother_tau, smoother_delta_t);

    // Find the neutron peak
    size_t neutron_peak_index = FindNeutronPeakIndex(waveform);

    // Estimate the baseline
    double baseline = EstimateBaseline(smoother, samples_for_baseline_);

    // Compute the extrema
    std::vector<Extreme> extrema = ComputeExtrema(smoother, baseline, neutron_peak_index, smoother_delta_t);
    SortExtrema(extrema);

    // Evaluate the noise
    std::vector<Extreme> noise_extrema = SelectNoiseExtrema(extrema, neutron_peak_index, noise_maximum_);
    double noise_average = AverageNoise(noise_extrema);
    double noise_max = MaxNoise(noise_extrema);

    // Find the start index for the neutron peak
    size_t neutron_start_index = FindStartIndex(smoother, neutron_peak_index, noise_max, baseline);

    // Find the gamma extrema
    std::vector<Extreme> gamma_extrema = SelectGammaExtrema(extrema, neutron_start_index, gamma_threshold_, noise_average);

    // Compute the gamma entries
    std::vector<GammaEntry> gamma_entries = ComputeGammaEntries(gamma_extrema, smoother, baseline, noise_average, noise_max, neutron_peak_index, neutron_start_index, noise_extrema.size(), smoother_delta_t);

    PushFrame(frame);
}
