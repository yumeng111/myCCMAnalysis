#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <string>
#include <iostream>
#include <algorithm>

#include <icetray/I3Frame.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/robust_statistics.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/geometry/CCMGeometry.h>


// Identifies waveform regions that may correspond to potential peaks
// Collects statistics about these regions
/*
class BCMInfo(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)

    def Configure(self):
        self.output = open("bcm.txt", "a")
        pass


    def check_beam_trigger(self,samplevec):
         beamtrig_waveform = samplevec[self.beamtrig_channel].waveform
         return findFirstNIMSample(beamtrig_waveform)

    def DAQ(self, frame):
        samplevec = list(frame['CCMWaveforms'])
        t0 = self.check_beam_trigger(samplevec)
        if t0 == -1:
            return False

        bcm_start, bcm_peak = self.get_bcm_time(frame)

        frame["BCMStartTime"] = dataclasses.I3UInt64(int(bcm_start))
        frame["BCMPeakTime"] = dataclasses.I3UInt64(int(bcm_peak))

        return True

	def get_bcm_window(data, der, thresh = 0.3):
		i_max = np.argmax(data)
		d_max = data[i_max]
		d_mode = mode(data[i_max-1000:i_max])
		dev = median_absolute_deviation(data[i_max-1000:i_max])
		i_first = i_max - np.argmax((np.logical_or(data[:i_max+1] <= d_mode+dev, np.logical_and(np.abs(der[:i_max+1]) < thresh, data[:i_max+1] <= (d_max - d_mode)/2.0 + d_mode)))[::-1])
		i_last = i_max + np.argmax(data[i_max:] <= d_mode+dev)
		return (i_first, i_last)
*/

class BeamCurrentMonitorSummary : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string waveforms_name_;

    double time_before_peak_;
    double exp_smoothing_tau_;
    double derivative_threshold_;

    // Internal state
    bool geo_seen;
    CCMPMTKey bcm_key;
    size_t bcm_channel;

    // Constants
    constexpr static double ns_per_sample = 2.0;

public:
    BeamCurrentMonitorSummary(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void DAQ(I3FramePtr frame);
    void Finish();

    CCMBCMSummary GetBCMSummary(CCMWaveformUInt16 const & bcm_waveform);
    std::vector<double> SmoothWaveform(std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end);
    std::vector<double> ComputeDerviative(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end);
    size_t FindBaselineRegionLastPos(
            double baseline,
            double peak_value,
            size_t last_pos,
            double deriv_threshold,
            std::vector<double>::const_iterator smoothed_begin, std::vector<double>::const_iterator smoothed_end,
            std::vector<double>::const_iterator deriv_begin, std::vector<double>::const_iterator deriv_end);
    size_t FindBCMFirstPos(
            double baseline,
            double baseline_stddev,
            double peak_value,
            size_t last_pos,
            double deriv_threshold,
            std::vector<double>::const_iterator smoothed_begin, std::vector<double>::const_iterator smoothed_end,
            std::vector<double>::const_iterator deriv_begin, std::vector<double>::const_iterator deriv_end);
    size_t FindBCMLastPos(
            double baseline,
            double baseline_stddev,
            size_t first_pos,
            std::vector<double>::const_iterator smoothed_begin, std::vector<double>::const_iterator smoothed_end);
};

I3_MODULE(BeamCurrentMonitorSummary);

BeamCurrentMonitorSummary::BeamCurrentMonitorSummary(const I3Context& context) : I3Module(context),
    geometry_name_(""), waveforms_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("TimeBeforePeak", "Time in ns before the BCM peak to consider when computing the baseline and looking for the BCM start time", double(2.0));
    AddParameter("ExpSmoothingTau", "Time constant for exponential smoothing", double(10.0));
    AddParameter("DerivativeThreshold", "Theshold below which derivativ is considered to be zero", double(0.3));
}

void BeamCurrentMonitorSummary::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("TimeBeforePeak", time_before_peak_);
    GetParameter("DerivativeThreshold", derivative_threshold_);
}

std::vector<double> BeamCurrentMonitorSummary::SmoothWaveform(std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end) {

    size_t N = std::distance(begin, end);

    std::vector<double> smoothed_wf(N);

    // Average the first N samples to get the starting value for the exponential smoothing
    int init_avg_end_idx = std::min(3 * std::max(int(exp_smoothing_tau_ / ns_per_sample), 1), int(N));
    double exp_start = 0.0;
    std::vector<uint16_t>::const_iterator x_it = begin;
    for(size_t i=0; i < init_avg_end_idx; ++i, ++x_it) {
        exp_start += -double(*x_it);
    }
    exp_start /= init_avg_end_idx;

    double alpha = exp(-ns_per_sample / exp_smoothing_tau_);
    double y_i = exp_start;
    x_it = begin;
    for(size_t i=0; i < smoothed_wf.size(); ++i, ++x_it) {
        double x = -double(*x_it);
        y_i += (1.0 - alpha) * (x - y_i);
        smoothed_wf[i] = y_i;
    }

    // Keep track of previous value before box smoothing
    double unsmoothed_prev = smoothed_wf[0];

    // Box smoothing for the first index
    smoothed_wf[0] = (smoothed_wf[0] + smoothed_wf[1]) / 2.0;

    // Box smoothing for the last index
    double smoothed_wf_last_value = (smoothed_wf[smoothed_wf.size() - 2] + smoothed_wf[smoothed_wf.size() - 1]) / 2.0;

    // Box smoothing for everything in between
    for(size_t i=1; i < smoothed_wf.size() - 1; ++i) {
        double current = smoothed_wf[i];
        double next = smoothed_wf[i + 1];
        smoothed_wf[i] = (unsmoothed_prev + current + next) / 3.0;
        unsmoothed_prev = current;
    }
    smoothed_wf[smoothed_wf.size() - 1] = smoothed_wf_last_value;

    return smoothed_wf;
}

std::vector<double> BeamCurrentMonitorSummary::ComputeDerviative(std::vector<double>::const_iterator begin, std::vector<double>::const_iterator end) {

    size_t N = std::distance(begin, end);
    std::vector<double> derivative(N);
    double prevprev, prev, current, next, nextnext;

    // Use the forward difference approximation for the first point
    current = *begin;
    next = *(begin + 1);
    nextnext = *(begin + 2);
    derivative[0] = (-3.0 * current + 4.0 * next - nextnext) / (2.0 * ns_per_sample);

    // Use the central difference approximation for interior points
    for(size_t i=1; i<N-1; ++i) {
        prev = current;
        current = next;
        next = *(begin + i + 1);
        derivative[i] = (next - prev) / (2.0 * ns_per_sample);
    }

    // Use the backward difference approximation for the last point
    prevprev = prev;
    prev = current;
    current = next;
    derivative[N-1] = (3.0 * current - 4.0 * prev + prevprev) / (2.0 * ns_per_sample);

    return derivative;
}

void BeamCurrentMonitorSummary::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    geo_seen = true;
    bool found_bcm_key = false;
    for(auto const & key_om : geo.pmt_geo) {
        if(key_om.second.omtype == CCMOMGeo::OMType::BeamCurrentMonitor) {
            bcm_key = key_om.first;
            found_bcm_key = true;
        }
    }
    if(not found_bcm_key) {
        log_fatal("CCMGeometry does not contain a channel corresponding to a BeamCurrentMonitor");
    }
    bcm_channel = geo.pmt_channel_map.at(bcm_key);
}

size_t BeamCurrentMonitorSummary::FindBCMFirstPos(
        double baseline,
        double baseline_stddev,
        double peak_value,
        size_t last_pos,
        double deriv_threshold,
        std::vector<double>::const_iterator smoothed_begin, std::vector<double>::const_iterator smoothed_end,
        std::vector<double>::const_iterator deriv_begin, std::vector<double>::const_iterator deriv_end) {
    std::vector<double>::const_iterator smoothed_rbegin = smoothed_end - 1;
    std::vector<double>::const_iterator deriv_rbegin = deriv_end - 1;
    std::vector<double>::const_iterator smoothed_rend = smoothed_begin - 1;
    std::vector<double>::const_iterator deriv_rend = deriv_begin - 1;

    double min_val = baseline + baseline_stddev;
    double max_val_for_threshold = (peak_value - baseline) / 2.0 + baseline;
    size_t pos = last_pos;
    while(smoothed_rbegin != smoothed_rend and deriv_rbegin != deriv_rend) {
        double val = *smoothed_rbegin;
        double deriv = *deriv_rbegin;

        if(val <= min_val or (deriv <= deriv_threshold and val <= max_val_for_threshold)) {
            return pos;
        }

        if(pos == 0)
            break;
        --smoothed_rbegin;
        --deriv_rbegin;
        --pos;
    }
    return 0;
}

size_t BeamCurrentMonitorSummary::FindBCMLastPos(
        double baseline,
        double baseline_stddev,
        size_t first_pos,
        std::vector<double>::const_iterator smoothed_begin, std::vector<double>::const_iterator smoothed_end) {
    double min_val = baseline + baseline_stddev;
    size_t pos = first_pos;
    while(smoothed_begin != smoothed_end) {
        double val = *smoothed_begin;

        if(val <= min_val) {
            return pos;
        }

        if(pos == 0)
            break;
        --smoothed_begin;
        --pos;
    }
    return 0;
}

size_t BeamCurrentMonitorSummary::FindBaselineRegionLastPos(
        double baseline,
        double peak_value,
        size_t last_pos,
        double deriv_threshold,
        std::vector<double>::const_iterator smoothed_begin, std::vector<double>::const_iterator smoothed_end,
        std::vector<double>::const_iterator deriv_begin, std::vector<double>::const_iterator deriv_end) {
    std::vector<double>::const_iterator smoothed_rbegin = smoothed_end - 1;
    std::vector<double>::const_iterator deriv_rbegin = deriv_end - 1;
    std::vector<double>::const_iterator smoothed_rend = smoothed_begin - 1;
    std::vector<double>::const_iterator deriv_rend = deriv_begin - 1;

    double max_val_for_threshold = (peak_value - baseline) / 2.0 + baseline;
    size_t pos = last_pos;
    while(smoothed_rbegin != smoothed_rend and deriv_rbegin != deriv_rend) {
        double val = *smoothed_rbegin;
        double deriv = *deriv_rbegin;

        if(deriv <= deriv_threshold and val <= max_val_for_threshold) {
            return pos;
        }

        if(pos == 0)
            break;
        --smoothed_rbegin;
        --deriv_rbegin;
        --pos;
    }
    return 0;
}

CCMBCMSummary BeamCurrentMonitorSummary::GetBCMSummary(CCMWaveformUInt16 const & bcm_waveform) {
    std::vector<uint16_t> const & wf = bcm_waveform.GetWaveform();
    std::vector<uint16_t>::const_iterator peak_elem = std::min_element(wf.begin(), wf.end());

    // Find the peak
    int peak_pos = std::distance(wf.begin(), peak_elem);
    double peak_value = wf[peak_pos];

    // Define the search region to start a fixed time before the peak
    int first_search_begin_idx = std::max(peak_pos - int(time_before_peak_ / ns_per_sample), 0);
    // Define the search region to end at the peak
    int first_search_end_idx = std::min(peak_pos + 1, int(wf.size()));
    int first_search_length = first_search_end_idx - first_search_begin_idx;

    int last_search_begin_idx = peak_pos;
    int last_search_end_idx = int(wf.size());

    // Smooth the waveform in the search region
    // Applies box_smooth(exponential_smooth(-waveform))
    std::vector<double> smoothed_wf = SmoothWaveform(wf.begin() + first_search_begin_idx, wf.begin() + last_search_end_idx);

    // Compute the derivative in the search region
    std::vector<double> derivative = ComputeDerviative(smoothed_wf.begin(), smoothed_wf.end());

    // Get an initial estimate for the baseline
    double baseline = robust_stats::Mode(smoothed_wf.begin(), smoothed_wf.begin() + first_search_length);

    // Try to get away from the region with the BCM waveform
    size_t baseline_last_idx = FindBaselineRegionLastPos(
            baseline,
            peak_value,
            peak_pos,
            derivative_threshold_,
            smoothed_wf.begin(), smoothed_wf.begin() + first_search_length,
            derivative.begin(), derivative.begin() + first_search_length);

    // Re-estimate the baseline
    baseline = robust_stats::Mode(smoothed_wf.begin(), smoothed_wf.begin() + (baseline_last_idx - first_search_begin_idx + 1));

    // Estimate the stddev of the baseline
    double baseline_stddev = robust_stats::MedianAbsoluteDeviation(smoothed_wf.begin(), smoothed_wf.begin() + (baseline_last_idx - first_search_begin_idx + 1), baseline);

    size_t bcm_first_pos = FindBCMFirstPos(
            baseline,
            baseline_stddev,
            peak_value,
            peak_pos,
            derivative_threshold_,
            smoothed_wf.begin(), smoothed_wf.begin() + first_search_length,
            derivative.begin(), derivative.begin() + first_search_length);

    size_t bcm_last_pos = FindBCMLastPos(
            baseline,
            baseline_stddev,
            peak_pos,
            smoothed_wf.begin() + (last_search_begin_idx - first_search_begin_idx), smoothed_wf.begin() + (last_search_begin_idx - first_search_begin_idx) + (last_search_end_idx - last_search_begin_idx));

    CCMBCMSummary bcm;
    bcm.bcm_start_time = bcm_first_pos * ns_per_sample;
    bcm.bcm_end_time = bcm_last_pos * ns_per_sample;
    bcm.bcm_peak_time = peak_pos * ns_per_sample;
    bcm.bcm_peak_value = peak_value - baseline;
    bcm.bcm_integral = 0;
    bcm.bcm_baseline = baseline;
    bcm.bcm_baseline_stddev = baseline_stddev;

    for(size_t pos=bcm_first_pos; pos <= bcm_last_pos; ++pos) {
        bcm.bcm_integral += (double(wf[pos]) - baseline);
    }
    bcm.bcm_integral *= ns_per_sample;

    return bcm;
}

void BeamCurrentMonitorSummary::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("Geometry not seen yet!");
    }
    I3Vector<CCMWaveformUInt16> const & waveforms = frame->Get<I3Vector<CCMWaveformUInt16> const>(waveforms_name_);
    CCMWaveformUInt16 const & bcm_waveform = waveforms.at(bcm_channel);

    boost::shared_ptr<CCMBCMSummary> bcm = boost::make_shared<CCMBCMSummary>(GetBCMSummary(bcm_waveform));

    frame->Put("BCMSummary", bcm);
}

void BeamCurrentMonitorSummary::Finish() {
}
