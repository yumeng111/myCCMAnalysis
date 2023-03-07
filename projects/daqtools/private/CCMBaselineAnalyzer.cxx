#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Orientation.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"


namespace {

class SourceExhausted {};

template<class T>
class Source {
public:
    virtual T output() = 0;
};

template<class T>
class WaveformSource : public Source<T> {
protected:
    std::vector<T> const & source;
    size_t current_pos;
public:
    WaveformSource(CCMWaveform<T> const & source_wf) :
        source(source_wf.GetWaveform()), current_pos(0) {}
    inline T output() override {
        if(current_pos < source.size()) {
            T const & x = source[current_pos];
            current_pos += 1;
            return x;
        } else {
            throw SourceExhausted();
            return 0;
        }
    }
    size_t Size() const {
        return source.size();
    }
};

template<class T, class U>
class Filter : public Source<U> {
protected:
    Source<T> * input;
public:
    Filter() {}
    Filter(Source<T> * input) : input(input) {}
    virtual inline U filter(T const & t) = 0;
    inline U output() override {
        return filter(input->output());
    }
};

template<class T, class U>
class ExponentialFilter : public Filter<T, U> {
protected:
    U tau;
    U delta_t;
    U y_i;
    U alpha;
    bool init;
public:
    ExponentialFilter(U tau, U delta_t, Source<T> * input) :
        tau(tau), delta_t(delta_t), init(true), Filter<T, U>(input) {
        alpha = exp(-delta_t / tau);
    }

    inline U filter(T const & t) override {
        if(init) {
            y_i = t;
            init = false;
        }
        y_i += (1.0 - alpha) * (t - y_i);
        return y_i;
    }
};

template<class T, class U>
class MultiplicativeFilter : public Filter<T, U> {
protected:
    U c;
public:
    MultiplicativeFilter(U c, Source<T> * input) :
        c(c), Filter<T, U>(input) {}
    inline U filter(T const & t) override {
        U x = U(t)*c;
        return x;
    }
};

template<class T, class U>
class Buffer : public Filter<T, U> {
protected:
    std::deque<T> buffer;
    size_t max_size;
    bool buffer_filling;
    bool source_exhausted;
public:
    Buffer(std::deque<T> init_buffer, size_t max_size, Source<T> * source_input) :
        buffer(init_buffer), max_size(max_size), buffer_filling(true), source_exhausted(false), Filter<T, U>(source_input) {
        if(buffer.size() >= max_size)
            buffer_filling = false;
        size_t to_pop = buffer.size() - std::min(max_size, buffer.size());
        for(size_t i=0; i<to_pop; ++i) {
            buffer.pop_front();
        }
    }
    inline U output() override {
        if(buffer_filling) {
            try {
                T t = Filter<T, U>::input->output();
                buffer.push_back(t);
                if(buffer.size() == max_size)
                    buffer_filling = false;
            } catch(SourceExhausted const & e) {
                buffer.pop_front();
                if(buffer.size() == 0)
                    throw SourceExhausted();
            }
        } else if(source_exhausted) {
            buffer.pop_front();
            if(buffer.size() == 0)
                throw SourceExhausted();
        } else {
            try {
                T t = Filter<T, U>::input->output();
                buffer.push_back(t);
                if(buffer.size() > max_size)
                    buffer.pop_front();
            } catch(SourceExhausted const & e) {
                buffer.pop_front();
                if(buffer.size() == 0)
                    throw SourceExhausted();
            }
        }
        return filter(buffer);
    }
    virtual inline U filter(std::deque<T> const & buffer) = 0;
};

template<class T, class U>
class BufferView : public Filter<T, U> {
protected:
    std::deque<T> buffer;
    size_t buffer_begin;
    size_t buffer_end;
public:
    BufferView(Source<T> * input) : Filter<T, U>(input) {
        buffer_begin = 0;
        buffer_end = 0;
    }

    inline U filter(T const & t) override {
        buffer.push_back(t);
        buffer_end += 1;
        return t;
    }

    virtual inline std::deque<T> const & GetBuffer() const {
        return buffer;
    }

    virtual inline std::tuple<size_t, size_t> GetWindow() const {
        return {buffer_begin, buffer_end};
    }

    virtual inline void Reset() {
        if(buffer.size() > 0) {
            buffer_begin = buffer_end - 1;
            buffer = {buffer.back()};
        } else {
            buffer_begin = buffer_end;
        }
    }
};

template<class T, class U>
class SmoothingFilter : public Buffer<T, U> {
protected:
    size_t zero_pos;
    std::deque<U> coeff;
    U total;
public:
    SmoothingFilter(size_t zero_pos, std::deque<U> coeff, Source<T> * input) :
        zero_pos(zero_pos), coeff(coeff), total(0), Buffer<T, U>({}, coeff.size(), input) {
        for(size_t i=0; i<coeff.size(); ++i) {
            total += coeff[i];
        }
        size_t n_preload = zero_pos;
        for(size_t i=0; i<n_preload; ++i)
            Buffer<T, U>::buffer.push_back(input->output());
    }
    inline U filter(std::deque<T> const & buffer) override {
        U x;
        if(Buffer<T, U>::buffer_filling) {
            size_t n_excluded = coeff.size() - buffer.size();
            U this_total = 0;
            U sum = 0;
            for(size_t i=0; i<buffer.size(); ++i) {
                size_t j = i + n_excluded;
                this_total += coeff[j];
                sum += coeff[j] * buffer[i];
            }
            if(this_total == 0)
                x = 0;
            else
                x = sum / this_total;
        } else if(Buffer<T, U>::source_exhausted) {
            if(buffer.size() <= zero_pos)
                throw SourceExhausted();
            U this_total = 0;
            U sum = 0;
            for(size_t i=0; i<buffer.size(); ++i) {
                this_total += coeff[i];
                sum += coeff[i] * buffer[i];
            }
            if(this_total == 0)
                x = 0;
            else
                x = sum / this_total;
        } else {
            U sum = 0;
            for(size_t i=0; i<coeff.size(); ++i) {
                sum += coeff[i] * buffer[i];
            }
            x = sum / total;
        }
        return x;
    }
    inline U filter(T const & t) override {return t;};
};

template<class T, class U>
class BackDerivativeFilter : public Filter<T, U> {
protected:
    U bin_width;
    U t_last;
    bool init;
public:
    BackDerivativeFilter(U bin_width, Source<T> * input) :
        bin_width(bin_width), init(true), Filter<T, U>(input) {}
    inline U filter(T const & t) override {
        U x;
        if(init) {
            x = 0;
            init = false;
        } else {
            x = (t - t_last) / (bin_width);
        }
        t_last = t;
        return x;
    }
};

template<class T, class U>
class FrontDerivativeFilter : public Filter<T, U> {
protected:
    U bin_width;
    U t_last;
    bool ending;
    bool end;
public:
    FrontDerivativeFilter(U bin_width, Source<T> * input) :
        bin_width(bin_width), ending(false), end(false), Filter<T, U>(input) {
        t_last = Filter<T, U>::input->output();
    }
    inline U filter(T const & t) override {
        U x;
        if(end) {
            throw SourceExhausted();
        } else if(ending) {
            x = 0;
            end = true;
        } else {
            x = (t - t_last) / (bin_width);
        }
        t_last = t;
        return x;
    }
};

template<class T, class U, class K>
struct DerivativeWindowCapture {
    Source<K> & d_buffer;
    BufferView<T, U> & raw_buffer;
    K max_derivative;
    int state;
    size_t pos;
    // States:
    // 0 - outside of region without activity
    // 1 - inside of region without activity

    DerivativeWindowCapture(Source<K> & d_buffer, BufferView<T, U> & raw_buffer, K max_derivative) :
        d_buffer(d_buffer), raw_buffer(raw_buffer), max_derivative(max_derivative), state(0), pos(0) {}

    bool Inside(T const & t) const {
        return t <= max_derivative and t >= -max_derivative;
    }

    bool Next() {
        pos += 1;
        K derivative = d_buffer.output();
        bool inside = Inside(derivative);
        bool found_region = false;
        if(state == 0) { // currently outside region
            if(inside) { // moving inside region
                raw_buffer.Reset();
            } else { // still outside region
            }
        } else { // already inside
            if(not inside) { // moving outside region
                found_region = true;
            } else { // still inside region
            }
        }
        state = inside;
        return found_region;
    }

    std::tuple<typename std::deque<U>::const_iterator, typename std::deque<U>::const_iterator, size_t> GetBuffer() const {
        return {raw_buffer.GetBuffer().cbegin(), raw_buffer.GetBuffer().cend() - 1, raw_buffer.GetBuffer().size()};
    }

    std::tuple<size_t, size_t> GetWindow() const {
        return raw_buffer.GetWindow();
    }
};

struct WindowStats {
    double mean = 0;
    double variance = 0;
    double k_ = 0;
    double M_ = 0;
    double S_ = 0;

    long double time;

    template<typename T>
    void AddSample(T const & x) {
        k_ += 1;
        double delta = x - M_;
        M_ += delta / k_;
        double delta2 = x - M_;
        S_ += delta * delta2;
    }

    void Finalize() {
        mean = M_;
        if(k_ > 1)
            variance = S_ / (k_ - 1);
        else
            variance = 0;
    }

    template<typename T>
    void AddSamples(std::deque<T> const & buffer) {
        for(size_t i=0; i<buffer.size(); ++i) {
            AddSample<T>(buffer[i]);
        }
        Finalize();
    }

    template<typename Iterator>
    void AddSamples(Iterator begin, Iterator end) {
        using value_type = typename std::iterator_traits<Iterator>::value_type;
        while(begin != end) {
            AddSample<value_type>(*begin);
            ++begin;
        }
        Finalize();
    }

    WindowStats() :
        k_(0), M_(0), S_(0), mean(0), variance(0) {
    }

    template<typename Iterator>
    WindowStats(Iterator begin, Iterator end) :
        k_(0), M_(0), S_(0), mean(0), variance(0) {
        AddSamples<Iterator>(begin, end);
    }

    template<typename T>
    WindowStats(std::deque<T> const & buffer) :
        k_(0), M_(0), S_(0), mean(0), variance(0) {
        AddSamples<T>(buffer);
    }
};

} // anonymous namespace


// Identifies waveform regions that may correspond to regions without any activity
// Collects statistics about these regions for analyzing the PMT baselines

class CCMBaselineAnalyzer : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string daq_config_name_;
    std::string waveforms_name_;
    std::string abs_time_name_;

    std::string baseline_fit_output_name_;

    double initial_derivative_threshold_;
    size_t minimum_sample_length_;
    double baseline_minimum_window_fraction_;
    double baseline_maximum_window_fraction_;
    size_t baseline_sample_edge_cut_;
    size_t num_triggers_for_threshold_;
    size_t max_waveform_sample_;

    std::map<size_t, double> current_derivative_threshold_;
    std::map<size_t, size_t> current_window_size_;
    std::vector<I3FramePtr> cached_frames;
    bool thresholds_tuned_;

    size_t n_daq_frames = 0;

    std::tuple<WindowStats, std::vector<WindowStats>, size_t, std::vector<std::tuple<size_t, size_t>>> GetBaselineStats(CCMWaveformUInt16 const & wf, double derivative_threshold, size_t min_window_size);
    void FindThreshold(std::vector<I3FramePtr> const & frames, size_t wf_idx);
    void FindThresholds(std::vector<I3FramePtr> const & frames);
    void AddBaselineStats(I3FramePtr frame);

public:
    CCMBaselineAnalyzer(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMBaselineAnalyzer);

CCMBaselineAnalyzer::CCMBaselineAnalyzer(const I3Context& context) : I3Module(context),
    geometry_name_(""), daq_config_name_(""), waveforms_name_(""), abs_time_name_(""), baseline_fit_output_name_(""),
    initial_derivative_threshold_(0.3), minimum_sample_length_(30),
    baseline_minimum_window_fraction_(0.005), baseline_maximum_window_fraction_(0.9),
    baseline_sample_edge_cut_(3), num_triggers_for_threshold_(100),
    max_waveform_sample_(8000),
    current_derivative_threshold_(), current_window_size_(), cached_frames(),
    thresholds_tuned_(false) {

    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMAnalysis::Binary::CCMDAQConfig>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("InitialDerivativeThreshold", "The threshold first used to find regions of inactivity. This is iteratively increased until regions of inactivity represent a certain fraction of the waveform.",initial_derivative_threshold_);
    AddParameter("MinimumSampleLength", "The smallest number of consecutive samples to use when finding regions of inactivity.", minimum_sample_length_);
    AddParameter("BaselineMinimumWindowFraction", "The smallest allowed fraction that regions of inactivity may occupy across the sampled waveforms. This is needed to ensure we have enough samples to do the calculation.", baseline_minimum_window_fraction_);
    AddParameter("BaselineMaximumWindowFraction", "The largest allowed fraction that regions of inactivity may occupy across the sampled waveforms. This is only present to alert us to the unlikely scenario where a PMT has too much activity to get a measurement of the baseline.", baseline_maximum_window_fraction_);
    AddParameter("BaselineSampleEdgeCut", "The number of samples to cut from the edges of each region of inactivity. This helps us avoid contamination of the baseline measurement from the adjacent regions of activity.", baseline_sample_edge_cut_);
    AddParameter("NumTriggersForThreshold","The number of triggers to use when tuning the derivative threshold.", num_triggers_for_threshold_);
    AddParameter("BaselineFitOutputName", "The output key of the baseline fit.", std::string("BaselineFit"));
    AddParameter("AbsoluteTimeName", "Key for absoluting timing information of frame", std::string("FirstTriggerTime"));
    AddParameter("MaxWaveformSample", "The maximum sample number to consider when finding baseline samples.", max_waveform_sample_);
}

void CCMBaselineAnalyzer::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMDAQConfigName", daq_config_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("InitialDerivativeThreshold", initial_derivative_threshold_);
    GetParameter("MinimumSampleLength", minimum_sample_length_);
    GetParameter("BaselineMinimumWindowFraction", baseline_minimum_window_fraction_);
    GetParameter("BaselineMaximumWindowFraction", baseline_maximum_window_fraction_);
    GetParameter("BaselineSampleEdgeCut", baseline_sample_edge_cut_);
    GetParameter("NumTriggersForThreshold", num_triggers_for_threshold_);
    GetParameter("BaselineFitOutputName", baseline_fit_output_name_);
    GetParameter("AbsoluteTimeName", abs_time_name_);
    GetParameter("MaxWaveformSample", max_waveform_sample_);
}


std::tuple<WindowStats, std::vector<WindowStats>, size_t, std::vector<std::tuple<size_t, size_t>>> CCMBaselineAnalyzer::GetBaselineStats(CCMWaveformUInt16 const & wf, double derivative_threshold, size_t min_window_size) {
    std::vector<WindowStats> stats;
    std::vector<std::tuple<size_t, size_t>> positions;
    WindowStats total;

    if(wf.GetWaveform().size() == 0)
        return {total, stats, 0, positions};

    WaveformSource<uint16_t> wf_source(wf);
    MultiplicativeFilter<uint16_t, double> m_filter(-1.0, &wf_source);
    BufferView<double, double> buffer_m(&m_filter);
    ExponentialFilter<double, double> exp_filter(10.0, 2.0, &buffer_m);
    std::deque<double> smoothing_kernel = {1,1,1};
    SmoothingFilter<double, double> smoothing_filter(1, smoothing_kernel, &exp_filter);
    BackDerivativeFilter<double, double> derivative_filter(2.0, &smoothing_filter);
    DerivativeWindowCapture<double, double, double> window(derivative_filter, buffer_m, derivative_threshold);


    size_t all_samples = wf_source.Size();
    size_t inactive_samples = 0;

    try {
        size_t sample_number = 0;
        while(true) {
            bool have_window = window.Next();
            if(have_window or (sample_number == max_waveform_sample_ and window.state == 1)) {
                std::tuple<std::deque<double>::const_iterator, std::deque<double>::const_iterator, size_t> window_buffer = window.GetBuffer();
                std::tuple<size_t, size_t> window_bounds = window.GetWindow();
                std::deque<double>::const_iterator begin = std::get<0>(window_buffer);
                std::deque<double>::const_iterator end = std::get<1>(window_buffer);
                size_t size = std::get<2>(window_buffer);
                if(size >= min_window_size) {
                    if(2 * baseline_sample_edge_cut_ + size < min_window_size)
                        continue;
                    positions.push_back(window_bounds);
                    stats.emplace_back(begin + baseline_sample_edge_cut_, end - baseline_sample_edge_cut_);
                    stats.back().time = ((long double)(window.pos) - 1)*2.0 + wf.GetStartTime();
                    total.AddSamples(begin + baseline_sample_edge_cut_, end - baseline_sample_edge_cut_);
                    inactive_samples += size;
                }
            }
            if(sample_number == max_waveform_sample_) {
                break;
            }
            sample_number += 1;
        }
    } catch(SourceExhausted const & e) {
    }

    return {total, stats, inactive_samples, positions};
}

void CCMBaselineAnalyzer::FindThreshold(std::vector<I3FramePtr> const & frames, size_t wf_idx) {
    double threshold = initial_derivative_threshold_;
    size_t min_window_size = minimum_sample_length_;

    double delta = threshold / 2.0;

    double last_ratio = 0;
    size_t n_last_was_below = 0;
    size_t max_last_below = 10;

    while(true) {
        size_t all_samples = 0;
        size_t inactive_samples = 0;

        for(size_t i=0; i<frames.size(); ++i) {
            CCMWaveformUInt16Series const & waveforms = frames[i]->Get<CCMWaveformUInt16Series>(waveforms_name_);
            std::tuple<WindowStats, std::vector<WindowStats>, size_t, std::vector<std::tuple<size_t, size_t>>> window_stats =
                GetBaselineStats(waveforms[wf_idx], threshold, min_window_size);
            size_t wf_inactive_samples = std::get<2>(window_stats);
            inactive_samples += wf_inactive_samples;
            all_samples += waveforms[wf_idx].GetWaveform().size();
        }

        double inactive_ratio = double(inactive_samples) / double(all_samples);
        if(inactive_ratio >= baseline_minimum_window_fraction_ and inactive_ratio <= baseline_maximum_window_fraction_) {
            break;
        }
        if(inactive_ratio > baseline_maximum_window_fraction_) {
            threshold -= delta;
            delta /= 2.0;
            n_last_was_below = 0;
        } else if(inactive_ratio < baseline_minimum_window_fraction_) {
            if(inactive_ratio == last_ratio and n_last_was_below >= max_last_below) {
                size_t reduce_by = std::max(size_t(1), size_t(min_window_size * 0.1));
                if(min_window_size <= (std::max(baseline_sample_edge_cut_ * 2 + 1, size_t(10)) + reduce_by)) {
                    if(inactive_ratio == 0) {
                        threshold = initial_derivative_threshold_;
                        min_window_size = minimum_sample_length_;
                        break;
                    }
                    throw std::runtime_error("Could not find derivative threshold that meets criteria.");
                }
                threshold = initial_derivative_threshold_;
                delta = threshold / 2.0;
                last_ratio = 0;
                n_last_was_below = 0;
                min_window_size -= reduce_by;
                continue;
            }
            threshold += delta;
            n_last_was_below += 1;
        }
    }
    current_derivative_threshold_[wf_idx] = threshold;
    current_window_size_[wf_idx] = min_window_size;
}

void CCMBaselineAnalyzer::FindThresholds(std::vector<I3FramePtr> const & frames) {
    size_t n_waveforms = frames[0]->Get<CCMWaveformUInt16Series>(waveforms_name_).size();
    for(size_t i=0; i<n_waveforms; ++i) {
        std::cout << "Finding thresholds for waveform number " << i << "/" << n_waveforms << std::endl;
        FindThreshold(frames, i);
    }
}

void CCMBaselineAnalyzer::AddBaselineStats(I3FramePtr frame) {
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series>(waveforms_name_);
    long double frame_time = (long double)(frame->Get<I3PODHolder<int64_t>>(abs_time_name_).value) * 2.0;
    size_t size = waveforms.size();
    boost::shared_ptr<I3Vector<double>> total_means(new I3Vector<double>(size));
    boost::shared_ptr<I3Vector<double>> total_variances(new I3Vector<double>(size));
    boost::shared_ptr<I3Vector<I3Vector<double>>> window_means(new I3Vector<I3Vector<double>>(size));
    boost::shared_ptr<I3Vector<I3Vector<double>>> window_variances(new I3Vector<I3Vector<double>>(size));
    boost::shared_ptr<I3Vector<I3Vector<long double>>> sample_times(new I3Vector<I3Vector<long double>>(size));
    boost::shared_ptr<I3Vector<I3Vector<size_t>>> sample_begin(new I3Vector<I3Vector<size_t>>(size));
    boost::shared_ptr<I3Vector<I3Vector<size_t>>> sample_end(new I3Vector<I3Vector<size_t>>(size));

    for(size_t i=0; i<waveforms.size(); ++i) {
        std::tuple<WindowStats, std::vector<WindowStats>, size_t, std::vector<std::tuple<size_t, size_t>>> window_stats =
            GetBaselineStats(waveforms[i], current_derivative_threshold_[i], current_window_size_[i]);
        WindowStats const & total_stats = std::get<0>(window_stats);
        std::vector<WindowStats> stats = std::get<1>(window_stats);
        std::vector<std::tuple<size_t, size_t>> positions = std::get<3>(window_stats);
        I3Vector<size_t> begins;
        I3Vector<size_t> ends;
        double total_mean = total_stats.mean;
        double total_variance = total_stats.variance;
        I3Vector<double> window_mean;
        I3Vector<double> window_variance;
        I3Vector<long double> times;
        for(size_t j=0; j<stats.size(); ++j) {
            window_mean.push_back(stats[j].mean);
            window_variance.push_back(stats[j].variance);
            times.push_back(stats[j].time + frame_time);
            begins.push_back(std::get<0>(positions[j]));
            ends.push_back(std::get<1>(positions[j]));
        }
        total_means->operator[](i) = total_mean;
        total_variances->operator[](i) = total_variance;
        window_means->operator[](i) = window_mean;
        window_variances->operator[](i) = window_variance;
        sample_times->operator[](i) = times;
        sample_begin->operator[](i) = begins;
        sample_end->operator[](i) = ends;
    }

    frame->Put("BaselineEstimates", total_means);
    frame->Put("BaselineEstimateVariances", total_variances);
    frame->Put("BaselineSamples", window_means);
    frame->Put("BaselineSampleVariances", window_variances);
    frame->Put("BaselineSampleTimes", sample_times);
    frame->Put("BaselineSampleBeginPos", sample_begin);
    frame->Put("BaselineSampleEndPos", sample_end);
}

void CCMBaselineAnalyzer::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() == I3Frame::Geometry) {
        PushFrame(frame);
        return;
    }

    if(frame->GetStop() != I3Frame::DAQ) {
        PushFrame(frame);
        return;
    }

    n_daq_frames += 1;

    if(thresholds_tuned_) {
        std::cout << "Computing baseline for frame " << n_daq_frames << std::endl;
        AddBaselineStats(frame);
        PushFrame(frame);
    } else {
        if(cached_frames.size() >= num_triggers_for_threshold_ - 1) {
            cached_frames.push_back(frame);
            std::cout << "Have " << cached_frames.size() << " frames" << std::endl;
            std::cout << "Finding thresholds..." << std::endl;
            FindThresholds(cached_frames);
            thresholds_tuned_ = true;
            std::cout << "Found thresholds!" << std::endl;
            for(size_t i=0; i<current_window_size_.size(); ++i) {
                std::cout << "Threshold " << i << "/" << current_window_size_.size() << ": " << current_derivative_threshold_[i] << " Window size: " << current_window_size_[i] << std::endl;
            }
            std::cout << "Pushing frames..." << std::endl;
            for(size_t i=0; i<cached_frames.size(); ++i) {
                std::cout << "Computing baseline for frame " << i + 1 << std::endl;
                AddBaselineStats(cached_frames[i]);
                PushFrame(cached_frames[i]);
            }
            cached_frames.clear();
        } else {
            cached_frames.push_back(frame);
        }
    }
}

void CCMBaselineAnalyzer::Finish() {
    // log_notice_stream(
    //         "Triggers seen: " << triggers_seen << "\n" <<
    //         "Merged triggers ouput: " << merged_triggers_output << "\n" <<
    //         "Total triggers output: " << total_triggers_output
    // );
    // log_notice_stream(
    //     "Merged " << offsets.size() << " DAQ streams into " << counter << " separate triggers. Encountered "
    //     << counter-incomplete_counter << " complete triggers and " << incomplete_counter << " incomplete triggers");
}
