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
    // States:
    // 0 - outside of region without activity
    // 1 - inside of region without activity

    DerivativeWindowCapture(Source<K> & d_buffer, BufferView<T, U> & raw_buffer, K max_derivative) :
        d_buffer(d_buffer), raw_buffer(raw_buffer), max_derivative(max_derivative), state(0) {}

    bool Inside(T const & t) const {
        return t <= max_derivative and t >= -max_derivative;
    }

    bool Next() {
        K derivative = d_buffer.output();
        bool inside = Inside(derivative);
        // std::cout << "Currently inside (" << (state ? "True" : "False") << ") will be (" << (inside ? "True" : "False") << ") with derivative = " << derivative << std::endl;
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
};

struct WindowStats {
    double mean;
    double variance;
    double k_;
    double M_;
    double S_;

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
        k_(0), M_(0), S_(0) {
    }

    template<typename Iterator>
    WindowStats(Iterator begin, Iterator end) :
        k_(0), M_(0), S_(0) {
        AddSamples<Iterator>(begin, end);
    }

    template<typename T>
    WindowStats(std::deque<T> const & buffer) :
        k_(0), M_(0), S_(0) {
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

    std::string baseline_fit_output_name_;

    double initial_derivative_threshold_;
    size_t minimum_sample_length_;
    double baseline_minimum_window_fraction_;
    double baseline_maximum_window_fraction_;
    size_t baseline_sample_edge_cut_;
    size_t num_triggers_for_threshold_;

    std::map<size_t, double> curent_derivative_threshold_;

public:
    CCMBaselineAnalyzer(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMBaselineAnalyzer);

CCMBaselineAnalyzer::CCMBaselineAnalyzer(const I3Context& context) : I3Module(context) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMAnalysis::Binary::CCMDAQConfig>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("InitialDerivativeThreshold", "The threshold first used to find regions of inactivity. This is iteratively increased until regions of inactivity represent a certain fraction of the waveform.",initial_derivative_threshold_);
    AddParameter("MinimumSampleLength", "The smallest number of consecutive samples to use when finding regions of inactivity.", minimum_sample_length_);
    AddParameter("BaselineMinimumWindowFraction", "The smallest allowed fraction that regions of inactivity may occupy across the sampled waveforms. This is needed to ensure we have enough samples to do the calculation.", baseline_minimum_window_fraction_);
    AddParameter("BaselineMaximumWindowFraction", "The largest allowed fraction that regions of inactivity may occupy across the sampled waveforms. This is only present to alert us to the unlikely scenario where a PMT has too much activity to get a measurement of the baseline.", baseline_maximum_window_fraction_);
    AddParameter("BaselineSampleEdgeCut", "The number of samples to cut from the edges of each region of inactivity. This helps us avoid contamination of the baseline measurement from the adjacent regions of activity.", baseline_sample_edge_cut_);
    AddParameter("NumTriggersForThreshold","The number of triggers to use when tuning the derivative threshold.", num_triggers_for_threshold_);
    AddParameter("BaselineFitOutputName", "The output key of the baseline fit.", baseline_fit_output_name_);
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

    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series>(waveforms_name_);

    double max_derivative = initial_derivative_threshold_;

    for(size_t i=0; i<waveforms.size(); ++i) {
        CCMWaveformUInt16 const & wf = waveforms[i];
        WaveformSource<uint16_t> wf_source(wf);
        MultiplicativeFilter<uint16_t, double> m_filter(-1.0, &wf_source);
        BufferView<double, double> buffer_m(&m_filter);
        ExponentialFilter<double, double> exp_filter(10.0, 2.0, &buffer_m);
        std::deque<double> smoothing_kernel = {1,1,1};
        SmoothingFilter<double, double> smoothing_filter(1, smoothing_kernel, &exp_filter);
        BackDerivativeFilter<double, double> derivative_filter(2.0, &smoothing_filter);
        DerivativeWindowCapture<double, double, double> window(derivative_filter, buffer_m, max_derivative);

        std::vector<WindowStats> stats;
        WindowStats total;

        try {
            while(true) {
                bool have_window = window.Next();
                if(have_window) {
                    std::tuple<std::deque<double>::const_iterator, std::deque<double>::const_iterator, size_t> window_buffer = window.GetBuffer();
                    std::deque<double>::const_iterator begin = std::get<0>(window_buffer);
                    std::deque<double>::const_iterator end = std::get<1>(window_buffer);
                    size_t size = std::get<2>(window_buffer);
                    if(size >= minimum_sample_length_) {
                        stats.emplace_back(begin, end);
                        total.AddSamples(begin, end);
                        std::cout << "Window: avg(" << stats.back().mean << ") first(" << *begin << ") len(" << size << ")" << std::endl;
                    }
                }
            }
        } catch(SourceExhausted const & e) {
        }

        /*
         * Do something with this information...
         */
        double max_length = 0;
        double tot_length = 0;
        for(size_t j=0; j<stats.size(); ++j) {
            max_length = std::max(max_length, stats[j].k_);
            tot_length += stats[j].k_;
        }
        double avg_length = tot_length / stats.size();

        std::cout << "N regions found: " << stats.size() << std::endl;
        std::cout << "Max length: " << max_length << std::endl;
        std::cout << "Avg length: " << avg_length << std::endl;
        std::cout << "Mean ADC: " << total.mean << std::endl;
        std::cout << "Std dev: " << sqrt(total.variance) << std::endl;
        std::cout << std::endl;
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
