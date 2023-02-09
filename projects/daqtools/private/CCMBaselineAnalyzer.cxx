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
    virtual T output();
};

template<class T>
class WaveformSource : virtual public Source<T> {
protected:
    std::vector<T> const & source;
    size_t current_pos;
public:
    WaveformSource(CCMWaveform<T> const & source_wf) :
        source(source_wf.GetWaveform()), current_pos(0) {}
    virtual inline T output() {
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
class Filter : virtual public Source<U> {
protected:
    Source<T> * input;
public:
    Filter() {}
    Filter(Source<T> * input) : input(input) {}
    virtual inline U filter(T const & t);
    virtual inline U output() {
        return filter(input->output());
    }
};

template<class T, class U>
class ExponentialFilter : virtual public Filter<T, U> {
protected:
    U tau;
    U delta_t;
    U y_i;
    U alpha;
public:
    ExponentialFilter(U tau, U delta_t, U y_0, Source<T> * input) :
        tau(tau), delta_t(delta_t), y_i(y_0), Filter<T, U>(input) {
        alpha = exp(-delta_t / tau);
    }

    virtual inline U filter(T const & t) {
        y_i = y_i + (1.0 - alpha) * (t - y_i);
        return y_i;
    }
};

template<class T, class U>
class MultiplicativeFilter : virtual public Filter<T, U> {
protected:
    U c;
public:
    MultiplicativeFilter(U c, Source<T> * input) :
        c(c), Filter<T, U>(input) {}
    virtual inline U filter(T const & t) {
        return U(t)*c;
    }
};

template<class T, class U>
class Buffer : virtual public Filter<T, U> {
protected:
    std::deque<T> buffer;
    size_t max_size;
    bool buffer_filling;
    bool source_exhausted;
public:
    Buffer(std::deque<T> init_buffer, size_t max_size, Source<T> * input) :
        buffer(init_buffer), max_size(max_size), buffer_filling(true), source_exhausted(false), Filter<T, U>(input) {
        if(buffer.size() >= max_size)
            buffer_filling = false;
        size_t to_pop = buffer.size() - std::max(max_size, buffer.size());
        for(size_t i=0; i<to_pop; ++i) {
            buffer.pop_front();
        }
    }
    virtual inline U output() {
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
    virtual inline U filter(std::deque<T> const & buffer);
};

template<class T, class U>
class BufferView : virtual public Filter<T, U> {
protected:
    std::deque<T> buffer;
    size_t buffer_begin;
    size_t buffer_end;
public:
    BufferView(Source<T> * input) : Filter<T, U>(input) {
        buffer_begin = 0;
        buffer_end = 0;
    }
    virtual inline U filter(T const & t) {
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
class SmoothingFilter : virtual public Buffer<T, U> {
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
    virtual inline U filter(std::deque<T> const & buffer) {
        if(Buffer<T, U>::buffer_filling) {
            size_t n_excluded = coeff.size() - buffer.size();
            U this_total = 0;
            U sum = 0;
            for(size_t i=0; i<buffer.size(); ++i) {
                size_t j = i + n_excluded;
                this_total += coeff[j];
                sum += coeff[j] * buffer[i];
            }
            if(this_total == 0) {
                return 0;
            }
            return sum / this_total;
        } else if(Buffer<T, U>::source_exhausted) {
            if(buffer.size() <= zero_pos)
                throw SourceExhausted();
            U this_total = 0;
            U sum = 0;
            for(size_t i=0; i<buffer.size(); ++i) {
                this_total += coeff[i];
                sum += coeff[i] * buffer[i];
            }
            if(this_total == 0) {
                return 0;
            }
            return sum / this_total;

        } else {
            U sum = 0;
            for(size_t i=0; i<coeff.size(); ++i) {
                sum += coeff[i] * buffer[i];
            }
            return sum / total;
        }
    }
};

template<class T, class U>
class DerivativeFilter : virtual public Filter<T, U> {
protected:
    U bin_width;
    size_t n_seen;
    U t_m2;
    U t_m1;
public:
    DerivativeFilter(U bin_width, Source<T> * input) :
        bin_width(bin_width), n_seen(0), Filter<T, U>(input) {}
    virtual inline U filter(T const & t) {
        U x;
        if(n_seen >= 3) {
            x = (t - t_m2) / (bin_width * 2.0);
        } else {
            x = 0;
        }
        t_m2 = t_m1;
        t_m1 = t;
        return x;
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
}

void CCMBaselineAnalyzer::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMDAQConfigName", daq_config_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
}

template<class T, class U, class K>
struct DerivativeWindowCapture {
    Source<K> & d_buffer;
    BufferView<T, U> & raw_buffer;
    K max_derivative;
    int state;
    // States:
    // 0 - outside of region without activity
    // 1 - inside of region without activity

    DerivativeWindowCapture(Source<K> & d_buffer, BufferView<T, U> raw_buffer, K max_derivative) :
        d_buffer(d_buffer), raw_buffer(raw_buffer), max_derivative(max_derivative), state(0) {}

    bool Inside(T const & t) const {
        return t <= max_derivative and t >= -max_derivative;
    }

    bool Next() {
        K derivative = d_buffer.output();
        bool inside = Inside(derivative);
        if(state == 0) { // currently outside region
            if(inside) { // moving inside region
                raw_buffer.Reset();
                return false;
            } else { // still outside region
                return false;
            }
        } else { // already inside
            if(not inside) { // moving outside region
                return true;
            } else { // still inside region
                return false;
            }
        }
    }

    std::deque<K> const & GetBuffer() const {
        return raw_buffer.GetBuffer();
    }
};

struct WindowStats {
    double mean;
    double variance;
    double k_;
    double M_;
    double S_;
    WindowStats() {}
    template<typename T>
    WindowStats(std::deque<T> const & buffer) {
        double k = 0;
        double M = 0;
        double S = 0;
        for(size_t i=0; i<buffer.size(); ++i) {
            k += 1;
            double const & x = buffer[i];
            double new_M = M + (x + M)* 1.0 / k;
            double new_S = S + (x - M) * (x - new_M);
            M = new_M;
            S = new_S;
        }
        mean = M;
        variance = S / (k - 1);
        k_ = k; M_ = M; S_ = M;
    }

    WindowStats(std::vector<WindowStats> const & stats) {
        double k = 0;
        double M = 0;
        double S = 0;
        for(size_t i=0; i<stats.size(); ++i) {
            WindowStats const & stat = stats[i];
            double k_a = k;
            k += stat.k_;
            double delta = stat.S_ - S;
            M = M + delta * stat.k_ / k;
            S += stat.M_ + delta * delta * k_a * stat.k_ / k;
        }
        mean = M;
        variance = S / (k - 1);
        k_ = k; M_ = M; S_ = M;
    }
};

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

    double max_derivative = 0.65;

    for(size_t i=0; i<waveforms.size(); ++i) {
        CCMWaveformUInt16 const & wf = waveforms[i];
        WaveformSource<uint16_t> wf_source(wf);
        MultiplicativeFilter<uint16_t, double> m_filter(-1.0, &wf_source);
        BufferView<double, double> buffer_m(&m_filter);
        ExponentialFilter<double, double> exp_filter(10.0, 2.0, 0.0, &buffer_m);
        std::deque<double> smoothing_kernel = {1,1,1};
        SmoothingFilter<double, double> smoothing_filter(1, smoothing_kernel, &exp_filter);
        DerivativeFilter<double, double> derivative_filter(2.0, &smoothing_filter);
        DerivativeWindowCapture<double, double, double> window(derivative_filter, buffer_m, max_derivative);

        std::vector<WindowStats> stats;

        try {
            while(true) {
                bool have_window = window.Next();
                if(have_window) {
                    stats.emplace_back(window.GetBuffer());
                }
            }
        } catch(SourceExhausted const & e) {
        }

        WindowStats stat(stats);

        /*
         * Do something with this information...
         */
        std::cout << "N regions found: " << stats.size() << std::endl;
        std::cout << "Mean ADC: " << stat.mean << std::endl;
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
