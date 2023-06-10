#ifndef I3_WaveformSmoother_H
#define I3_WaveformSmoother_H

#include <vector>
#include <utility>
#include <algorithm>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>

static const unsigned waveformsmoother_version_ = 1;
class WaveformSmoother : public I3FrameObject {
    std::vector<uint16_t>::const_iterator begin;
    std::vector<uint16_t>::const_iterator end;
    size_t N;
    std::vector<double> smoothed_wf;
    std::vector<double> derivative;
    size_t index;
    size_t max_computed;
    double delta_t;
    double tau;
    double current_value;
    double current_derivative;
    double alpha;
    double y_i;
    double exp_prevprev, exp_prev, exp_current, exp_next, exp_nextnext;
    double x;
public:
    WaveformSmoother(std::vector<uint16_t>::const_iterator begin, std::vector<uint16_t>::const_iterator end, double delta_t, double tau) :
        begin(begin), end(end), N(std::distance(begin, end)), smoothed_wf(3), derivative(1), index(0), max_computed(0), delta_t(delta_t), tau(tau) {
        int init_avg_end_idx = std::min(3 * std::max(int(tau / delta_t), 1), int(N));
        double exp_start = 0.0;
        std::vector<uint16_t>::const_iterator x_it = begin;
        for(size_t i=0; i < init_avg_end_idx and x_it != end; ++i, ++x_it) {
            exp_start += -double(*x_it);
        }
        exp_start /= init_avg_end_idx;
        alpha = exp(-delta_t / tau);
        y_i = exp_start;

        // Exponentially smooth the first element
        x = -double(begin[0]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_prevprev = y_i;

        // Exponentially smooth the second element
        x = -double(begin[1]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_prev = y_i;

        // Exponentially smooth the third element
        x = -double(begin[2]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_current = y_i;

        // Exponentially smooth the fourth element
        x = -double(begin[3]);
        y_i += (1.0 - alpha) * (x - y_i);
        exp_next = y_i;

        // Box smoothing for the first element
        smoothed_wf[0] = (exp_prevprev + exp_prev) / 2.0;

        // Box smoothing for the second element
        smoothed_wf[1] = (exp_prevprev + exp_prev + exp_current) / 3.0;

        // Box smoothing for the third element
        smoothed_wf[2] = (exp_prev + exp_current + exp_next) / 3.0;

        // Derivative for the first element
        derivative[0] = (-3.0 * smoothed_wf[0] + 4.0 * smoothed_wf[1] - smoothed_wf[2]) / (2.0 * delta_t);

        index = 0;

        current_value = smoothed_wf[index];
        current_derivative = derivative[index];
    }

    void Reset() {
        index = 0;
    }

    void Reset(size_t reset_index) {
        index = std::min(index, reset_index);
    }

    void Reset(std::vector<double>::const_iterator end) {
        index = std::min(ptrdiff_t(index), std::distance(smoothed_wf.cbegin(), end) - 1);
    }

    std::pair<std::vector<uint16_t>::const_iterator, std::vector<uint16_t>::const_iterator> GetRawWaveform() const {
        return {begin, begin + index+1};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetSmoothedWaveform() const {
        return {smoothed_wf.cbegin(), smoothed_wf.cbegin() + index+1};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetFullSmoothedWaveform() const {
        return {smoothed_wf.cbegin(), smoothed_wf.cend()};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetDerivative() const {
        return {derivative.cbegin(), derivative.cbegin() + index+1};
    }

    std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> GetFullDerivative() const {
        return {derivative.cbegin(), derivative.cend()};
    }

    int CurrentIndex() const {
        return index;
    }

    uint16_t RawValue() const {
        return begin[index];
    }

    double Value() const {
        return smoothed_wf[index];
    }

    double Derivative() const {
        return derivative[index];
    }

    void Next() {
        if(index < N-1)
            index += 1;
        if(index <= max_computed)
            return;
        size_t idx = index;
        // Exponential smoothing for the nextnext element
        // Box smoothing for the next element
        // Derivative for the current element
        if(idx < N-2) {
            x = -double(begin[idx+2]);
            y_i += (1.0 - alpha) * (x - y_i);
            exp_nextnext = y_i;

            // Box smoothing for everything in between
            //smoothed_wf[idx+1] = (exp_current + exp_next + exp_nextnext) / 3.0;
            smoothed_wf.push_back((exp_current + exp_next + exp_nextnext) / 3.0);
            exp_current = exp_next;
            exp_next = exp_nextnext;

            // Finite difference derivative
            //derivative[idx] = (smoothed_wf[idx+1] - smoothed_wf[idx-1]) / (2.0 * delta_t);
            derivative.push_back((smoothed_wf[idx+1] - smoothed_wf[idx-1]) / (2.0 * delta_t));
        } else if(idx == N-2) {
            smoothed_wf.resize(N);
            derivative.resize(N);
            // Box smoothing for the last element
            smoothed_wf[N-1] = (exp_current + exp_next) / 2.0;

            // Derivative for the second to last element
            derivative[N-2] = (smoothed_wf[N-1] - smoothed_wf[N-3]) / (2.0 * delta_t);

            // Derivative for the last element
            derivative[N-1] = (3.0 * smoothed_wf[N-1] - 4.0 * smoothed_wf[N-2] + smoothed_wf[N-3]) / (2.0 * delta_t);
        }
        current_value = smoothed_wf[idx];
        current_derivative = derivative[idx];
        max_computed = idx;
    }

    size_t Size() {
        return N;
    }
};

I3_POINTER_TYPEDEFS(WaveformSmoother);

#endif // I3_WaveformSmoother_H
