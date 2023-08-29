#ifndef I3_WaveformDerivative_H
#define I3_WaveformDerivative_H

#include <vector>
#include <utility>
#include <algorithm>
#include <icetray/I3DefaultName.h>
#include <icetray/I3FrameObject.h>
#include <icetray/I3PointerTypedefs.h>
#include <icetray/serialization.h>
#include <icetray/I3Logging.h>

template<typename T>
class WaveformDerivative {
    typename std::vector<T>::const_iterator begin;
    typename std::vector<T>::const_iterator end;
    typedef std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator> PairResult;
    size_t N;
    std::vector<double> derivative;
    size_t index;
    size_t max_computed;
    double delta_t;
    double current_value;
    double current_derivative;
    double x;
public:
    template<typename Iter>
    WaveformDerivative(Iter begin, Iter end, double delta_t) :
        begin(begin), end(end), N(std::distance(begin, end)), derivative(1), index(0), max_computed(0), delta_t(delta_t) {

        // Derivative for the first element
        derivative[0] = (-3.0 * begin[0] + 4.0 * begin[1] - begin[2]) / (2.0 * delta_t);

        index = 0;

        current_value = begin[index];
        current_derivative = derivative[index];
    }

    void Reset() {
        index = 0;
    }

    void Reset(size_t reset_index) {
        if(reset_index > max_computed) {
            Reset(max_computed);
            while(max_computed < reset_index) {
                Next();
                if(max_computed == N-1)
                    break;
            }
        }
        index = std::min(reset_index, max_computed);
    }

    PairResult GetRawWaveform() const {
        return {begin, begin + index+1};
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

    double Value() const {
        return begin[index];
    }

    double Baseline() const {
        return 0.0;
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
        if(idx < N-1) {
            // Finite difference derivative
            derivative.push_back((begin[idx+1] - begin[idx-1]) / (2.0 * delta_t));
        } else if(idx == N-1) {
            derivative.resize(N);

            // Derivative for the last element
            derivative[N-1] = (3.0 * begin[N-1] - 4.0 * begin[N-2] + begin[N-3]) / (2.0 * delta_t);
        }
        current_value = begin[idx];
        current_derivative = derivative[idx];
        max_computed = idx;
    }

    size_t Size() {
        return N;
    }
};

template<>
class WaveformDerivative<uint16_t> {
    std::vector<uint16_t>::const_iterator begin;
    std::vector<uint16_t>::const_iterator end;
    typedef std::pair<std::vector<uint16_t>::const_iterator, std::vector<uint16_t>::const_iterator> PairResult;
    size_t N;
    std::vector<double> derivative;
    std::vector<double> baselines;
    size_t index;
    size_t max_computed;
    double delta_t;
    double current_value;
    double current_derivative;
    double current_baseline;
    double x;
public:
    template<typename Iter>
    WaveformDerivative(Iter begin, Iter end, double delta_t) :
        begin(begin), end(end), N(std::distance(begin, end)), derivative(1), baselines(1), index(0), max_computed(0), delta_t(delta_t) {

        // Derivative for the first element
        derivative[0] = (-3.0 * begin[0] + 4.0 * begin[1] - begin[2]) / (2.0 * delta_t);
        baselines[0] = 0.0;
        for(size_t i=0; i<5; ++i) {
            baselines[0] += begin[i];
        }
        baselines[0] /= 5;

        index = 0;

        current_value = begin[index];
        current_derivative = derivative[index];
        current_baseline = baselines[index];
    }

    void Reset() {
        index = 0;
    }

    void Reset(size_t reset_index) {
        if(reset_index > max_computed) {
            Reset(max_computed);
            while(max_computed < reset_index) {
                Next();
                if(max_computed == N-1)
                    break;
            }
        }
        index = std::min(reset_index, max_computed);
    }

    PairResult GetRawWaveform() const {
        return {begin, begin + index+1};
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

    double Value() const {
        return -begin[index];
    }

    double Baseline() const {
        return -baselines[index];
    }

    double Derivative() const {
        return -derivative[index];
    }

    void Next() {
        if(index < N-1)
            index += 1;
        if(index <= max_computed)
            return;
        size_t idx = index;
        if(idx < N-1) {
            // Finite difference derivative
            derivative.push_back((begin[idx+1] - begin[idx-1]) / (2.0 * delta_t));
            if(idx >= 5) {
                current_baseline -= begin[idx-1] / 5.0;
                current_baseline += begin[idx] / 5.0;
            }
        } else if(idx == N-1) {
            derivative.resize(N);

            // Derivative for the last element
            derivative[N-1] = (3.0 * begin[N-1] - 4.0 * begin[N-2] + begin[N-3]) / (2.0 * delta_t);
            current_baseline -= begin[N-2] / 5.0;
            current_baseline += begin[N-1] / 5.0;
        }
        current_value = begin[idx];
        current_derivative = derivative[idx];
        baselines.push_back(current_baseline);
        max_computed = idx;
    }

    size_t Size() {
        return N;
    }
};

#endif // I3_WaveformDerivative_H
