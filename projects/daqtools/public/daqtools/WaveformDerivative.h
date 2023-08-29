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

    PairResult GetDerivative() const {
        return {derivative.cbegin(), derivative.cbegin() + index+1};
    }

    PairResult GetFullDerivative() const {
        return {derivative.cbegin(), derivative.cend()};
    }

    int CurrentIndex() const {
        return index;
    }

    double Value() const {
        return begin[index];
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

#endif // I3_WaveformDerivative_H
