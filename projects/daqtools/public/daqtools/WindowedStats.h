#ifndef I3_WindowedStats_H
#define I3_WindowedStats_H

#include <set>
#include <deque>
#include <algorithm>

#include <icetray/robust_statistics.h>

/*
def weighted_incremental_variance(data_weight_pairs):
    w_sum = w_sum2 = mean = S = 0

    for x, w in data_weight_pairs:
        w_sum = w_sum + w
        w_sum2 = w_sum2 + w**2
        mean_old = mean
        mean = mean_old + (w / w_sum) * (x - mean_old)
        S = S + w * (x - mean_old) * (x - mean)

    population_variance = S / w_sum
    # Bessel's correction for weighted samples
    # Frequency weights
    sample_frequency_variance = S / (w_sum - 1)
    # Reliability weights
    sample_reliability_variance = S / (w_sum - w_sum2 / w_sum)
 */

namespace {
    class KahanAccumulator {
        double sum = 0;
        double c = 0;
    public:
        void Add(double x) {
            double const y = x - c;
            double const t = sum + y;
            c = (t - sum) - y;
            sum = t;
        }
        double Sum() const { return sum; }
        operator double() const { return sum; }
        operator double&() { return sum; }
        double operator+=(double x) { Add(x); return sum; }
        double operator-=(double x) { Add(-x); return sum; }
    };
}

class WindowedStats {
    KahanAccumulator w_sum;
    KahanAccumulator w_sum2;
    KahanAccumulator mean;
    KahanAccumulator S;
    std::deque<double> buffer;
    std::multiset<double> sorted_samples;
public:
    WindowedStats() {}

    void AddValue(double x) {
        buffer.push_back(x);
        sorted_samples.insert(x);

        w_sum += 1;
        w_sum2 += 1;
        double const mean_old = mean;
        mean += (1 / w_sum) * (x - mean_old);
        S += 1 * (x - mean_old) * (x - mean);
    }

    void RemoveValue() {
        if(buffer.size() == 0)
            return;
        {
            double const & to_remove = buffer.front();
            sorted_samples.erase(to_remove);
            w_sum -= 1;
            w_sum2 -= 1;
            double const mean_old = mean;
            mean -= (1 / w_sum) * (to_remove - mean_old);
            S -= 1 * (to_remove - mean_old) * (to_remove - mean);
        }
        buffer.pop_front();
    }

    double Median() {
        const size_t N = sorted_samples.size();
        const size_t half = N / 2;
        std::multiset<double>::iterator it = sorted_samples.begin();
        for(size_t i=0; i<(half-1); ++i)
            ++it;

        // Odd count: return middle
        if(N % 2) {
            ++it;
            return *it;
        }

        // Even count: return average of middle two.
        double ret = *it;
        ++it;
        ret += *it;
        ret /= 2;
        return ret;
    }

    double Mode() {
        return robust_stats::Mode(sorted_samples.begin(), sorted_samples.end());
    }

    double MedianAbsoluteDeviation(double median) {
        return robust_stats::MedianAbsoluteDeviation(
                sorted_samples.begin(),
                sorted_samples.end(),
                median);
    }

    double Mean() {
        return mean;
    }

    double Variance() {
        return S / w_sum;
    }

    double SampleVariance() {
        return S / (w_sum - 1);
    }

    double ReliabilityVariance() {
        return S / (w_sum - w_sum2 / w_sum);
    }

    double StandardDeviation() {
        return sqrt(Variance());
    }

    double SampleStandardDeviation() {
        return sqrt(SampleVariance());
    }

    double ReliabilityStandardDeviation() {
        return sqrt(ReliabilityVariance());
    }
};

#endif // I3_WindowedStats_H
