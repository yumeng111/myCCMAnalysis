#ifndef I3_OnlineRobustStats_H
#define I3_OnlineRobustStats_H

#include <set>
#include <deque>
#include <algorithm>

#include <icetray/robust_statistics.h>

class OnlineRobustStats {
    std::deque<double> buffer;
    std::multiset<double> sorted_samples;
public:
    OnlineRobustStats() {}

    void AddValue(double x) {
        buffer.push_back(x);
        sorted_samples.insert(x);
    }

    void RemoveValue() {
        if(buffer.size() == 0)
            return;
        {
            double const & to_remove = buffer.front();
            sorted_samples.erase(to_remove);
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

    double Stddev(double median) {
        return robust_stats::MedianAbsoluteDeviation(
                sorted_samples.begin(),
                sorted_samples.end(),
                median);
    }
};

class OnlineRobustStatsBatched {
    std::deque<std::vector<double>> buffer;
    std::multiset<double> sorted_samples;
    bool mode_computed = false;
    bool median_computed = false;
    bool stddev_computed = false;
    bool mean_computed = false;
    double mode;
    double median;
    double stddev;
    double stddev_cached_x;
    double mean;
    void ResetState() {
        mode_computed = false;
        median_computed = false;
        stddev_computed = false;
        mean_computed = false;
    }
public:
    OnlineRobustStatsBatched() {}

    void AddValue(double x) {
        buffer.push_back({x});
        sorted_samples.insert(x);
        ResetState();
    }

    void AddValues(std::vector<double> const & x) {
        buffer.push_back(x);
        sorted_samples.insert(x.begin(), x.end());
        ResetState();
    }

    void AddValues(I3Vector<double> const & x) {
        buffer.push_back(x);
        sorted_samples.insert(x.begin(), x.end());
        ResetState();
    }

    void RemoveValue() {
        if(buffer.size() == 0)
            return;
        {
            std::vector<double> const & to_remove = buffer.front();
            for(double const & x : to_remove) {
                sorted_samples.erase(x);
            }
        }
        buffer.pop_front();
        ResetState();
    }

    double Median() {
        if(median_computed)
            return median;
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
        median_computed = true;
        median = median;
        return ret;
    }

    double Mode() {
        if(mode_computed)
            return mode;
        double ret = robust_stats::Mode(sorted_samples.begin(), sorted_samples.end());
        mode_computed = true;
        mode = ret;
        return ret;
    }

    double Stddev(double x) {
        if(stddev_computed and x == stddev_cached_x)
            return stddev;
        double ret = robust_stats::MedianAbsoluteDeviation(
                sorted_samples.begin(),
                sorted_samples.end(),
                x);
        stddev_computed = true;
        stddev = ret;
        stddev_cached_x = x;
        return ret;
    }

    double Mean() {
        if(mean_computed)
            return mean;
        double ret = 0.0;
        std::multiset<double>::iterator it = sorted_samples.begin();
        std::multiset<double>::iterator end = sorted_samples.end();
        size_t N = sorted_samples.size();
        while(it != end) {
            ret += *it;
            ++it;
        }
        mean_computed = true;
        ret /= N;
        mean = ret;
        return ret;
    }

    size_t NBatches() {
        return buffer.size();
    }

    size_t NSamples() {
        return sorted_samples.size();
    }
};

#endif // I3_OnlineRobustStats_H
