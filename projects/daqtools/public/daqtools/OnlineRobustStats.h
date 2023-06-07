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
public:
    OnlineRobustStatsBatched() {}

    void AddValue(double x) {
        buffer.push_back({x});
        sorted_samples.insert(x);
    }

    void AddValues(std::vector<double> const & x) {
        buffer.push_back(x);
        sorted_samples.insert(x.begin(), x.end());
    }

    void AddValues(I3Vector<double> const & x) {
        buffer.push_back(x);
        sorted_samples.insert(x.begin(), x.end());
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

    size_t NBatches() {
        return buffer.size();
    }

    size_t NSamples() {
        return sorted_samples.size();
    }
};

#endif // I3_OnlineRobustStats_H
