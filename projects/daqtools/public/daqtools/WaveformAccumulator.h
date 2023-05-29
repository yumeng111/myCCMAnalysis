#ifndef I3_WaveformAccumulator_H
#define I3_WaveformAccumulator_H

#include <deque>
#include <string>
#include <vector>

class WaveformAccumulator {
private:
    bool use_kahan_summation;
    bool first;
    std::deque<double> sum;
    std::deque<unsigned int> counts;
    std::deque<double> compensation;
    int sum_fixed_position;

    void AddElementBasic(size_t i, double val, unsigned int count) {
        sum[i] += val;
        counts[i] += count;
    }
    void AddElementKahan(size_t i, double val, unsigned int count) {
        double y = val - compensation[i];
        double t = sum[i] + y;
        compensation[i] = (t - sum[i]) - y;
        sum[i] = t;
        counts[i] += count;
    }
    void AddElement(size_t i, double val, unsigned int count) {
        if(use_kahan_summation)
            AddElementKahan(i, val, count);
        else
            AddElementBasic(i, val, count);
    }
    void AddElement(size_t i, double val) {
        AddElement(i, val, 1);
    }
public:
    WaveformAccumulator(bool use_kahan_summation = true) :
        use_kahan_summation(use_kahan_summation),
        first(true),
        sum_fixed_position(0) {}

    void AddWaveform(std::vector<double> const & wf, int wf_fixed_position) {
        AddWaveform(wf, wf_fixed_position, std::vector<unsigned int>(wf.size(), 1));
    }
    void AddWaveform(std::vector<double> const & wf, int wf_fixed_position, std::vector<unsigned int> wf_counts) {
        if(first) {
            std::copy(wf.begin(), wf.end(), std::back_inserter(sum));
            counts = std::deque<unsigned int>(sum.size(), 1);
            if(use_kahan_summation) {
                compensation = std::deque<double>(sum.size(), 0);
                sum_fixed_position = wf_fixed_position;
            }
            first = false;
        } else {
            size_t sum_length = sum.size();
            size_t wf_length = wf.size();

            int sum_start = 0;
            int wf_start = 0 + int(sum_fixed_position) - int(wf_fixed_position);
            size_t start_extension = 0;
            if(wf_start < sum_start) {
                start_extension = sum_start - wf_start;
                for(size_t i=0; i<start_extension; ++i) {
                    sum.push_front(0);
                    compensation.push_front(0);
                    counts.push_front(0);
                }
                wf_start = 0;
                sum_start = start_extension;
                sum_fixed_position = wf_fixed_position;
            }

            int sum_end = sum_start + sum_length;
            int wf_end = wf_start + wf_length;
            size_t end_extension = 0;
            if(wf_end > sum_end) {
                end_extension = wf_end - sum_end;
                for(size_t i=0; i<end_extension; ++i) {
                    sum.push_back(0);
                    compensation.push_back(0);
                    counts.push_back(0);
                }
            }

            for(size_t i=0; i<wf_length; ++i) {
                AddElement(i + wf_start, wf[i], wf_counts[i]);
            }
        }
    }

    std::deque<double> GetSummedWaveform() const {
        return sum;
    }

    std::deque<unsigned int> GetCounts() const {
        return counts;
    }

    int GetFixedPosition() const {
        return sum_fixed_position;
    }
};

#endif // I3_WaveformAccumulator_H
