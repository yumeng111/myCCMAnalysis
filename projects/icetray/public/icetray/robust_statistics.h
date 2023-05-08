// Copyright 2017 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef ROBUST_STATISTICS_H_
#define ROBUST_STATISTICS_H_

// Robust statistics: Mode, Median, MedianAbsoluteDeviation.

#include <stddef.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <vector>

namespace robust_stats {

// @return i in [idx_begin, idx_begin + half_count) that minimizes
// sorted[i + half_count] - sorted[i].
template <typename T>
size_t MinRange(const T* const sorted, const size_t idx_begin,
        const size_t half_count) {
    T min_range = std::numeric_limits<T>::max();
    size_t min_idx = 0;

    for (size_t idx = idx_begin; idx < idx_begin + half_count; ++idx) {
        assert(sorted[idx] <= sorted[idx + half_count]);
        const T range = sorted[idx + half_count] - sorted[idx];

        if (range < min_range) {
            min_range = range;
            min_idx = idx;
        }
    }

    return min_idx;
}

template <class Iterator, class U = typename std::iterator_traits<Iterator>::value_type>
size_t MinRange(Iterator begin, Iterator end, const size_t idx_begin) {
    Iterator half_begin = end;
    U min_range = std::numeric_limits<U>::max();
    size_t min_idx = 0;
    size_t idx = idx_begin;

    while(begin != end) {
        assert(*begin <= *half_begin);
        const U range = *half_begin - *begin;

        if (range < min_range) {
            min_range = range;
            min_idx = idx;
        }
        ++begin;
        ++half_begin;
        ++idx;
    }

    return min_idx;
}

// Returns an estimate of the mode by calling MinRange on successively
// halved intervals. "sorted" must be in ascending order. This is the
// Half Sample Mode estimator proposed by Bickel in "On a fast, robust
// estimator of the mode", with complexity O(N log N). The mode is less
// affected by outliers in highly-skewed distributions than the median.
// The averaging operation below assumes "T" is an unsigned integer type.
template <typename T>
T Mode(const T* const sorted, const size_t num_values) {
    size_t idx_begin = 0;
    size_t half_count = num_values / 2;

    while (half_count > 1) {
        idx_begin = MinRange(sorted, idx_begin, half_count);
        half_count >>= 1;
    }

    const T x = sorted[idx_begin + 0];

    if (half_count == 0) {
        return x;
    }

    assert(half_count == 1);
    const T average = (x + sorted[idx_begin + 1] + 1) / 2;
    return average;
}

template <class Iterator, class U = typename std::iterator_traits<Iterator>::value_type>
U Mode(Iterator begin, Iterator end) {
    size_t num_values = std::distance(begin, end);
    size_t idx_begin = 0;
    size_t half_count = num_values / 2;

    while (half_count > 1) {
        idx_begin = MinRange(begin + idx_begin, begin + idx_begin + half_count, idx_begin);
        half_count >>= 1;
    }

    const U x = *(begin + idx_begin);

    if (half_count == 0) {
        return x;
    }

    assert(half_count == 1);
    const U average = (x + *(begin + idx_begin + 1) + 1) / 2;
    return average;
}

// Sorts integral values in ascending order. About 3x faster than std::sort for
// input distributions with very few unique values.
template <class T>
void CountingSort(T* begin, T* end) {
    // Unique values and their frequency (similar to flat_map).
    using Unique = std::pair<T, int>;
    std::vector<Unique> unique;

    for (const T* p = begin; p != end; ++p) {
        const T value = *p;
        const auto pos =
            std::find_if(unique.begin(), unique.end(),
                    [value](const Unique & u) {
                    return u.first == value;
                    });

        if (pos == unique.end()) {
            unique.push_back(std::make_pair(*p, 1));
        } else {
            ++pos->second;
        }
    }

    // Sort in ascending order of value (pair.first).
    std::sort(unique.begin(), unique.end());
    // Write that many copies of each unique value to the array.
    T* p = begin;

    for (const auto& value_count : unique) {
        std::fill(p, p + value_count.second, value_count.first);
        p += value_count.second;
    }

    assert(p == end);
}

// Returns the median value. Side effect: sorts "samples".
template <typename T>
T Median(std::vector<T>* samples) {
    assert(!samples->empty());
    std::sort(samples->begin(), samples->end());
    const size_t half = samples->size() / 2;

    // Odd count: return middle
    if (samples->size() % 2) {
        return (*samples)[half];
    }

    // Even count: return average of middle two.
    return ((*samples)[half] + (*samples)[half - 1]) / 2;
}

// Returns a robust measure of variability.
template <typename T>
T MedianAbsoluteDeviation(const std::vector<T>& samples, const T median) {
    assert(!samples.empty());
    std::vector<T> abs_deviations;
    abs_deviations.reserve(samples.size());

    for (const T sample : samples) {
        abs_deviations.push_back(std::abs(sample - median));
    }

    return Median(&abs_deviations);
}

template <class Iterator, class U = typename std::iterator_traits<Iterator>::value_type>
U MedianAbsoluteDeviation(Iterator begin, Iterator end, const U median) {
    assert(begin != end);
    size_t size = std::distance(begin, end);
    std::vector<U> abs_deviations;
    abs_deviations.reserve(size);

    while(begin != end) {
        abs_deviations.push_back(std::abs(*begin - median));
        ++begin;
    }

    return Median(&abs_deviations);
}

}  // namespace robust_stats

#endif  // ROBUST_STATISTICS_H_

