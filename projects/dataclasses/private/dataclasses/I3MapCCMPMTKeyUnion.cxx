/**
 *  $Id$
 *
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#include <queue>
#include <memory>
#include <algorithm>
#include <functional>

#include "dataclasses/I3MapCCMPMTKeyUnion.h"

#include "icetray/CCMPMTKey.h"
#include "dataclasses/physics/CCMRecoPulse.h"

CCMRecoPulseSeriesMapUnion::CCMRecoPulseSeriesMapUnion() : keys_(), unified_() {}

CCMRecoPulseSeriesMapUnion::CCMRecoPulseSeriesMapUnion(const I3Frame &frame,
    const std::vector<std::string> &keys) : keys_(keys), unified_() {}

inline bool time_less(CCMRecoPulse const & a, CCMRecoPulse const & b) {
    return a.GetTime() < b.GetTime();
}

inline bool is_sorted(std::vector<CCMRecoPulse> const & vec) {
    for(size_t i = 1; i < vec.size(); ++i) {
        if(time_less(vec[i], vec[i-1]))
            return false;
    }
    return true;
}

inline bool time_equal(CCMRecoPulse const & a, CCMRecoPulse const & b) {
    return a.GetTime() == b.GetTime();
}

inline void compact_by_time(std::vector<CCMRecoPulse> & vec) {
    if(vec.empty()) return;
    size_t dest = 0;
    for(size_t i = 1; i < vec.size(); ++i) {
        if(time_equal(vec[i], vec[dest])) {
            vec[dest].SetCharge(vec[dest].GetCharge() + vec[i].GetCharge());
            vec[dest].SetWidth (vec[dest].GetWidth()  + vec[i].GetWidth());
        } else {
            ++dest;
            if(i != dest)
                vec[dest] = std::move(vec[i]);
        }
    }
    vec.resize(dest + 1);
}

void HeapMerge(
        std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>> const & srcs,
        std::vector<CCMRecoPulse> & dest
        ) {

    struct View {
        const CCMRecoPulse* data;
        size_t size;
        std::unique_ptr<std::vector<CCMRecoPulse>> owned_buf;
    };
    std::vector<View> views;
    views.reserve(srcs.size());
    for(std::vector<CCMRecoPulse> const & vec : srcs) {
        if(vec.empty()) {
            views.push_back(View{nullptr, 0, nullptr});
            continue;
        }
        if(vec.size() < 2 or is_sorted(vec)) {
            views.push_back(View{vec.data(), vec.size(), nullptr});
        } else {
            std::unique_ptr<std::vector<CCMRecoPulse>> buf = std::make_unique<std::vector<CCMRecoPulse>>(vec);
            std::sort(buf->begin(), buf->end(), time_less);
            views.push_back(View{buf->data(), buf->size(), std::move(buf)});
        }
    }

    // Min-heap of (time, source index, element index)
    struct Node {
        double t; // time
        size_t s; // source index
        size_t i; // element index within source
    };
    struct Cmp {
        bool operator()(Node const& a, Node const& b) const {
            if (a.t != b.t) return a.t > b.t; // min-heap by time
            // tie-breaker to keep it strict-weak-ordered
            if (a.s != b.s) return a.s > b.s;
            return a.i > b.i;
        }
    };
    std::priority_queue<Node, std::vector<Node>, Cmp> pq;

    // Seed heap
    for(size_t s = 0; s < srcs.size(); ++s) {
        if(views[s].size > 0) {
            pq.push(Node{views[s].data[0].GetTime(), s, 0});
        }
    }

    while(not pq.empty()) {
        Node current = pq.top(); pq.pop();
        CCMRecoPulse const & p = views[current.s].data[current.i];

        if(not dest.empty() and time_equal(p, dest.back())) {
            // accumulate into last
            dest.back().SetCharge(dest.back().GetCharge() + p.GetCharge());
            dest.back().SetWidth (dest.back().GetWidth()  + p.GetWidth());
        } else {
            dest.push_back(p);
        }

        size_t next_i = current.i + 1;
        if(next_i < views[current.s].size) {
            pq.push(Node{views[current.s].data[next_i].GetTime(), current.s, next_i});
        }
    }
}

inline void SimpleMerge(
        std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>> const & srcs,
        std::vector<CCMRecoPulse> & dest
        ) {
    for(std::vector<CCMRecoPulse> const & vec : srcs) {
        dest.insert(dest.end(), vec.begin(), vec.end());
    }
    std::sort(dest.begin(), dest.end(), time_less);
    compact_by_time(dest);
}

inline void TwoWayMerge(
        std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>> const & srcs,
        std::vector<CCMRecoPulse> & dest
        ) {

    std::vector<CCMRecoPulse> const * a_ptr = &srcs[0].get();
    std::vector<CCMRecoPulse> const * b_ptr = &srcs[1].get();
    std::vector<CCMRecoPulse> a_buf;
    std::vector<CCMRecoPulse> b_buf;

    // Check if both inputs are sorted
    if(not is_sorted(*a_ptr)) {
        a_buf = srcs[0].get();
        std::sort(a_buf.begin(), a_buf.end(), time_less);
        a_ptr = &a_buf;
    }
    if(not is_sorted(*b_ptr)) {
        b_buf = srcs[1].get();
        std::sort(b_buf.begin(), b_buf.end(), time_less);
        b_ptr = &b_buf;
    }

    std::vector<CCMRecoPulse> const & a = *a_ptr;
    std::vector<CCMRecoPulse> const & b = *b_ptr;

    size_t i = 0, j = 0;
    while(i < a.size() and j < b.size()) {
        if(time_less(a[i], b[j])) {
            dest.push_back(a[i]);
            ++i;
        } else if(time_less(b[j], a[i])) {
            dest.push_back(b[j]);
            ++j;
        } else {
            // equal time
            CCMRecoPulse merged = a[i];
            merged.SetCharge(merged.GetCharge() + b[j].GetCharge());
            merged.SetWidth (merged.GetWidth()  + b[j].GetWidth());
            dest.push_back(merged);
            ++i; ++j;
        }
    }
    while(i < a.size()) {
        dest.push_back(a[i]);
        ++i;
    }
    while(j < b.size()) {
        dest.push_back(b[j]);
        ++j;
    }
}

CCMRecoPulseSeriesMapConstPtr CCMRecoPulseSeriesMapUnion::Apply(const I3Frame &frame) const {
    typedef I3Map<CCMPMTKey, std::vector<CCMRecoPulse>> MapType;
	if(unified_) return unified_;
    if(keys_.size() == 1) {
        boost::shared_ptr<MapType const> pmap = frame.Get<boost::shared_ptr<MapType const>>(keys_.front());
        if(!pmap) log_fatal("Couldn't find '%s' in the frame!", keys_.front().c_str());
        return pmap;
    }

    std::vector<boost::shared_ptr<const MapType>> pointers;
    pointers.reserve(keys_.size());
    std::map<CCMPMTKey, std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>>> sources_by_key;
    std::map<CCMPMTKey, size_t> total_sizes;
    for(std::string const & key: keys_) {
        boost::shared_ptr<MapType const> pmap = frame.Get<boost::shared_ptr<MapType const>>(key);
        if(!pmap) log_fatal("Couldn't find '%s' in the frame!", key.c_str());
        pointers.push_back(pmap);
        for(auto const & [k, pulses] : *pmap) {
            sources_by_key[k].push_back(std::cref(pulses));
            total_sizes[k] += pulses.size();
        }
    }

    if(keys_.size() == 2) {
        unified_ = boost::make_shared<MapType>();
        for(auto & entry : sources_by_key) {
            CCMPMTKey const & key = entry.first;
            std::vector<CCMRecoPulse> & dest = (*unified_)[entry.first];
            std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>> const & srcs = entry.second;
            if(srcs.size() == 2) {
                dest.reserve(total_sizes[key]);
                TwoWayMerge(srcs, dest);
            } else if(srcs.size() == 1) {
                dest = srcs.front().get();
            } // else leave empty
        }
        return unified_;
    }

    unified_ = MergePulseSeries(sources_by_key, total_sizes);
    return unified_;
}

CCMRecoPulseSeriesMapPtr CCMRecoPulseSeriesMapUnion::MergePulseSeries(
    std::map<CCMPMTKey, std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>>> const & sources_by_key,
    std::map<CCMPMTKey, size_t> total_sizes) {

    typedef I3Map<CCMPMTKey, std::vector<CCMRecoPulse>> MapType;
    boost::shared_ptr<MapType> unified = boost::make_shared<MapType>();

    for(auto & entry : sources_by_key) {
        CCMPMTKey const & key = entry.first;
        std::vector<CCMRecoPulse> & dest = (*unified)[entry.first];
        std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>> const & srcs = entry.second;
        if(srcs.size() == 1) {
            dest = srcs.front().get();
        } else if(srcs.size() == 2) {
            dest.reserve(total_sizes[key]);
            TwoWayMerge(srcs, dest);
        } else if(total_sizes[key] < 96) {
            dest.reserve(total_sizes[key]);
            SimpleMerge(srcs, dest);
        } else {
            dest.reserve(total_sizes[key]);
            HeapMerge(srcs, dest);
        }
    }

    return unified;
}

bool CCMRecoPulseSeriesMapUnion::operator==(const CCMRecoPulseSeriesMapUnion& other) const {
	return keys_ == other.keys_;
}

bool CCMRecoPulseSeriesMapUnion::operator!=(const CCMRecoPulseSeriesMapUnion& other) const {
	return keys_ != other.keys_;
}

template <class Archive>
void CCMRecoPulseSeriesMapUnion::serialize(Archive& ar, unsigned version) {
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Keys", keys_);
}

std::ostream& CCMRecoPulseSeriesMapUnion::Print(std::ostream& oss) const {
	oss << "CCMRecoPulseSeriesMapUnion Keys:";
	for(const auto& key : keys_)
		oss << "\n  " << key;
	return oss;
}

std::ostream& operator<<(std::ostream& os, const CCMRecoPulseSeriesMapUnion& un) {
	return(un.Print(os));
}

I3_SERIALIZABLE(CCMRecoPulseSeriesMapUnion);
