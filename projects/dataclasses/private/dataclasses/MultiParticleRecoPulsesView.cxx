/**
 *  $Id$
 *
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#include <functional>

#include "dataclasses/MultiParticleRecoPulsesView.h"
#include "dataclasses/I3MapCCMPMTKeyUnion.h"

#include "icetray/CCMPMTKey.h"
#include "dataclasses/physics/CCMRecoPulse.h"

MultiParticleRecoPulsesView::MultiParticleRecoPulsesView() : key_(), unified_() {}

MultiParticleRecoPulsesView::MultiParticleRecoPulsesView(I3Frame const & frame,
    std::string const & key) : key_(key), unified_() {}

CCMRecoPulseSeriesMapConstPtr MultiParticleRecoPulsesView::Apply(const I3Frame &frame) const {
    typedef I3Map<CCMPMTKey, std::vector<CCMRecoPulse>> MapType;
    typedef I3Map<I3ParticleID, std::map<CCMPMTKey, std::vector<CCMRecoPulse>>> InputType;
	if(unified_) return unified_;

    boost::shared_ptr<InputType const> input = frame.Get<boost::shared_ptr<InputType const>>(key_);
    if(!input) log_fatal("Couldn't find '%s' in the frame!", key_.c_str());
    if(input->size() <= 1) {
        unified_ = boost::make_shared<MapType>(input->begin()->second);
        return unified_;
    }

    std::map<CCMPMTKey, std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>>> sources_by_key;
    std::map<CCMPMTKey, size_t> total_sizes;
    for(auto const & [id, pmap] : *input) {
        for(auto const & [k, pulses] : pmap) {
            sources_by_key[k].push_back(std::cref(pulses));
            total_sizes[k] += pulses.size();
        }
    }
	unified_ = CCMRecoPulseSeriesMapUnion::MergePulseSeries(sources_by_key, total_sizes);

    return unified_;
}

bool MultiParticleRecoPulsesView::operator==(const MultiParticleRecoPulsesView& other) const {
	return key_ == other.key_;
}

bool MultiParticleRecoPulsesView::operator!=(const MultiParticleRecoPulsesView& other) const {
	return key_ != other.key_;
}

template <class Archive>
void MultiParticleRecoPulsesView::serialize(Archive& ar, unsigned version) {
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Key", key_);
}

std::ostream& MultiParticleRecoPulsesView::Print(std::ostream& oss) const {
	oss << "MultiParticleRecoPulsesView Key: " << key_;
	return oss;
}

std::ostream& operator<<(std::ostream& os, const MultiParticleRecoPulsesView& un) {
	return(un.Print(os));
}

I3_SERIALIZABLE(MultiParticleRecoPulsesView);
