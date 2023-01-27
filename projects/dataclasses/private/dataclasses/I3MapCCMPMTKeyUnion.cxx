/**
 *  $Id$
 *  
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */

#include "dataclasses/I3MapCCMPMTKeyUnion.h"
#include "dataclasses/physics/CCMRecoPulse.h"
#include "boost/make_shared.hpp"
#include "boost/foreach.hpp"

CCMRecoPulseSeriesMapUnion::CCMRecoPulseSeriesMapUnion() : keys_(), unified_() {}

CCMRecoPulseSeriesMapUnion::CCMRecoPulseSeriesMapUnion(const I3Frame &frame,
    const std::vector<std::string> &keys) : keys_(keys), unified_() {}

static bool PulseCompare(const CCMRecoPulse &p1, const CCMRecoPulse &p2)
{
	return p1.GetTime() < p2.GetTime();
}

CCMRecoPulseSeriesMapConstPtr
CCMRecoPulseSeriesMapUnion::Apply(const I3Frame &frame) const
{
	typedef CCMRecoPulseSeriesMap Map;
	typedef boost::shared_ptr<const Map> MapConstPtr;
	typedef Map::value_type Pair;
	typedef Pair::second_type Series;
	typedef Series::value_type Element;
	
	if (unified_)
		return unified_;
	
	unified_ = boost::make_shared<Map>();
	
	BOOST_FOREACH(const std::string &key, keys_) {
		MapConstPtr pmap = frame.Get<MapConstPtr>(key);
		if (!pmap)
			log_fatal("Couldn't find '%s' in the frame!",
			    key.c_str());
		
		BOOST_FOREACH(const Pair &pair, *pmap) {
			Series &univec = (*unified_)[pair.first];
			BOOST_FOREACH(const Element &element, pair.second)
				univec.push_back(element);
		}
	}
	
	BOOST_FOREACH(Pair &pair, *unified_)
		std::sort(pair.second.begin(), pair.second.end(),
		    PulseCompare);
	
	return unified_;
}

bool
CCMRecoPulseSeriesMapUnion::operator==(const CCMRecoPulseSeriesMapUnion& other) const
{
	return keys_ == other.keys_;
}

bool
CCMRecoPulseSeriesMapUnion::operator!=(const CCMRecoPulseSeriesMapUnion& other) const
{
	return keys_ != other.keys_;
}

template <class Archive>
void
CCMRecoPulseSeriesMapUnion::serialize(Archive& ar, unsigned version)
{
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Keys", keys_);
}

std::ostream& CCMRecoPulseSeriesMapUnion::Print(std::ostream& oss) const{
	oss << "CCMRecoPulseSeriesMapUnion Keys:";
	for(const auto& key : keys_)
		oss << "\n  " << key;
	return oss;
}

std::ostream& operator<<(std::ostream& os, const CCMRecoPulseSeriesMapUnion& un){
	return(un.Print(os));
}

I3_SERIALIZABLE(CCMRecoPulseSeriesMapUnion);
