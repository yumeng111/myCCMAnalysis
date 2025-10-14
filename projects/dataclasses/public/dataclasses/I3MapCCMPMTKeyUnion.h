/**
 *  $Id$
 *
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#ifndef DATACLASSES_I3MAPCCMPMTKEYUNION_H_INCLUDED
#define DATACLASSES_I3MAPCCMPMTKEYUNION_H_INCLUDED

#include <string>
#include <functional>
#include "icetray/I3FrameObject.h"
#include "icetray/I3Frame.h"
#include "icetray/serialization.h"
#include "dataclasses/physics/CCMRecoPulse.h"

static const unsigned ccmrecopulseseriesmapunion_version_ = 0;

class CCMRecoPulseSeriesMapUnion : public I3FrameObject {
public:
	/*
	 * Construct a union for the map stored at "key." All bits are set.
	 */
    CCMRecoPulseSeriesMapUnion(const I3Frame&, const std::vector<std::string> &keys);
    CCMRecoPulseSeriesMapUnion();

	std::ostream& Print(std::ostream&) const override;

	CCMRecoPulseSeriesMapConstPtr Apply(const I3Frame&) const;
	std::vector<std::string> GetSources() const { return keys_; }

	bool operator==(const CCMRecoPulseSeriesMapUnion&) const;
	bool operator!=(const CCMRecoPulseSeriesMapUnion&) const;

    static CCMRecoPulseSeriesMapPtr MergePulseSeries(
        std::map<CCMPMTKey, std::vector<std::reference_wrapper<std::vector<CCMRecoPulse> const>>> const & sources_by_key,
        std::map<CCMPMTKey, size_t> total_sizes
        );
private:
	std::vector<std::string> keys_;
	mutable CCMRecoPulseSeriesMapPtr unified_;

	friend class icecube::serialization::access;
	template <class Archive> void serialize(Archive& ar, unsigned version);

	SET_LOGGER("CCMRecoPulseSeriesMapUnion");
};

std::ostream& operator<<(std::ostream&, const CCMRecoPulseSeriesMapUnion&);

I3_CLASS_VERSION(CCMRecoPulseSeriesMapUnion, ccmrecopulseseriesmapunion_version_);
I3_POINTER_TYPEDEFS(CCMRecoPulseSeriesMapUnion);

#endif
