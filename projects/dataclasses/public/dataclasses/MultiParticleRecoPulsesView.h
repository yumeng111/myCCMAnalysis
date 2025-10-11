/**
 *  $Id$
 *
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *
 */

#ifndef DATACLASSES_MultiParticleRecoPulsesView_H_INCLUDED
#define DATACLASSES_MultiParticleRecoPulsesView_H_INCLUDED

#include <string>
#include "icetray/I3FrameObject.h"
#include "icetray/I3Frame.h"
#include "icetray/serialization.h"
#include "dataclasses/physics/CCMRecoPulse.h"

static const unsigned multiparticlerecopulseview_version_ = 0;

class MultiParticleRecoPulsesView : public I3FrameObject {
public:
	/*
	 * Construct a union of pulses across particles for the map stored at "key.".
	 */
    MultiParticleRecoPulsesView(I3Frame const &, std::string const & key);
    MultiParticleRecoPulsesView();

	std::ostream& Print(std::ostream&) const override;

	CCMRecoPulseSeriesMapConstPtr Apply(const I3Frame&) const;
    std::string GetSource() const { return key_; }

	bool operator==(const MultiParticleRecoPulsesView&) const;
	bool operator!=(const MultiParticleRecoPulsesView&) const;
private:
	std::string key_;
	mutable CCMRecoPulseSeriesMapPtr unified_;

	friend class icecube::serialization::access;
	template <class Archive> void serialize(Archive& ar, unsigned version);

	SET_LOGGER("MultiParticleRecoPulsesView");
};

std::ostream& operator<<(std::ostream&, const MultiParticleRecoPulsesView&);

I3_CLASS_VERSION(MultiParticleRecoPulsesView, multiparticlerecopulseview_version_);
I3_POINTER_TYPEDEFS(MultiParticleRecoPulsesView);

#endif // DATACLASSES_MultiParticleRecoPulsesView_H_INCLUDED
