/**
 *  $Id$
 *  
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */

#ifndef DATACLASSES_I3MAPCCMPMTKEYMASK_H_INCLUDED
#define DATACLASSES_I3MAPCCMPMTKEYMASK_H_INCLUDED

#include <functional>
#include <string>
#include <list>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <boost/dynamic_bitset.hpp> 

#include "icetray/I3FrameObject.h"
#include "icetray/CCMPMTKey.h"
#include "icetray/I3Frame.h"
#include "icetray/serialization.h"
#include "dataclasses/physics/CCMRecoPulse.h"

static const unsigned ccmrecopulseseriesmapmask_version_ = 1;

class CCMRecoPulseSeriesMapMask : public I3FrameObject {
public:
	/* 
	 * Construct a mask for the map stored at "key." All bits are set.
	 */
	CCMRecoPulseSeriesMapMask(const I3Frame&, const std::string &key);
	/*
	 * Construct a mask for the map stored at "key," corresponding to
	 * the map `subset.' The subset must be an ordered subset of the map
	 * stored at `key'.
	 */
	CCMRecoPulseSeriesMapMask(const I3Frame&, const std::string &key,
	    const CCMRecoPulseSeriesMap &subset);
	/*
	 * Construct a mask by predicate. This may be either a free function
	 * or an instance of a class that implements
	 * bool operator()(const CCMPMTKey&, size_t, const CCMRecoPulse&).
	 */
 	CCMRecoPulseSeriesMapMask(const I3Frame&, const std::string &key,
 	    boost::function<bool (const CCMPMTKey&, size_t, const CCMRecoPulse&)> predicate);
	
	CCMRecoPulseSeriesMapMask();
	 
	std::ostream& Print(std::ostream&) const override;
  
	/*
	 * Set/unset all elements for an CCMPMTKey.
	 */
	void Set(const CCMPMTKey&, bool);
	/*
	 * Set/unset an element for CCMPMTKey by index.
	 */
	void Set(const CCMPMTKey&, const unsigned idx, bool);
	/*
	 * Set/unset an element for CCMPMTKey by value.
	 */
	void Set(const CCMPMTKey&, const CCMRecoPulse&, bool);
	/*
	 * Clear all elements of the mask.
	 */
	void SetNone();
	
	/*
	 * Apply the mask to the target map in the frame.
	 */
	boost::shared_ptr<const CCMRecoPulseSeriesMap> Apply(const I3Frame &frame) const;
	
	/**
	 * Return true if this mask is derived from key
	 */
	bool HasAncestor(const I3Frame &frame, const std::string &key) const;
	
	/**
	 * Convert this mask into a form that can be applied directly to key
	 * without resort to intermediate masks.
	 */
	boost::shared_ptr<CCMRecoPulseSeriesMapMask> Repoint(const I3Frame &frame, const std::string &key) const;
	
	/*
	 * Get the name of the frame object the mask was made from.
	 */
	std::string GetSource() const { return key_; }
	 
	/*
	 * Get the number of set bits in the mask.
	 */
	unsigned GetSum() const;
	
	/*
	 * Are any bits set?
	 */
	bool GetAnySet() const;
	
	/*
	 * Are all bits set?
	 */
	bool GetAllSet() const;
	 
	 
	std::vector<boost::dynamic_bitset<uint8_t> > GetBits() const;
	
	/*
	 * Logical operators, applied elementwise.
	 */
	CCMRecoPulseSeriesMapMask operator&(const CCMRecoPulseSeriesMapMask&) const;
	CCMRecoPulseSeriesMapMask operator|(const CCMRecoPulseSeriesMapMask&) const;
	CCMRecoPulseSeriesMapMask operator^(const CCMRecoPulseSeriesMapMask&) const;
	/** Equivalent to this & ~other */
	CCMRecoPulseSeriesMapMask Remove(const CCMRecoPulseSeriesMapMask&) const;
	
	CCMRecoPulseSeriesMapMask& operator&=(const CCMRecoPulseSeriesMapMask&);
	CCMRecoPulseSeriesMapMask& operator|=(const CCMRecoPulseSeriesMapMask&);
	CCMRecoPulseSeriesMapMask& operator^=(const CCMRecoPulseSeriesMapMask&);
		
	bool operator==(const CCMRecoPulseSeriesMapMask&) const;
	bool operator!=(const CCMRecoPulseSeriesMapMask&) const;
	
private:
	typedef uint8_t mask_t;
	
	struct bitmask {
		uint16_t size_;
		uint8_t padding_;
		mask_t *mask_;
		
		bitmask() : size_(0), padding_(0), mask_(NULL) {};
		bitmask(unsigned length, bool set=true);
		bitmask(const bitmask& other);
		bitmask& operator=(const bitmask& other);
		~bitmask();
		void set_all();
		void unset_all();
		inline bool any() const;
		inline bool all() const;
		inline void set(const unsigned, bool);
		
		inline bool get(const unsigned) const;
		unsigned sum() const;
		size_t size() const;
		
		bool operator==(const bitmask&) const;
		
		friend class icecube::serialization::access;
		
		template <class Archive> void load(Archive & ar, unsigned version);
		template <class Archive> void save(Archive & ar, unsigned version) const;
		
		I3_SERIALIZATION_SPLIT_MEMBER();
	};
	
	friend std::ostream& operator<<(std::ostream&, const bitmask&);
		
	std::string key_;
	bitmask omkey_mask_;
	std::list<bitmask> element_masks_;
	CCMRecoPulseSeriesMapConstPtr source_;
	mutable CCMRecoPulseSeriesMapPtr masked_;
	
	inline void ResetCache() { masked_.reset(); }

	int FindKey(const CCMPMTKey &key, std::list<bitmask>::iterator &list_it,
	    const CCMRecoPulseSeriesMap::mapped_type **vec);
	
	static bool IsOrderedSubset(const CCMRecoPulseSeriesMap&, const CCMRecoPulseSeriesMap&);
	static void FillSubsetMask(bitmask&, const CCMRecoPulseSeriesMap::mapped_type&,
	    const CCMRecoPulseSeriesMap::mapped_type&);
	
	/**
	 * Collapse this mask with its source, making it depend only on its grandparent.
	 */
	boost::shared_ptr<CCMRecoPulseSeriesMapMask> CollapseLevel(const I3Frame &frame) const;
	
	template <typename BinaryOperator>
	CCMRecoPulseSeriesMapMask ApplyBinaryOperator(const CCMRecoPulseSeriesMapMask&) const;
	
	struct operator_and : public std::binary_function<mask_t, mask_t, mask_t> {
		inline mask_t operator()(mask_t lhs, mask_t rhs) { return lhs & rhs; }
	};
	
	struct operator_andnot : public std::binary_function<mask_t, mask_t, mask_t> {
		inline mask_t operator()(mask_t lhs, mask_t rhs) { return lhs & ~rhs; }
	};
	
	struct operator_or : public std::binary_function<mask_t, mask_t, mask_t> {
		inline mask_t operator()(mask_t lhs, mask_t rhs) { return lhs | rhs; }
	};
	
	struct operator_xor : public std::binary_function<mask_t, mask_t, mask_t> {
		inline mask_t operator()(mask_t lhs, mask_t rhs) { return lhs ^ rhs; }
	};
	
	friend class icecube::serialization::access;
	template <class Archive> void load(Archive & ar, unsigned version);
	template <class Archive> void save(Archive & ar, unsigned version) const;
	
	I3_SERIALIZATION_SPLIT_MEMBER();
	
	SET_LOGGER("CCMRecoPulseSeriesMapMask");
};

std::ostream& operator<<(std::ostream&, const CCMRecoPulseSeriesMapMask&);

template<> void CCMRecoPulseSeriesMapMask::bitmask::load(icecube::archive::xml_iarchive& ar, unsigned version);
template<> void CCMRecoPulseSeriesMapMask::bitmask::save(icecube::archive::xml_oarchive& ar, unsigned version) const;

I3_CLASS_VERSION(CCMRecoPulseSeriesMapMask, ccmrecopulseseriesmapmask_version_);
I3_POINTER_TYPEDEFS(CCMRecoPulseSeriesMapMask);

#endif /* DATACLASSES_I3MAPOMKEYMASK_H_INCLUDED */

