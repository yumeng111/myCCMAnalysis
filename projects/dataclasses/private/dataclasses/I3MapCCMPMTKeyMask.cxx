/**
 *  $Id$
 *  
 *  Copyright (C) 2011
 *  Jakob van Santen <vansanten@wisc.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 */

#include <algorithm>
#include "dataclasses/I3MapCCMPMTKeyMask.h"
#include "dataclasses/physics/CCMRecoPulse.h"
#include "boost/make_shared.hpp"
#include <serialization/binary_object.hpp>

CCMRecoPulseSeriesMapMask::CCMRecoPulseSeriesMapMask() {}

CCMRecoPulseSeriesMapMask::CCMRecoPulseSeriesMapMask(const I3Frame &frame, const std::string &key)
    : key_(key)
{
	source_ = frame.Get<boost::shared_ptr<const CCMRecoPulseSeriesMap> >(key_);

	if (!source_)
		log_fatal("The map named '%s' doesn't exist in the frame!\n", key_.c_str());

	CCMRecoPulseSeriesMap::const_iterator it = source_->begin();
	if (source_->size() == 0)
		return;

	omkey_mask_ = bitmask(source_->size());
	for ( ; it != source_->end(); it++) {
		if (it->second.size() > 0)
			element_masks_.push_back(bitmask(it->second.size()));
		else
			omkey_mask_.set(std::distance(source_->begin(), it),
			    false);
	}
}

void
CCMRecoPulseSeriesMapMask::FillSubsetMask(bitmask &mask,
    const CCMRecoPulseSeriesMap::mapped_type &superset,
    const CCMRecoPulseSeriesMap::mapped_type &subset)
{
	unsigned idx = 0;
	unsigned midx = 0;
	CCMRecoPulseSeriesMap::mapped_type::const_iterator sup_vit = superset.begin();
	CCMRecoPulseSeriesMap::mapped_type::const_iterator sub_vit = subset.begin();
	
	for ( ; sup_vit != superset.end(); sup_vit++, idx++) {
		if (idx > 0 && (idx % (8*sizeof(mask_t)) == 0))
			midx++;
		
		assert(midx < mask.size_);
		if ((sub_vit != subset.end()) && ((*sub_vit) == (*sup_vit))) {
			sub_vit++; /* We have a match, do nothing. */
		} else {
			/* Unset the corresponding bit. */
			mask.mask_[midx] &= ~(1 << (idx % (8*sizeof(mask_t))));
		}
	}
}

CCMRecoPulseSeriesMapMask::CCMRecoPulseSeriesMapMask(const I3Frame &frame,
    const std::string &key, const CCMRecoPulseSeriesMap &subset)
    : key_(key)
{
	source_ = frame.Get<boost::shared_ptr<const CCMRecoPulseSeriesMap> >(key_);
	
	if (!source_)
		log_fatal("The map named '%s' doesn't exist in the frame!", key_.c_str());
	
	if (!IsOrderedSubset(*source_, subset))
		log_fatal("The passed map is not an ordered subset of map '%s'!", key_.c_str());
	
	omkey_mask_ = bitmask(source_->size());
	assert(omkey_mask_.size() == source_->size());
	
	CCMRecoPulseSeriesMap::const_iterator sup_mit = source_->begin();
	CCMRecoPulseSeriesMap::const_iterator sub_mit = subset.begin();
	unsigned idx = 0;
	for ( ; sup_mit != source_->end(); sup_mit++, idx++) {
		/* 
		 * NB: since the subset is ordered, we can wait until the
		 * superset catches up to increment its iterator.
		 */
		if ((sub_mit != subset.end()) && (sub_mit->first == sup_mit->first)) {
			bitmask mask = bitmask(sup_mit->second.size(), sub_mit->second.size() != 0);
			
			FillSubsetMask(mask, sup_mit->second, sub_mit->second);
			element_masks_.push_back(mask);
			
			sub_mit++;
		} else {
			omkey_mask_.set(idx, false); /* Note absence. */
		}
	}
}

CCMRecoPulseSeriesMapMask::CCMRecoPulseSeriesMapMask(const I3Frame &frame,
    const std::string &key, boost::function<bool (const CCMPMTKey&, size_t, const CCMRecoPulse&)> predicate)
    : key_(key)
{
	source_ = frame.Get<boost::shared_ptr<const CCMRecoPulseSeriesMap> >(key_);
	
	if (!source_)
		log_fatal("The map named '%s' doesn't exist in the frame!\n", key_.c_str());
	
	omkey_mask_ = bitmask(source_->size(), false);
	assert(omkey_mask_.size() == source_->size());
	
	size_t om_idx(0);
	BOOST_FOREACH(const CCMRecoPulseSeriesMap::value_type &pair, *source_) {
		bitmask mask = bitmask(pair.second.size(), true);
		
		size_t pulse_idx(0);
		BOOST_FOREACH(const CCMRecoPulseSeriesMap::mapped_type::value_type &pulse, pair.second) {
			mask.set(pulse_idx, predicate(pair.first, pulse_idx, pulse));
			pulse_idx++;
		}
		
		if (mask.sum() > 0) {
			omkey_mask_.set(om_idx, true);
			element_masks_.push_back(mask);
		}
		om_idx++;
	}
}

void
CCMRecoPulseSeriesMapMask::SetNone()
{
	omkey_mask_.unset_all();
	element_masks_.clear();
	
	ResetCache();
}


unsigned
CCMRecoPulseSeriesMapMask::GetSum() const
{
	unsigned sum = 0;
	std::list<bitmask>::const_iterator list_it = element_masks_.begin();
	
	for ( ; list_it != element_masks_.end(); list_it++)
		sum += list_it->sum();
	
	return sum;
}

bool
CCMRecoPulseSeriesMapMask::GetAnySet() const
{
	return (omkey_mask_.any());
}

bool
CCMRecoPulseSeriesMapMask::GetAllSet() const
{
	if (omkey_mask_.all()) {
		bool all = true;
		std::list<bitmask>::const_iterator list_it = element_masks_.begin();
		for ( ; list_it != element_masks_.end(); list_it++)
			all = (all && list_it->all());
		
		return all;
	} else {
		return false;
	}
}

std::vector<boost::dynamic_bitset<uint8_t> >
CCMRecoPulseSeriesMapMask::GetBits() const
{
	typedef boost::dynamic_bitset<uint8_t> bitset;
	std::vector<bitset> bits;
	
	std::list<bitmask>::const_iterator list_it = element_masks_.begin();	
	for (unsigned i = 0; i < omkey_mask_.size(); i++)
		if (omkey_mask_.get(i)) {
			const bitmask &privmask = *(list_it++);
			bitset pubmask(privmask.size());
			for (unsigned j = 0; j < privmask.size(); j++)
				pubmask.set(j, privmask.get(j));
			bits.push_back(pubmask);
		} else {
			bits.push_back(bitset(0));
		}
	
	return bits;
}

int
CCMRecoPulseSeriesMapMask::FindKey(const CCMPMTKey &key,
    std::list<bitmask>::iterator &list_it, const CCMRecoPulseSeriesMap::mapped_type **vec)
{
	if (!source_) {
		log_error("This mask cannot be modified.");
		return -1;
	}
	
	CCMRecoPulseSeriesMap::const_iterator source_it = source_->find(key);
	
	if (source_it == source_->end())
		log_fatal("Key doesn't exist in the source map!");
		
	*vec = &(source_it->second);
	
	/* Find the bitmask corresponding to this OM. */
	list_it = element_masks_.begin();
	unsigned omkey_idx = std::distance(source_->begin(), source_it);
	for (unsigned i = 0; i < omkey_idx; i++)
		if (omkey_mask_.get(i))
			list_it++;
	
	return omkey_idx;
}

void
CCMRecoPulseSeriesMapMask::Set(const CCMPMTKey &key,
    const CCMRecoPulseSeriesMap::mapped_type::value_type &target, bool set_it)
{
	int omkey_idx;
	std::list<bitmask>::iterator list_it;
	const CCMRecoPulseSeriesMap::mapped_type *vec;
	
	ResetCache();
	
	if ((omkey_idx = FindKey(key, list_it, &vec)) < 0)
		return;
	
	/* Insert a new bitmask if necessary. */
	if (!omkey_mask_.get(omkey_idx)) {
		if (set_it) {
			list_it = element_masks_.insert(list_it,
			     bitmask(vec->size(), false));
			omkey_mask_.set(omkey_idx, true);
		} else {
			return;
		}
	}
		
	CCMRecoPulseSeriesMap::mapped_type::const_iterator vec_it = vec->begin();
	unsigned idx = 0;
	for ( ; vec_it != vec->end(); vec_it++, idx++)
		if (*vec_it == target) {
			list_it->set(idx, set_it);
			break;
		}	
}


void
CCMRecoPulseSeriesMapMask::Set(const CCMPMTKey &key, const unsigned idx, bool set_it)
{
	int omkey_idx;
	std::list<bitmask>::iterator list_it;
	const CCMRecoPulseSeriesMap::mapped_type *vec;
	
	ResetCache();
	
	if ((omkey_idx = FindKey(key, list_it, &vec)) < 0)
		return;
	
	/* Insert a new bitmask if necessary. */
	if (!omkey_mask_.get(omkey_idx)) {
		if (set_it) {
			list_it = element_masks_.insert(list_it,
			     bitmask(vec->size(), false));
			omkey_mask_.set(omkey_idx, true);
		} else {
			return;
		}
	}
			
	list_it->set(idx, set_it);
}

void
CCMRecoPulseSeriesMapMask::Set(const CCMPMTKey &key, bool set_it)
{
	int omkey_idx;
	std::list<bitmask>::iterator list_it;
	const CCMRecoPulseSeriesMap::mapped_type *vec;
	
	ResetCache();
	
	if ((omkey_idx = FindKey(key, list_it, &vec)) < 0)
		return;
	
	/* Insert a new bitmask if necessary. */
	if (!omkey_mask_.get(omkey_idx)) {
		if (set_it && vec->size() > 0) {
			list_it = element_masks_.insert(list_it,
			     bitmask(vec->size(), true));
			omkey_mask_.set(omkey_idx, true);
		}
		return;
	}
	
	if (set_it && vec->size() > 0) {
		list_it->set_all();
	} else {
		omkey_mask_.set(omkey_idx, false);

		if (!set_it)
			element_masks_.erase(list_it);
	}
}

bool
CCMRecoPulseSeriesMapMask::IsOrderedSubset(const CCMRecoPulseSeriesMap &super,
    const CCMRecoPulseSeriesMap &sub)
{
	CCMRecoPulseSeriesMap::const_iterator sup_mit = super.begin();
	CCMRecoPulseSeriesMap::const_iterator sub_mit = sub.begin();
	for ( ; sub_mit != sub.end(); sub_mit++) {
		while (sup_mit != super.end() && sup_mit->first != sub_mit->first)
			sup_mit++;
		if (sup_mit == super.end())
			return false;
		
		/* 
		 * NB: even in the vector portion, we expect elements in the subset
		 * to appear in order relative to their order in the superset.
		 */
		CCMRecoPulseSeriesMap::mapped_type::const_iterator sup_vit = sup_mit->second.begin();
		CCMRecoPulseSeriesMap::mapped_type::const_iterator sub_vit = sub_mit->second.begin();
		for ( ; sub_vit != sub_mit->second.end(); sub_vit++) {
			while (sup_vit == sup_mit->second.end() && !(*sup_vit == *sub_vit))
				sup_vit++;
			if (sup_vit == sup_mit->second.end())
				return false;
		}
	}
	
	return true;
}

CCMRecoPulseSeriesMapMask
CCMRecoPulseSeriesMapMask::operator&(const CCMRecoPulseSeriesMapMask &other) const
{
	return ApplyBinaryOperator<operator_and>(other);
}

CCMRecoPulseSeriesMapMask
CCMRecoPulseSeriesMapMask::operator|(const CCMRecoPulseSeriesMapMask &other) const
{
	return ApplyBinaryOperator<operator_or>(other);
}

CCMRecoPulseSeriesMapMask
CCMRecoPulseSeriesMapMask::operator^(const CCMRecoPulseSeriesMapMask &other) const
{
	return ApplyBinaryOperator<operator_xor>(other);
}

CCMRecoPulseSeriesMapMask
CCMRecoPulseSeriesMapMask::Remove(const CCMRecoPulseSeriesMapMask &other) const
{
	return ApplyBinaryOperator<operator_andnot>(other);
}

CCMRecoPulseSeriesMapMask&
CCMRecoPulseSeriesMapMask::operator&=(const CCMRecoPulseSeriesMapMask &other)
{
	(*this) = ApplyBinaryOperator<operator_and>(other);
	return *this;
}

CCMRecoPulseSeriesMapMask&
CCMRecoPulseSeriesMapMask::operator|=(const CCMRecoPulseSeriesMapMask &other)
{
	(*this) = ApplyBinaryOperator<operator_or>(other);
	return *this;
}

CCMRecoPulseSeriesMapMask&
CCMRecoPulseSeriesMapMask::operator^=(const CCMRecoPulseSeriesMapMask &other)
{
	(*this) = ApplyBinaryOperator<operator_xor>(other);
	return *this;
}

template <typename BinaryOperator>
CCMRecoPulseSeriesMapMask
CCMRecoPulseSeriesMapMask::ApplyBinaryOperator(const CCMRecoPulseSeriesMapMask &other) const
{
	if (key_ != other.key_)
		log_fatal("Masks cannot be combined!");
	
	CCMRecoPulseSeriesMapMask newmask;
	BinaryOperator op;
	
	newmask.key_ = key_;
	newmask.source_ = source_;
	newmask.omkey_mask_ = omkey_mask_;
	newmask.omkey_mask_.unset_all();
	
	const unsigned n_omkeys = omkey_mask_.size_*8*sizeof(mask_t) - omkey_mask_.padding_;
	std::list<bitmask>::const_iterator lit = element_masks_.begin();
	std::list<bitmask>::const_iterator rit = other.element_masks_.begin();
	
	for (unsigned omkey_idx = 0; omkey_idx < n_omkeys; omkey_idx++) {
		bitmask lhs, rhs;
		
		const bool l_exists = omkey_mask_.get(omkey_idx);
		const bool r_exists = other.omkey_mask_.get(omkey_idx);
		if (!(l_exists || r_exists)) {
			lhs = rhs = bitmask(8, false);
		} else if (l_exists && r_exists) {
			lhs = *lit;
			rhs = *rit;
			lit++; rit++;
		} else if (l_exists) {
			lhs = rhs = *lit;
			rhs.unset_all();
			lit++;
		} else {
			lhs = rhs = *rit;
			lhs.unset_all();
			rit++;
		}
		
		assert(lhs.size_ == rhs.size_);
		assert(lhs.padding_ == rhs.padding_);
		for (unsigned i = 0; i < lhs.size_; i++)
			lhs.mask_[i] = op(lhs.mask_[i], rhs.mask_[i]);
			
		if (lhs.sum() != 0) {
			newmask.element_masks_.push_back(lhs);
			newmask.omkey_mask_.set(omkey_idx, true);
		}		
	}
	
	return newmask;
}

bool
CCMRecoPulseSeriesMapMask::operator==(const CCMRecoPulseSeriesMapMask &other) const
{
	return (key_ == other.key_ &&
			omkey_mask_ == other.omkey_mask_ &&
			element_masks_ == other.element_masks_);
}

bool
CCMRecoPulseSeriesMapMask::operator!=(const CCMRecoPulseSeriesMapMask &other) const
{
	return !operator==(other);
}

boost::shared_ptr<const CCMRecoPulseSeriesMap>
CCMRecoPulseSeriesMapMask::Apply(const I3Frame &frame) const
{
	if (masked_)
		return masked_;
	
	boost::shared_ptr<const CCMRecoPulseSeriesMap> source =
	    frame.Get<boost::shared_ptr<const CCMRecoPulseSeriesMap> >(key_);
	
	if (!source)
		log_fatal("The map named '%s' doesn't exist in the frame!\n", key_.c_str());
	if (source->size() != omkey_mask_.size())
		log_fatal("This mask was made from a map with %zu keys, but "
		    "the map named '%s' has %zu keys.", omkey_mask_.size(),
		    key_.c_str(), source->size());
	
	masked_ = boost::make_shared<CCMRecoPulseSeriesMap>();
	
	CCMRecoPulseSeriesMap::const_iterator source_it = source->begin();
	CCMRecoPulseSeriesMap::iterator inserter = masked_->begin();
	std::list<bitmask>::const_iterator list_it = element_masks_.begin();
	unsigned omkey_idx = 0;
	
	for ( ; source_it != source->end(); source_it++, omkey_idx++) {
		
		if (!omkey_mask_.get(omkey_idx))
			continue;
		else if (list_it->sum() == 0) {
			list_it++;
			continue;
		}
		
		unsigned idx = 0;
		CCMRecoPulseSeriesMap::mapped_type target_vec;
		CCMRecoPulseSeriesMap::mapped_type::const_iterator source_vit =
		    source_it->second.begin();

		if (source_it->second.size() != list_it->size())
			log_fatal("The mask for OM(%d,%d) has %zu entries, but source "
			    "pulse vector has %zu entries!", source_it->first.GetRegion(),
			    source_it->first.GetSensor(), list_it->size(),
			    source_it->second.size());

		for ( ; source_vit != source_it->second.end(); source_vit++, idx++)
			if (list_it->get(idx))
				target_vec.push_back(*source_vit);
				
		inserter = masked_->insert(inserter,
		    std::make_pair(source_it->first, target_vec));
			
		list_it++;
	}
	
	return masked_;
}

struct null_deleter
{
	void operator()(void const *) const {}
};

bool 
CCMRecoPulseSeriesMapMask::HasAncestor(const I3Frame &frame, const std::string &key) const
{
	if (key == GetSource())
		return true;
	CCMRecoPulseSeriesMapMaskConstPtr source(this, null_deleter());
	while ((source = frame.Get<CCMRecoPulseSeriesMapMaskConstPtr>(source->GetSource()))) {
		if (source->GetSource() == key)
			return true;
	}
	return false;
}

boost::shared_ptr<CCMRecoPulseSeriesMapMask>
CCMRecoPulseSeriesMapMask::Repoint(const I3Frame &frame, const std::string &key) const
{
	if (!HasAncestor(frame, key))
		log_fatal_stream(key << " is not an ancestor of this mask");
	if (key == GetSource())
		return boost::make_shared<CCMRecoPulseSeriesMapMask>(*this);
	
	CCMRecoPulseSeriesMapMaskConstPtr source(this, null_deleter());
	while (source->GetSource() != key)
		source = source->CollapseLevel(frame);
	return boost::const_pointer_cast<CCMRecoPulseSeriesMapMask>(source);
}

boost::shared_ptr<CCMRecoPulseSeriesMapMask>
CCMRecoPulseSeriesMapMask::CollapseLevel(const I3Frame &frame) const
{
	boost::shared_ptr<const CCMRecoPulseSeriesMapMask> source =
	    frame.Get<boost::shared_ptr<const CCMRecoPulseSeriesMapMask> >(key_);
	if (!source)
		log_fatal_stream(key_ << " is not a mask in the frame");
	
	boost::shared_ptr<CCMRecoPulseSeriesMapMask> collapsed
	    = boost::make_shared<CCMRecoPulseSeriesMapMask>(*source);
	collapsed->ResetCache();
	
	std::list<bitmask>::const_iterator source_element = source->element_masks_.begin();
	std::list<bitmask>::const_iterator element = element_masks_.begin();
	std::list<bitmask>::iterator dest_element = collapsed->element_masks_.begin();
	
	// Perform an unaligned logical AND between the parent and daughter masks
	unsigned idx = 0;
	for (unsigned source_idx = 0; source_idx < source->omkey_mask_.size();
	    source_idx++) {
		if (source->omkey_mask_.get(source_idx)) {
			if (source_element->any()) {
				if (omkey_mask_.get(idx)) {
					unsigned pidx = 0;
					for (unsigned source_pidx = 0; source_pidx < source_element->size(); source_pidx++) {
						if (source_element->get(source_pidx)) {
							if (!element->get(pidx))
								dest_element->set(source_pidx, false);
							pidx++;
						}
					}
					element++;
					dest_element++;
				} else {
					collapsed->omkey_mask_.set(source_idx, false);
					dest_element = collapsed->element_masks_.erase(dest_element);
				}
				
				idx++;
			// If no bits are set in the parent mask, the daughter 
			// mask may not event exist. Compactify the output.
			} else {
				collapsed->omkey_mask_.set(source_idx, false);
				dest_element = collapsed->element_masks_.erase(dest_element);
			}
			
			source_element++;
		}
	}
	
	return collapsed;
}

CCMRecoPulseSeriesMapMask::bitmask::bitmask(unsigned length, bool set)
{
	size_ = (length != 0) ? (length-1u)/(8*sizeof(mask_t)) + 1 : 1;
	padding_ = size_*8*sizeof(mask_t) - length;
	mask_ = (mask_t*)malloc(size_*sizeof(mask_t));
	
	if (set)
		set_all();
	else
		unset_all();
};

CCMRecoPulseSeriesMapMask::bitmask::bitmask(const bitmask& other) : size_(other.size_), padding_(other.padding_), mask_(NULL)
{
	if (size_ != 0) {
		mask_ = (mask_t*)malloc(size_*sizeof(mask_t));
		memcpy(mask_, other.mask_, size_*sizeof(mask_t));
	}
}

CCMRecoPulseSeriesMapMask::bitmask&
CCMRecoPulseSeriesMapMask::bitmask::operator=(const bitmask& other)
{
	size_ = other.size_;
	padding_ = other.padding_;
	if (mask_)
		mask_ = (mask_t*)realloc(mask_, size_*sizeof(mask_t));
	else
		mask_ = (mask_t*)malloc(size_*sizeof(mask_t));
	memcpy(mask_, other.mask_, size_*sizeof(mask_t));
	
	return *this;
}

CCMRecoPulseSeriesMapMask::bitmask::~bitmask()
{
	if (mask_)
		free(mask_);
};

inline bool
CCMRecoPulseSeriesMapMask::bitmask::any() const
{
	mask_t test = 0;
	for (unsigned i = 0; i < size_; i++)
		test |= mask_[i];
	
	return (test != 0);
}

inline bool
CCMRecoPulseSeriesMapMask::bitmask::all() const
{
	bool test = true;
	const mask_t all = std::numeric_limits<mask_t>::max();
	
	const unsigned max_i = (size_ == 0) ? 0 : size_-1;
	for (unsigned i = 0; i < max_i; i++)
		test = test && (mask_[i] == all);
	
	return (test && ((mask_[size_-1] == (all >> padding_))));
}

void
CCMRecoPulseSeriesMapMask::bitmask::set_all()
{
	/* All bits set. */
	const mask_t all = std::numeric_limits<mask_t>::max();
	
	assert(size_ > 0);
	memset(mask_, all, (size_-1)*sizeof(mask_t));
	mask_[size_-1] = (all >> padding_); /* Unset unused top bits. */
}

void
CCMRecoPulseSeriesMapMask::bitmask::unset_all()
{
	memset(mask_, 0, size_*sizeof(mask_t));
}

inline void
CCMRecoPulseSeriesMapMask::bitmask::set(const unsigned idx, bool set_it)
{
	assert(mask_ && idx < (8*sizeof(mask_t)*size_ - padding_));
	if (set_it)
		mask_[idx/(8*sizeof(mask_t))] |= (1 << (idx % (8*sizeof(mask_t))));
	else
		mask_[idx/(8*sizeof(mask_t))] &= ~(1 << (idx % (8*sizeof(mask_t))));
}

inline bool
CCMRecoPulseSeriesMapMask::bitmask::get(const unsigned idx) const
{
	assert(mask_ && idx < (8*sizeof(mask_t)*size_ - padding_));
	
	return mask_[idx/(8*sizeof(mask_t))] & (1 << (idx % (8*sizeof(mask_t))));
}

unsigned
CCMRecoPulseSeriesMapMask::bitmask::sum() const
{
	unsigned sum = 0;
	
	for (unsigned i = 0; i < size_; i++)
		for (unsigned j = 0; j < 8*sizeof(mask_t); j++)
			if (mask_[i] & (1 << j))
				sum++;
	
	return sum;
}

size_t
CCMRecoPulseSeriesMapMask::bitmask::size() const
{
	return 8*sizeof(mask_t)*size_ - padding_;
}

bool
CCMRecoPulseSeriesMapMask::bitmask::operator==(const bitmask& other) const
{
	return (size_ == other.size_ &&
			padding_ == other.padding_ &&
			std::equal(mask_,mask_+sizeof(mask_t)*size_,other.mask_));
}

std::ostream& operator<<(std::ostream& os, const CCMRecoPulseSeriesMapMask::bitmask& mask){
	for(std::size_t i=0; i<mask.size(); i++)
		os << (mask.get(i) ? '1' : '0');
	return os;
}

std::ostream& CCMRecoPulseSeriesMapMask::Print(std::ostream& oss) const{
	oss << "CCMRecoPulseSeriesMapMask:\n"
	    << "  Key: " << key_ << '\n'
	    << "  CCMPMTKeyMask: " << omkey_mask_ << '\n'
	    << "  ElementMasks:";
	std::size_t idx=0;
	for(const bitmask& mask : element_masks_)
		oss << "\n    " << idx++ << ": " << mask;
	return oss;
}

std::ostream& operator<<(std::ostream& os, const CCMRecoPulseSeriesMapMask& mask){
	return(mask.Print(os));
}

template <class Archive>
void
CCMRecoPulseSeriesMapMask::bitmask::load(Archive & ar, unsigned version)
{
	ar & make_nvp("Size", size_);
	ar & make_nvp("Padding", padding_);
	
	mask_ = (mask_t*)malloc(size_*sizeof(mask_t));
	ar & make_nvp("Bitmask", icecube::serialization::make_binary_object(
	    mask_, size_*sizeof(mask_t)));
}

template <class Archive>
void
CCMRecoPulseSeriesMapMask::bitmask::save(Archive & ar, unsigned version) const
{
	ar & make_nvp("Size", size_);
	ar & make_nvp("Padding", padding_);
	ar & make_nvp("Bitmask", icecube::serialization::make_binary_object(
	    mask_, size_*sizeof(mask_t)));
}

template<>
void
CCMRecoPulseSeriesMapMask::bitmask::save(icecube::archive::xml_oarchive& ar, unsigned version) const
{
	ar & make_nvp("Size", size_);
	ar & make_nvp("Padding", padding_);
	// Generate a convenient ASCII representation of the mask bits
	std::string s(this->size(),'0');
	unsigned i = 0;
	for (std::string::iterator it = s.begin(); it != s.end(); it++, i++)
		if (this->get(i))
			*it = '1';
	ar & make_nvp("Bitmask", s);
}

template<>
void
CCMRecoPulseSeriesMapMask::bitmask::load(icecube::archive::xml_iarchive& ar, unsigned version)
{
	ar & make_nvp("Size", size_);
	ar & make_nvp("Padding", padding_);
	std::string s;
	ar & make_nvp("Bitmask", s);
	mask_ = (mask_t*)malloc(size_*sizeof(mask_t));
	memset(mask_, 0, size_*sizeof(mask_t));
	// Parse the ASCII representation of the mask bits, treating any non-'0' character as true
	unsigned i = 0;
	for (std::string::iterator it = s.begin(); it != s.end(); it++, i++)
		if (*it != '0')
			this->set(i, true);
}

template <class Archive>
void
CCMRecoPulseSeriesMapMask::load(Archive & ar, unsigned version)
{
	std::vector<bitmask> elements;
	float time_reference;
	
	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Key", key_);
	ar & make_nvp("CCMPMTKeyMask", omkey_mask_);
	ar & make_nvp("ElementMasks", elements);
	
	std::copy(elements.begin(), elements.end(), std::back_inserter(element_masks_));
}

template <class Archive>
void
CCMRecoPulseSeriesMapMask::save(Archive & ar, unsigned version) const
{
	/* Remove trivial elements before serializing */
	std::vector<bitmask> elements;
	bitmask omkey_mask = omkey_mask_;
	std::list<bitmask>::const_iterator list_it = element_masks_.begin();
	for (unsigned idx = 0; idx != omkey_mask_.size(); idx++) {
		if (omkey_mask_.get(idx)) {
			if (!list_it->any()) {
				omkey_mask.set(idx, false);
				list_it++;
			} else {
				elements.push_back(*list_it++);
			}
		}
	}

	ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
	ar & make_nvp("Key", key_);
	ar & make_nvp("CCMPMTKeyMask", omkey_mask);
	ar & make_nvp("ElementMasks", elements);
}

I3_SPLIT_SERIALIZABLE(CCMRecoPulseSeriesMapMask);
