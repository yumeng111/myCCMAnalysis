/**
    copyright  (C) 2004
    the icecube collaboration
    @version $Id$
    @date    $Date$
*/

#ifndef DATACLASSES_I3MAP_H_INCLUDED
#define DATACLASSES_I3MAP_H_INCLUDED

#include <icetray/serialization.h>
#include <map>
#include <string>
#include <vector>

#include <dataclasses/Utility.h>
#include <icetray/I3Logging.h>
#include <icetray/I3FrameObject.h>
#include <icetray/has_operator.h>
#include <icetray/OMKey.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMPMTKeyPair.h>
#include "dataclasses/TriggerKey.h"
#include "dataclasses/ostream_overloads.hpp"

#ifndef __CINT__  // it puts the lotion in the basket
#include <boost/lexical_cast.hpp>
#endif

template <typename Key, typename Value>
struct I3Map : public I3FrameObject, public std::map<Key, Value>
{
  typedef std::map<Key, Value> base_t;
  typedef typename base_t::value_type value_type;
  typedef typename base_t::key_compare Compare;
  typedef typename base_t::allocator_type Allocator;

  I3Map() { }

  explicit I3Map( const Compare& comp,
              const Allocator& alloc = Allocator() ) :
      base_t(comp, alloc) {}

  explicit I3Map( const Allocator& alloc ) :
      base_t(alloc) {}

  template< class InputIt >
  I3Map( InputIt first, InputIt last,
     const Compare& comp = Compare(),
     const Allocator& alloc = Allocator() ) :
      base_t(first, last, comp, alloc) {}

  template< class InputIt >
  I3Map( InputIt first, InputIt last,
     const Allocator& alloc ) :
      base_t(first, last, alloc) {}

  I3Map( const I3Map<Key, Value>& other ) :
    base_t((base_t)other) {}
  I3Map( const base_t& other ) :
      base_t(other) {}

  I3Map( const I3Map<Key, Value>& other, const Allocator& alloc ) :
      base_t((base_t)other, alloc) {}
  I3Map( const base_t& other, const Allocator& alloc ) :
      base_t(other, alloc) {}

  I3Map( I3Map&& other ) :
      base_t((base_t)other) {}
  I3Map( base_t&& other ) :
      base_t(other) {}

  I3Map( I3Map&& other, const Allocator& alloc ) :
      base_t((base_t)other) {}
  I3Map( base_t&& other, const Allocator& alloc ) :
      base_t(other) {}

  I3Map( std::initializer_list<value_type> init,
     const Compare& comp = Compare(),
     const Allocator& alloc = Allocator() ) :
      base_t(init, comp, alloc) {}

  I3Map( std::initializer_list<value_type> init,
     const Allocator& alloc) :
      base_t(init, alloc) {}

  I3Map& operator=( const I3Map& other ) = default;
  I3Map& operator=( const base_t& other ) {
      return base_t::operator=(other);
  }

  I3Map& operator=( I3Map&& other ) = default;
  I3Map& operator=( base_t&& other ) {
      return base_t::operator=(other);
  }

  I3Map& operator=( std::initializer_list<value_type> ilist ) {
      return base_t::operator=(ilist);
  }

  template <class Archive>
  void serialize(Archive & ar, unsigned version)
  {
    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("map", base_object< std::map<Key, Value> >(*this));
  }

  ~I3Map();

  const Value&
  at(const Key& where) const
  {
    typename std::map<Key, Value>::const_iterator iter = this->find(where);
    if (iter == this->end())
      log_fatal("Map contains nothing at %s.", boost::lexical_cast<std::string>(where).c_str());

    return iter->second;
  }

  Value&
  at(const Key& where)
  {
    typename std::map<Key, Value>::iterator iter = this->find(where);
    if (iter == this->end())
      log_fatal("Map contains nothing at %s.", boost::lexical_cast<std::string>(where).c_str());

    return iter->second;
  }

  std::ostream& Print(std::ostream& os) const override{
    constexpr bool can_print=has_operator::insertion<std::ostream&,Key>::value
                             && has_operator::insertion<std::ostream&,Value>::value;
    return(PrintImpl(os,std::integral_constant<bool,can_print>()));
  }

private:
  std::ostream& PrintImpl(std::ostream& os, std::true_type can_print) const{
    os << '[';
    bool first=true;
    for(const auto& entry : *this){
      if(first)
        first=false;
      else
        os << ",\n";
      os << entry.first << " => " << entry.second;
    }
    os << ']';
    return os;
  }

  std::ostream& PrintImpl(std::ostream& os, std::false_type cant_print) const{
    os << '[' << this->size() << " element" << (this->size()==1?"":"s") << ']';
    return os;
  }
};

template <typename Key, typename Value>
I3Map<Key, Value> :: ~I3Map() { }

template <typename Key, typename Value>
std::ostream& operator<<(std::ostream& os, const I3Map<Key, Value> m){
  return(m.Print(os));
}

typedef I3Map<double, std::vector<double>> I3MapDoubleVectorDouble;
typedef I3Map<std::string, double> I3MapStringDouble;
typedef I3Map<std::string, int> I3MapStringInt;
typedef I3Map<std::string, bool> I3MapStringBool;
typedef I3Map<std::string, std::string> I3MapStringString;

typedef I3Map<std::string, std::vector<double> > I3MapStringVectorDouble;


typedef I3Map<std::string, I3MapStringDouble> I3MapStringStringDouble;

typedef I3Map<unsigned, unsigned> I3MapUnsignedUnsigned;
typedef I3Map<unsigned short, unsigned short> I3MapUShortUShort;
typedef I3Map<unsigned int, int> I3MapUIntInt;

typedef I3Map<unsigned int, std::vector<double> > I3MapUIntVectorDouble;
typedef I3Map<unsigned int, std::vector<unsigned int> > I3MapUIntVectorUInt;
typedef I3Map<int, std::vector<int> > I3MapIntVectorInt;
typedef I3Map<int, size_t> I3MapIntSizeT;
typedef I3Map<OMKey, std::vector<double> > I3MapKeyVectorDouble;
typedef I3Map<OMKey, std::vector<int> > I3MapKeyVectorInt;
typedef I3Map<OMKey, double > I3MapKeyDouble;
typedef I3Map<OMKey, unsigned int > I3MapKeyUInt;
typedef I3Map<TriggerKey, std::vector<unsigned int> > I3MapTriggerVectorUInt;
typedef I3Map<TriggerKey, double > I3MapTriggerDouble;
typedef I3Map<TriggerKey, int> I3MapTriggerInt;
typedef I3Map<TriggerKey, unsigned int> I3MapTriggerUInt;
typedef I3Map<CCMPMTKey, std::vector<std::vector<double>>> I3MapPMTKeyVectorVectorDouble;
typedef I3Map<CCMPMTKey, std::vector<std::vector<int>>> I3MapPMTKeyVectorVectorInt;
typedef I3Map<CCMPMTKey, std::vector<std::vector<unsigned int>>> I3MapPMTKeyVectorVectorUInt;
typedef I3Map<CCMPMTKey, std::vector<double> > I3MapPMTKeyVectorDouble;
typedef I3Map<CCMPMTKey, std::vector<int> > I3MapPMTKeyVectorInt;
typedef I3Map<CCMPMTKey, std::vector<unsigned int> > I3MapPMTKeyVectorUInt;
typedef I3Map<CCMPMTKey, double > I3MapPMTKeyDouble;
typedef I3Map<CCMPMTKey, int > I3MapPMTKeyInt;
typedef I3Map<CCMPMTKey, unsigned int > I3MapPMTKeyUInt;
typedef I3Map<CCMPMTKey, bool> I3MapPMTKeyBool;
typedef I3Map<CCMPMTKey, std::pair<uint32_t, uint32_t> > I3MapPMTKeyPairUInt32UInt32;
typedef I3Map<CCMPMTKey, std::pair<uint64_t, uint64_t> > I3MapPMTKeyPairUInt64UInt64;
typedef I3Map<std::tuple<CCMPMTKey, CCMPMTKey>, double> I3MapTuplePMTKeyPMTKeyDouble;
typedef I3Map<std::tuple<CCMPMTKey, CCMPMTKey>, std::vector<double>> I3MapTuplePMTKeyPMTKeyVectorDouble;
typedef I3Map<CCMPMTKeyPair, double> I3MapPMTKeyPairDouble;
typedef I3Map<CCMPMTKeyPair, std::vector<double>> I3MapPMTKeyPairVectorDouble;

I3_POINTER_TYPEDEFS(I3MapDoubleVectorDouble);
I3_POINTER_TYPEDEFS(I3MapStringDouble);
I3_POINTER_TYPEDEFS(I3MapStringInt);
I3_POINTER_TYPEDEFS(I3MapStringBool);
I3_POINTER_TYPEDEFS(I3MapStringString);
I3_POINTER_TYPEDEFS(I3MapStringVectorDouble);
I3_POINTER_TYPEDEFS(I3MapStringStringDouble);
I3_POINTER_TYPEDEFS(I3MapUnsignedUnsigned);
I3_POINTER_TYPEDEFS(I3MapUShortUShort);
I3_POINTER_TYPEDEFS(I3MapUIntInt);
I3_POINTER_TYPEDEFS(I3MapUIntVectorDouble);
I3_POINTER_TYPEDEFS(I3MapUIntVectorUInt);
I3_POINTER_TYPEDEFS(I3MapIntVectorInt);
I3_POINTER_TYPEDEFS(I3MapIntSizeT);
I3_POINTER_TYPEDEFS(I3MapKeyVectorDouble);
I3_POINTER_TYPEDEFS(I3MapKeyVectorInt);
I3_POINTER_TYPEDEFS(I3MapKeyDouble);
I3_POINTER_TYPEDEFS(I3MapKeyUInt);
I3_POINTER_TYPEDEFS(I3MapTriggerVectorUInt);
I3_POINTER_TYPEDEFS(I3MapTriggerDouble);
I3_POINTER_TYPEDEFS(I3MapTriggerInt);
I3_POINTER_TYPEDEFS(I3MapTriggerUInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyVectorVectorDouble);
I3_POINTER_TYPEDEFS(I3MapPMTKeyVectorVectorInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyVectorVectorUInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyVectorDouble);
I3_POINTER_TYPEDEFS(I3MapPMTKeyVectorInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyVectorUInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyDouble);
I3_POINTER_TYPEDEFS(I3MapPMTKeyInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyUInt);
I3_POINTER_TYPEDEFS(I3MapPMTKeyBool);
I3_POINTER_TYPEDEFS(I3MapPMTKeyPairUInt32UInt32);
I3_POINTER_TYPEDEFS(I3MapPMTKeyPairUInt64UInt64);
I3_POINTER_TYPEDEFS(I3MapTuplePMTKeyPMTKeyDouble);
I3_POINTER_TYPEDEFS(I3MapTuplePMTKeyPMTKeyVectorDouble);
I3_POINTER_TYPEDEFS(I3MapPMTKeyPairDouble);
I3_POINTER_TYPEDEFS(I3MapPMTKeyPairVectorDouble);

#endif // I3MAP_H_INCLUDED

