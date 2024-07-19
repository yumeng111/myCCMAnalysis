/**
 *  $Id$
 *  
 *  Copyright (C) 2003-2007
 *  The IceCube Collaboration <http://www.icecube.wisc.edu>
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *  1. Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 *  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 *  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 *  OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 *  SUCH DAMAGE.
 *  
 *  SPDX-License-Identifier: BSD-2-Clause
 *  
 */

#ifndef CCMPMTKEY_H_INCLUDED
#define CCMPMTKEY_H_INCLUDED

#include <utility>
#include "Utility.h"
#include <iostream>
#include <icetray/IcetrayFwd.h>
#include <icetray/serialization.h>

static const unsigned ccmpmtkey_version_ = 0;
#pragma GCC diagnostic push
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

/**
 * @brief A small class which is the region number, om number
 * and pmt number for a specific PMT inside a DOM.
 *
 * For IceCube, the PMT number will always be 0
 * and "PMT" is equivalent to "DOM". For IceTop, the PMT number
 * can be 0 or 1.
 */

class CCMPMTKey
{
  int regionNumber_;
  unsigned int sensorNumber_;
  unsigned char subsensorNumber_;

 public:

  CCMPMTKey() : regionNumber_(0), sensorNumber_(0), subsensorNumber_(0) {}

  CCMPMTKey(int region,unsigned int sensor) 
    : regionNumber_(region), sensorNumber_(sensor), subsensorNumber_(0) {}

  CCMPMTKey(int region,unsigned int sensor, unsigned char subsensor) 
    : regionNumber_(region), sensorNumber_(sensor), subsensorNumber_(subsensor) {}

  virtual ~CCMPMTKey(); 

  /**
   * retrieves the region number for this CCMPMTKey
   */
  int GetRegion() const { return regionNumber_; }

  /**
   * Sets the region number for this OM
   */
  void SetRegion(int region){ regionNumber_ = region; }

  /**
   * gets the sensor number
   */
  unsigned int GetSensor() const { return sensorNumber_; }

  /**
   * sets the sensor number in this region
   */
  void SetSensor(unsigned int sensor){ sensorNumber_ = sensor; }

  /**
   * gets the PMT number in the DOM
   */
  unsigned char GetSubsensor() const { return subsensorNumber_; }
    
  /**
   * sets the PMT number in the DOM
   */
  void SetSubsensor(unsigned char subsensor){ subsensorNumber_ = subsensor; }

  /**
   * equality operator.  
   * @return true if the string and om numbers of the two CCMPMTKey's match
   * @param rhs the CCMPMTKey to compare this one to.
   */
  bool operator==(const CCMPMTKey& rhs) const
    {
      return (rhs.sensorNumber_ == sensorNumber_ && 
         rhs.regionNumber_ == regionNumber_ && 
         rhs.subsensorNumber_ == subsensorNumber_);
    }

  /**
   * inequality operator
   * @return false if the string or om numbers are different
   * @param rhs the CCMPMTKey to compare this one to.
   */
  bool operator!=(const CCMPMTKey& rhs) const
    {
      return not (rhs == *this);
    }

  std::string str() const;

  struct hash
  {
    size_t operator()(const CCMPMTKey& key) const
    {
      return (((static_cast<uint64_t>(abs(key.GetRegion()+256)) * 256) + 
               static_cast<uint64_t>(key.GetSensor()))) * 256 + 
              static_cast<uint64_t>(key.GetSubsensor());
    }
  };

 private:
  friend class icecube::serialization::access;

  template <class Archive>
  void save(Archive& ar, unsigned version) const;

  template <class Archive>
  void load(Archive& ar, unsigned version);

  I3_SERIALIZATION_SPLIT_MEMBER()
};

I3_CLASS_VERSION(CCMPMTKey, ccmpmtkey_version_);

/**
 * comparison operator.  First compares the string numbers, then compares
 * the om numbers.  Required to put CCMPMTKeys as the key of a map
 * @param lhs the left-hand CCMPMTKey
 * @param rhs the right-hand CCMPMTKey
 * @return true if the lhs should be ordered before the rhs
 */
inline bool operator<(const CCMPMTKey& lhs,const CCMPMTKey& rhs)
{
  if(lhs.GetRegion() < rhs.GetRegion()) {
    return true;
  } else if(lhs.GetRegion() > rhs.GetRegion()) {
    return false;
  } else if(lhs.GetSensor() < rhs.GetSensor()) {
    return true;
  } else if(lhs.GetSensor() > rhs.GetSensor()) {
    return false;
  } else if(lhs.GetSubsensor() < rhs.GetSubsensor()) {
    return true;
  } else {
    return false;
  }
}

inline bool operator>(const CCMPMTKey& lhs,const CCMPMTKey& rhs)
{
  if(lhs.GetRegion() > rhs.GetRegion()) {
    return true;
  } else if(lhs.GetRegion() < rhs.GetRegion()) {
    return false;
  } else if(lhs.GetSensor() > rhs.GetSensor()) {
    return true;
  } else if(lhs.GetSensor() < rhs.GetSensor()) {
    return false;
  } else if(lhs.GetSubsensor() > rhs.GetSubsensor()) {
    return true;
  } else {
    return false;
  }
}

/**
 * streams an CCMPMTKey to an arbitrary ostream.  These are important,
 * the tray uses these conversion internally.
 */
std::ostream& operator<<(std::ostream&, const CCMPMTKey& key);
std::istream& operator>>(std::istream&,  CCMPMTKey&);
#pragma GCC diagnostic pop

I3_POINTER_TYPEDEFS(CCMPMTKey);

#endif //CCMPMTKEY_H_INCLUDED
