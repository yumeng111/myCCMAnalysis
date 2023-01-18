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

#ifndef CCMTriggerKey_H_INCLUDED
#define CCMTriggerKey_H_INCLUDED

#include <utility>
#include "Utility.h"
#include <iostream>
#include <icetray/IcetrayFwd.h>
#include <icetray/serialization.h>

static const unsigned ccmtriggerkey_version_ = 2;
#pragma GCC diagnostic push
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

/**
 * @brief A small class which is the string number, om number
 * and pmt number for a specific PMT inside a DOM.
 *
 * For IceCube, the PMT number will always be 0
 * and "PMT" is equivalent to "DOM". For IceTop, the PMT number
 * can be 0 or 1.
 */

#define CCMTriggerKey_H_CCMTriggerKey_TriggerType        \
  (UnknownType)(BoardTriggerCopy)(BeamTrigger)(StrobeTrigger)(CosmicTrigger)(LaserTrigger)(LEDTopTrigger)(LEDBottomTrigger)

class CCMTriggerKey {
public:
    enum TriggerType {UnknownType = 0, BoardTriggerCopy = 10, BeamTrigger = 20, StrobeTrigger = 30, CosmicTrigger = 40, LaserTrigger = 50, LEDTopTrigger = 60, LEDBottomTrigger = 70};
private:
    TriggerType triggerType_;
    unsigned int triggerNumber_;

    public:

    CCMTriggerKey() : triggerType_(UnknownType), triggerNumber_(0) {}

    CCMTriggerKey(TriggerType type, unsigned int number)
        : triggerType_(type), triggerNumber_(number) {}

    virtual ~CCMTriggerKey(); 

    /**
     * retrieves the string number for this CCMTriggerKey
     */
    TriggerType GetType() const { return triggerType_; }

    /**
     * Sets the string number for this OM
     */
    void SetType(TriggerType type){ triggerType_ = type; }

    /**
     * gets the OM number on the string
     */
    unsigned int GetNumber() const { return triggerNumber_; }

    /**
     * sets the OM number on the string
     */
    void SetNumber(unsigned int number){ triggerNumber_ = number; }

    /**
     * equality operator.  
     * @return true if the string and om numbers of the two CCMTriggerKey's match
     * @param rhs the CCMTriggerKey to compare this one to.
     */
    bool operator==(const CCMTriggerKey& rhs) const
    {
        return (rhs.triggerType_ == triggerType_ &&
                rhs.triggerNumber_ == triggerNumber_);
    }

    /**
     * inequality operator
     * @return false if the string or om numbers are different
     * @param rhs the CCMTriggerKey to compare this one to.
     */
    bool operator!=(const CCMTriggerKey& rhs) const
    {
        return not (rhs == *this);
    }

    std::string str() const;

    struct hash
    {
        size_t operator()(const CCMTriggerKey& key) const
        {
            return (((static_cast<size_t>(abs(static_cast<int>(key.GetType()) + 256)) * 1024) + 
                        static_cast<size_t>(key.GetNumber())));
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

I3_CLASS_VERSION(CCMTriggerKey, ccmtriggerkey_version_);

/**
 * comparison operator.  First compares the string numbers, then compares
 * the om numbers.  Required to put CCMTriggerKeys as the key of a map
 * @param lhs the left-hand CCMTriggerKey
 * @param rhs the right-hand CCMTriggerKey
 * @return true if the lhs should be ordered before the rhs
 */
inline bool operator<(const CCMTriggerKey& lhs,const CCMTriggerKey& rhs)
{
    if(static_cast<int>(lhs.GetType()) < static_cast<int>(rhs.GetType())) {
        return true;
    } else if(static_cast<int>(lhs.GetType()) > static_cast<int>(rhs.GetType())) {
        return false;
    } else if(lhs.GetNumber() < rhs.GetNumber()) {
        return true;
    } else {
        return false;
    }
}

/**
 * streams an CCMTriggerKey to an arbitrary ostream.  These are important,
 * the tray uses these conversion internally.
 */
std::ostream& operator<<(std::ostream&, const CCMTriggerKey& key);
std::istream& operator>>(std::istream&,  CCMTriggerKey&);
#pragma GCC diagnostic pop

I3_POINTER_TYPEDEFS(CCMTriggerKey);

#endif //CCMTriggerKey_H_INCLUDED
