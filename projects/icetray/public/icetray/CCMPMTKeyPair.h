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

#ifndef CCMPMTKEYPAIR_H_INCLUDED
#define CCMPMTKEYPAIR_H_INCLUDED

#include <utility>
#include "Utility.h"
#include <iostream>
#include <icetray/IcetrayFwd.h>
#include <icetray/serialization.h>
#include <icetray/CCMPMTKey.h>

static const unsigned ccmpmtkeypair_version_ = 0;
#pragma GCC diagnostic push
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

class CCMPMTKeyPair {
public:
    CCMPMTKey first;
    CCMPMTKey second;

    CCMPMTKeyPair() = default;

    CCMPMTKeyPair(CCMPMTKey first_, CCMPMTKey second_) : first(first_), second(second_) {}

    virtual ~CCMPMTKeyPair() = default;

    bool operator==(CCMPMTKeyPair const & rhs) const {
        return (rhs.first == first &&
                rhs.second == second);
    }

    bool operator!=(CCMPMTKeyPair const & rhs) const {
        return not (rhs == *this);
    }

    std::string str() const;

    struct hash {
        size_t operator()(CCMPMTKeyPair const & key) const {
            return CCMPMTKey::hash()(key.first) + (CCMPMTKey::hash()(key.second) << (8*3));
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

I3_CLASS_VERSION(CCMPMTKeyPair, ccmpmtkeypair_version_);

inline bool operator<(CCMPMTKeyPair const & lhs, CCMPMTKeyPair const & rhs) {
    if(lhs.first < rhs.first) {
        return true;
    } else if(lhs.first > rhs.first) {
        return false;
    } else if(lhs.second < rhs.second) {
        return true;
    } else {
        return false;
    }
}

inline bool operator>(CCMPMTKeyPair const & lhs, CCMPMTKeyPair const & rhs) {
    if(lhs.first > rhs.first) {
        return true;
    } else if(lhs.first < rhs.first) {
        return false;
    } else if(lhs.second > rhs.second) {
        return true;
    } else {
        return false;
    }
}

/**
 * streams an CCMPMTKeyPair to an arbitrary ostream. These are important,
 * the tray uses these conversion internally.
 */
std::ostream& operator<<(std::ostream&, CCMPMTKeyPair const & p);
std::istream& operator>>(std::istream&, CCMPMTKeyPair & p);
#pragma GCC diagnostic pop

I3_POINTER_TYPEDEFS(CCMPMTKeyPair);

#endif //CCMPMTKEYPAIR_H_INCLUDED
