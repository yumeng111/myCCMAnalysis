/**
 *  $Id$
 *
 *  Copyright (C) 2007
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
#include <icetray/serialization.h>
#include <icetray/CCMPMTKeyPair.h>
#include <sstream>

template <class Archive>
void CCMPMTKeyPair::save(Archive& ar, unsigned version) const {
    ar & make_nvp("first", first);
    ar & make_nvp("second", second);
}

template <class Archive>
void CCMPMTKeyPair::load(Archive& ar, unsigned version) {
    if (version>ccmpmtkeypair_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMPMTKeyPair class.",
                version,ccmpmtkeypair_version_);

    ar & make_nvp("first", first);
    ar & make_nvp("second", second);
}

I3_SPLIT_SERIALIZABLE(CCMPMTKeyPair);

#include <sstream>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

ostream& operator<<(ostream& os, CCMPMTKeyPair const & p) {
    os << "CCMPMTKeyPair(" << p.first << "," << p.second << ")";
    return os;
}

istream& operator>>(istream& is, CCMPMTKeyPair & p) {
    std::string s;
    is >> s;

    boost::regex reg_00("CCMPMTKeyPair\\(CCMPMTKey\\(([\\-\\d]+),(\\d+)\\),CCMPMTKey\\(([\\-\\d]+),(\\d+)\\)\\)");
    boost::regex reg_01("CCMPMTKeyPair\\(CCMPMTKey\\(([\\-\\d]+),(\\d+)\\),CCMPMTKey\\(([\\-\\d]+),(\\d+),(\\d+)\\)\\)");
    boost::regex reg_10("CCMPMTKeyPair\\(CCMPMTKey\\(([\\-\\d]+),(\\d+),(\\d+)\\),CCMPMTKey\\(([\\-\\d]+),(\\d+)\\)\\)");
    boost::regex reg_11("CCMPMTKeyPair\\(CCMPMTKey\\(([\\-\\d]+),(\\d+),(\\d+)\\),CCMPMTKey\\(([\\-\\d]+),(\\d+),(\\d+)\\)\\)");
    boost::smatch matches;
    CCMPMTKey newkey;
    if(boost::regex_search(s, matches, reg_11)) {
        // try the 3 entry version first
        log_trace("matches: %s %s %s, %s %s %s",
                matches.str(1).c_str(), matches.str(2).c_str(), matches.str(3).c_str(),
                matches.str(4).c_str(), matches.str(5).c_str(), matches.str(6).c_str()
        );

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(1)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(2)));
        newkey.SetSubsensor(static_cast<unsigned char>(boost::lexical_cast<unsigned int>(matches.str(3))));
        p.first = newkey;

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(4)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(5)));
        newkey.SetSubsensor(static_cast<unsigned char>(boost::lexical_cast<unsigned int>(matches.str(6))));
        p.second = newkey;
        return is;
    } else if(boost::regex_search(s, matches, reg_00)) {
        log_trace("matches: %s %s [old-style CCMPMTKey], %s %s [old-style CCMPMTKey]",
                matches.str(1).c_str(), matches.str(2).c_str(),
                matches.str(3).c_str(), matches.str(4).c_str()
        );

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(1)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(2)));
        newkey.SetSubsensor(0);
        p.first = newkey;

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(3)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(4)));
        newkey.SetSubsensor(0);
        p.second = newkey;
        return is;
    } else if(boost::regex_search(s, matches, reg_01)) {
        log_trace("matches: %s %s [old-style CCMPMTKey], %s %s %s",
                matches.str(1).c_str(), matches.str(2).c_str(),
                matches.str(3).c_str(), matches.str(4).c_str(), matches.str(5).c_str()
        );

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(1)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(2)));
        newkey.SetSubsensor(0);
        p.first = newkey;

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(3)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(4)));
        newkey.SetSubsensor(static_cast<unsigned char>(boost::lexical_cast<unsigned int>(matches.str(5))));
        p.second = newkey;
        return is;
    } else if(boost::regex_search(s, matches, reg_10)) {
        log_trace("matches: %s %s %s, %s %s [old-style CCMPMTKey]",
                matches.str(1).c_str(), matches.str(2).c_str(), matches.str(3).c_str(),
                matches.str(4).c_str(), matches.str(5).c_str()
        );

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(1)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(2)));
        newkey.SetSubsensor(static_cast<unsigned char>(boost::lexical_cast<unsigned int>(matches.str(3))));
        p.first = newkey;

        newkey.SetRegion(boost::lexical_cast<int>(matches.str(4)));
        newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(5)));
        p.second = newkey;
        return is;
    } else {
        log_fatal("Error parsing CCMPMTKeyPair value from string \"%s\"", s.c_str());
    }
}

std::string CCMPMTKeyPair::str() const {
    stringstream s;
    s<<*this;
    return s.str();
}



