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
#include <icetray/CCMPMTKey.h>
#include <sstream>
#include <icetray/I3FrameObject.h>

CCMPMTKey::~CCMPMTKey() { }

template <class Archive>
void CCMPMTKey::save(Archive& ar, unsigned version) const
{
  ar & make_nvp("RegionNumber",  regionNumber_);
  ar & make_nvp("SensorNumber",  sensorNumber_);
  ar & make_nvp("SubsensorNumber",  subsensorNumber_);
}

template <class Archive>
void CCMPMTKey::load(Archive& ar, unsigned version)
{
  if (version>ccmpmtkey_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMPMTKey class.",
              version,ccmpmtkey_version_);

  I3FrameObject fo;
  ar & make_nvp("RegionNumber",  regionNumber_);
  ar & make_nvp("SensorNumber",  sensorNumber_);
  ar & make_nvp("SubsensorNumber",  subsensorNumber_);
}

I3_SPLIT_SERIALIZABLE(CCMPMTKey);

#include <sstream>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

ostream& operator<<(ostream& os, const CCMPMTKey& key)
{
  os << "CCMPMTKey(" << key.GetRegion() << "," << key.GetSensor() << "," << static_cast<unsigned int>(key.GetSubsensor()) << ")";
  return os;
}

istream& operator>>(istream& is, CCMPMTKey& key)
{
  std::string s;
  is >> s;

  boost::regex reg("CCMPMTKey\\(([\\-\\d]+),(\\d+),(\\d+)\\)");
  boost::smatch matches;
  if (boost::regex_search(s, matches, reg)) {
    // try the 3 entry version first
    CCMPMTKey newkey;
    log_trace("matches: %s %s %s", matches.str(1).c_str(), matches.str(2).c_str(), matches.str(3).c_str());
    newkey.SetRegion(boost::lexical_cast<int>(matches.str(1)));
    newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(2)));
    newkey.SetSubsensor(static_cast<unsigned char>(boost::lexical_cast<unsigned int>(matches.str(3))));
    key = newkey;
    return is;
  } else {
    boost::regex reg("CCMPMTKey\\(([\\-\\d]+),(\\d+)\\)");
    boost::smatch matches;
    if (!boost::regex_search(s, matches, reg))
      log_fatal("Error parsing CCMPMTKey value from string \"%s\"", s.c_str());

    CCMPMTKey newkey;
    log_trace("matches: %s %s [old-style CCMPMTKey]", matches.str(1).c_str(), matches.str(2).c_str());
       
    newkey.SetRegion(boost::lexical_cast<int>(matches.str(1)));
    newkey.SetSensor(boost::lexical_cast<unsigned>(matches.str(2)));
    newkey.SetSubsensor(0);
    key = newkey;
    return is;
  }
}

string CCMPMTKey::str()const{
  stringstream s;
  s<<*this;
  return s.str();
}



