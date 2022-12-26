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
#include <icetray/CCMTriggerKey.h>
#include <sstream>
#include <icetray/I3FrameObject.h>

CCMTriggerKey::~CCMTriggerKey() { }

template <class Archive>
void CCMTriggerKey::save(Archive& ar, unsigned version) const
{
  ar & make_nvp("TriggerType",  triggerType_);
  ar & make_nvp("Number",  triggerNumber_);
}

template <class Archive>
void CCMTriggerKey::load(Archive& ar, unsigned version)
{
  if (version>ccmtriggerkey_version_)
    log_fatal("Attempting to read version %u from file but running version %u of CCMTriggerKey class.",
              version,ccmtriggerkey_version_);

  I3FrameObject fo;
  ar & make_nvp("TriggerType",  triggerType_);
  ar & make_nvp("Number",  triggerNumber_);
}

I3_SPLIT_SERIALIZABLE(CCMTriggerKey);

#include <sstream>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

ostream& operator<<(ostream& os, const CCMTriggerKey& key)
{
  os << "CCMTriggerKey(" << static_cast<int>(key.GetType()) << "," << key.GetNumber() << ")";
  return os;
}

istream& operator>>(istream& is, CCMTriggerKey& key)
{
  std::string s;
  is >> s;

  boost::regex reg("CCMTriggerKey\\(([\\-\\d]+),(\\d+)\\)");
  boost::smatch matches;
  if (boost::regex_search(s, matches, reg)) {
    CCMTriggerKey newkey;
    log_trace("matches: %s %s", matches.str(1).c_str(), matches.str(2).c_str());
    newkey.SetType(static_cast<CCMTriggerKey::TriggerType>(boost::lexical_cast<int>(matches.str(1))));
    newkey.SetNumber(boost::lexical_cast<unsigned int>(matches.str(2)));
    key = newkey;
    return is;
  } else {
      log_fatal("Error parsing CCMTriggerKey value from string \"%s\"", s.c_str());
  }
}

string CCMTriggerKey::str()const{
  stringstream s;
  s<<*this;
  return s.str();
}



