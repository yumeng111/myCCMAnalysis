/**
 *  $Id$
 *  
 *  Copyright (C) 2007
 *  Troy D. Straszheim  <troy@icecube.umd.edu>
 *  and the IceCube Collaboration <http://www.icecube.wisc.edu>
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
#include <icetray/CCMChannelKey.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>

using namespace boost::python;

inline static unsigned 
hash_ccmchannelkey (const CCMChannelKey& key)
{
  return CCMChannelKey::hash()(key); 
}

typedef CCMChannelKey value_type;
// make CCMChannelKey iterable: string,om,pmt = CCMChannelKey
static object ccmchannelkey_getitem(value_type const& x, int i) {
    if (i==0 || i==-2) return object(x.GetBoard());
    else if (i==1 || i==-1) return object(x.GetChannel()); 
    else {
        PyErr_SetString(PyExc_IndexError,"Index out of range.");
        throw_error_already_set();
        return object(); // None
    }
}
// __len__ = 3
static int ccmchannelkey_len(value_type const& x) { return 3; }

std::string repr(const CCMChannelKey& key){
  std::stringstream s;
  s << "CCMChannelKey(" << key.GetBoard() << "," << key.GetChannel() << ")";
  return s.str();
}

void
register_CCMChannelKey()
{
  class_<CCMChannelKey>("CCMChannelKey")
    .def(init<uint16_t, uint8_t>())
    PROPERTY(CCMChannelKey, board, Board)
    PROPERTY(CCMChannelKey, channel, Channel)
    .def("__str__", &CCMChannelKey::str)
    .def("__repr__", repr)
    .def("__hash__", hash_ccmchannelkey)
    .def("__getitem__", ccmchannelkey_getitem)
    .def("__len__", ccmchannelkey_len)
    .def(self == self)
    .def(self != self)
    .def(self < self)
    .def_pickle(boost_serializable_pickle_suite<CCMChannelKey>())
    ;

  from_python_sequence<std::vector<CCMChannelKey>, variable_capacity_policy>();
}
