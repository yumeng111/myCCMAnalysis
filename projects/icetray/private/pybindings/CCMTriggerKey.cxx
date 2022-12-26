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
#include <icetray/CCMTriggerKey.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>

using namespace boost::python;

#define ENUM_DEF(r,data,T) .value(BOOST_PP_STRINGIZE(T), data::T)

inline static unsigned 
hash_ccmtriggerkey (const CCMTriggerKey& key)
{
  return CCMTriggerKey::hash()(key); 
}

typedef CCMTriggerKey value_type;
// make CCMTriggerKey iterable: string,om,pmt = CCMTriggerKey
static object ccmtriggerkey_getitem(value_type const& x, int i) {
    if (i==0 || i==-2) return object(static_cast<int>(x.GetType()));
    else if (i==1 || i==-1) return object(x.GetNumber()); 
    else {
        PyErr_SetString(PyExc_IndexError,"Index out of range.");
        throw_error_already_set();
        return object(); // None
    }
}
// __len__ = 3
static int ccmtriggerkey_len(value_type const& x) { return 2; }

std::string repr(const CCMTriggerKey& key){
  std::stringstream s;
  s << "CCMTriggerKey(" << static_cast<int>(key.GetType()) << "," << key.GetNumber() << ")";
  return s.str();
}

void
register_CCMTriggerKey()
{

  enum_<CCMTriggerKey::TriggerType>("TriggerType")
    BOOST_PP_SEQ_FOR_EACH(ENUM_DEF,CCMTriggerKey,CCMTriggerKey_H_CCMTriggerKey_TriggerType)
    .export_values()
    ;
  ;

  def("identity", identity_<CCMTriggerKey::TriggerType>);
  class_<CCMTriggerKey>("CCMTriggerKey")
    .def(init<CCMTriggerKey::TriggerType, unsigned int>())
    PROPERTY(CCMTriggerKey, type, Type)
    PROPERTY(CCMTriggerKey, number, Number)
    .def("__str__", &CCMTriggerKey::str)
    .def("__repr__", repr)
    .def("__hash__", hash_ccmtriggerkey)
    .def("__getitem__", ccmtriggerkey_getitem)
    .def("__len__", ccmtriggerkey_len)
    .def(self == self)
    .def(self != self)
    .def(self < self)
    .def_pickle(boost_serializable_pickle_suite<CCMTriggerKey>())
    ;

  from_python_sequence<std::vector<CCMTriggerKey>, variable_capacity_policy>();
}
