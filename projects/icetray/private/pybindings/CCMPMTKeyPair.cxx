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
#include <icetray/CCMPMTKeyPair.h>
#include <icetray/python/boost_serializable_pickle_suite.hpp>

using namespace boost::python;

    inline static unsigned
hash_ccmpmtkeypair (CCMPMTKeyPair const & key) {
    return CCMPMTKeyPair::hash()(key);
}

typedef CCMPMTKeyPair value_type;
// make CCMPMTKeyPair iterable:
static object ccmpmtkeypair_getitem(value_type const & x, int i) {
    if (i==0 || i==-2) return object(x.first);
    else if (i==1 || i==-1) return object(x.second);
    else {
        PyErr_SetString(PyExc_IndexError,"Index out of range.");
        throw_error_already_set();
        return object(); // None
    }
}
// __len__ = 2
static int ccmpmtkeypair_len(value_type const& x) { return 2; }

std::string repr(CCMPMTKeyPair const & key) {
    std::stringstream s;
    s << key;
    return s.str();
}

    void
register_CCMPMTKeyPair() {
    class_<CCMPMTKeyPair>("CCMPMTKeyPair")
        .def(init<>())
        .def(init<CCMPMTKey, CCMPMTKey>())
        .def_readwrite("first", &CCMPMTKeyPair::first)
        .def_readwrite("second", &CCMPMTKeyPair::second)
        .def("__str__", &CCMPMTKeyPair::str)
        .def("__repr__", repr)
        .def("__hash__", hash_ccmpmtkeypair)
        .def("__getitem__", ccmpmtkeypair_getitem)
        .def("__len__", ccmpmtkeypair_len)
        .def(self == self)
        .def(self != self)
        .def(self < self)
        .def_pickle(boost_serializable_pickle_suite<CCMPMTKeyPair>())
        ;

    from_python_sequence<std::vector<CCMPMTKeyPair>, variable_capacity_policy>();
}
