//
//   Copyright (c) 2004, 2005, 2006, 2007   Troy D. Straszheim  
//   
//   $Id$
//
//   This file is part of IceTray.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions
//   are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//   
//   THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
//   ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//   ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
//   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
//   OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
//   OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
//   SUCH DAMAGE.
//   
//   SPDX-License-Identifier: BSD-2-Clause
//   
//

// this pragma has to go before any functions or whatever are defined.
#pragma GCC diagnostic ignored "-Wwrite-strings"

#include <vector>

#include <dataclasses/CCMTime.h>
#include <Python.h>
#include <datetime.h>
#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

#define HAVE_PYDATETIME_API

std::string repr(CCMTime t){
    std::stringstream out;
    out <<  "CCMTime(";
    t.Print(out);
    out << ")";
    return out.str();
}

void set_unix_time_default(CCMTime& t, time_t unixTime, double ns) {
    t.SetUnixTime(unixTime, ns);
}

void register_CCMTime() {
    scope i3time_scope = class_<CCMTime, bases<I3FrameObject>, 
          boost::shared_ptr<CCMTime> >("CCMTime")
              .def(init<int64_t, int64_t, double>())
              .def(init<int64_t, int64_t>())
              .def(init<const CCMTime&>())
              .def("__str__",&stream_to_string<CCMTime>)
              .def("__repr__",&repr)
              .def("get_utc_string",&CCMTime::GetUTCString)
              .def(self != self)
              .def(self == self)
              .def(self-self)
              .def(self-double())
              .def(self+double())
              .def(self<self)
              .def(self>self)
              .def(self<=self)
              .def(self>=self)
              .def(dataclass_suite<CCMTime>())
              ;

    register_pointer_conversions<CCMTime>();
}
