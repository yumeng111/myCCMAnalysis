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

#include <vector>

#include <dataclasses/physics/CCMWaveform.h>

#include <icetray/python/dataclass_suite.hpp>
#include <dataclasses/ostream_overloads.hpp>

using namespace boost::python;

std::string to_str(const CCMWaveform theWaveform){
    std::ostringstream oss;
    oss << theWaveform << std::flush;
    return oss.str();
}

void register_CCMWaveform()
{
  std::vector<double>& (CCMWaveform::*get_waveform)() = &CCMWaveform::GetWaveform;
  const std::vector<CCMWaveform::StatusCompound>& 
    (CCMWaveform::*get_waveform_information)() const = &CCMWaveform::GetWaveformInformation;
  object get_waveform_func = make_function(get_waveform, return_internal_reference<1>());
  object get_waveform_information_func = make_function(get_waveform_information, return_internal_reference<1>());
  unsigned (*get_status_static)(const std::vector<CCMWaveform::StatusCompound>&) = &CCMWaveform::GetStatus;
  unsigned (CCMWaveform::*get_status_member)() const = &CCMWaveform::GetStatus;

  {
    scope waveform_scope =
      class_<CCMWaveform, boost::shared_ptr<CCMWaveform> >("CCMWaveform")
      .def(copy_suite<CCMWaveform>())
      .add_property(snake_case("Source"), &CCMWaveform::GetSource)
      .add_property("status", get_status_member)
      .add_property("time", &CCMWaveform::GetStartTime, &CCMWaveform::SetStartTime)
      .add_property("binwidth", &CCMWaveform::GetBinWidth, &CCMWaveform::SetBinWidth)
      .add_property("waveform", get_waveform_func, &CCMWaveform::SetWaveform)
      .add_property("waveform_information", get_waveform_information_func, &CCMWaveform::SetWaveformInformation)

      // for static methods you need the both of these
      .def("get_status", get_status_static)
      .staticmethod("get_status")
      .def(self == self)
      .def(dataclass_suite<CCMWaveform>())
      .def("__str__", to_str)
     ;

    const std::pair<unsigned long long int, unsigned long long int>&
      (CCMWaveform::StatusCompound::* get_interval)() const 
      = &CCMWaveform::StatusCompound::GetInterval;

    class_<CCMWaveform::StatusCompound>("StatusCompound")
      #define PROPS (Status)
      BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, CCMWaveform::StatusCompound, PROPS)
      #undef PROPS
      .add_property("interval", make_function(get_interval, return_value_policy<copy_const_reference>()))
      .def( freeze() )
      ;

    enum_<CCMWaveform::Source>("Source")
      .value("V1730", CCMWaveform::V1730)
      .export_values()
      ;

    enum_<CCMWaveform::Status>("Status")
      .value("VIRGINAL", CCMWaveform::VIRGINAL)
      .value("COMBINED", CCMWaveform::COMBINED)
      .value("SATURATED", CCMWaveform::SATURATED)
      .value("UNDERSHOT", CCMWaveform::UNDERSHOT)
      .export_values()
      ;
  }

}
