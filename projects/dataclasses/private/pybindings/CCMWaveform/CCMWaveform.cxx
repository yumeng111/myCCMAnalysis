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

template<typename T>
std::string to_str(const CCMWaveform<T> theWaveform){
    std::ostringstream oss;
    oss << theWaveform << std::flush;
    return oss.str();
}

template<typename T>
void register_CCMWaveform_T(std::string name) {
  std::vector<T>& (CCMWaveform<T>::*get_waveform)() = &CCMWaveform<T>::GetWaveform;
  const std::vector<CCMStatusCompound>& 
    (CCMWaveform<T>::*get_waveform_information)() const = &CCMWaveform<T>::GetWaveformInformation;
  object get_waveform_func = make_function(get_waveform, return_internal_reference<1>());
  object get_waveform_information_func = make_function(get_waveform_information, return_internal_reference<1>());
  unsigned (*get_status_static)(const std::vector<CCMStatusCompound>&) = &CCMWaveform<T>::GetStatus;
  unsigned (CCMWaveform<T>::*get_status_member)() const = &CCMWaveform<T>::GetStatus;

  std::string class_name = "CCMWaveform" + name;

  {
    scope waveform_scope =
      class_<CCMWaveform<T>, boost::shared_ptr<CCMWaveform<T>> >(class_name.c_str())
      .def(copy_suite<CCMWaveform<T>>())
      .add_property(snake_case("Source"), &CCMWaveform<T>::GetSource)
      .add_property("status", get_status_member)
      .add_property("time", &CCMWaveform<T>::GetStartTime, &CCMWaveform<T>::SetStartTime)
      .add_property("binwidth", &CCMWaveform<T>::GetBinWidth, &CCMWaveform<T>::SetBinWidth)
      .add_property("waveform", get_waveform_func, &CCMWaveform<T>::SetWaveform)
      .add_property("waveform_information", get_waveform_information_func, &CCMWaveform<T>::SetWaveformInformation)

      // for static methods you need the both of these
      .def("get_status", get_status_static)
      .staticmethod("get_status")
      .def(self == self)
      .def(dataclass_suite<CCMWaveform<T>>())
      .def("__str__", to_str<T>)
     ;
  }
}

void register_CCMWaveform() {

    const std::pair<unsigned long long int, unsigned long long int>&
      (CCMStatusCompound::* get_interval)() const 
      = &CCMStatusCompound::GetInterval;

    class_<CCMStatusCompound>("CCMStatusCompound")
      #define PROPS (Status)
      BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, CCMStatusCompound, PROPS)
      #undef PROPS
      .add_property("interval", make_function(get_interval, return_value_policy<copy_const_reference>()))
      .def( freeze() )
      ;

    enum_<CCMSource>("CCMSource")
      .value("V1730", CCMSource::V1730)
      .export_values()
      ;

    enum_<CCMStatus>("CCMStatus")
      .value("VIRGINAL", CCMStatus::VIRGINAL)
      .value("COMBINED", CCMStatus::COMBINED)
      .value("SATURATED", CCMStatus::SATURATED)
      .value("UNDERSHOT", CCMStatus::UNDERSHOT)
      .export_values()
      ;

    register_CCMWaveform_T<uint16_t>("UInt16");
}
