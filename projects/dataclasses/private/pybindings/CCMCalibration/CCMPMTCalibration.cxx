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

#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/indexed_property.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>
#include <icetray/python/stream_to_string.hpp>

using namespace boost::python;

template<typename theType>
std::string to_str(const theType& theStr){
    std::ostringstream oss;
    oss << theStr << std::flush;
    return oss.str();
}

std::string to_str_CCMSPETemplate(const CCMSPETemplate& self){
    return to_str<CCMSPETemplate>(self); };
std::string to_str_CCMPMTCalibrationMap(const CCMPMTCalibrationMap& self){
    return to_str<CCMPMTCalibrationMap>(self); };

namespace {
    std::string to_str_SPEChargeDistribution(const SPEChargeDistribution& self){
        return to_str<SPEChargeDistribution>(self); };
};

void register_CCMPMTCalibration() {
{
    class_<CCMSinglePulseParameters, boost::shared_ptr<CCMSinglePulseParameters>>("CCMSinglePulseParameters")
        .def_readwrite("peak_height", &CCMSinglePulseParameters::peak_height)
        .def_readwrite("relative_peak_time", &CCMSinglePulseParameters::relative_peak_time)
        .def_readwrite("rise_time", &CCMSinglePulseParameters::rise_time)
        .def_readwrite("duration", &CCMSinglePulseParameters::duration)
        .def(dataclass_suite<CCMSinglePulseParameters>());

    class_<CCMSinglePulseParametersSeries, 
        CCMSinglePulseParametersSeriesPtr>("CCMSinglePulseParametersSeries")
            .def(dataclass_suite<CCMSinglePulseParametersSeries>())
            ;

    class_<CCMSPETemplate, boost::shared_ptr<CCMSPETemplate>>("CCMSPETemplate")
#define CCMSPETemplatePROPS (PulseParameters)
        BOOST_PP_SEQ_FOR_EACH(WRAP_PROP_INTERNAL_REFERENCE, CCMSPETemplate, CCMSPETemplatePROPS)
#undef CCMSPETemplatePROPS
        .def("EvaluateSinglePulse", &CCMSPETemplate::EvaluateSinglePulse)
        .staticmethod("EvaluateSinglePulse")
        .def("Evaluate", &CCMSPETemplate::Evaluate)
        .def(dataclass_suite<CCMSPETemplate>());

    scope outer = 
        class_<CCMPMTCalibration, boost::shared_ptr<CCMPMTCalibration> >("CCMPMTCalibration")
#define I3DOMCALPROPS (PulseStartTime)(PulseEndTime)(DroopTimeConstant)(PMTGain)(PMTDeltaT)(PMTRelEff)(MeanPMTCharge)(PMTMeanCharge)
        BOOST_PP_SEQ_FOR_EACH(WRAP_PROP, CCMPMTCalibration, I3DOMCALPROPS)
#undef I3DOMCALPROPS
#define I3DOMCALPROPS (SPEChargeDistribution)(SPETemplate)
        BOOST_PP_SEQ_FOR_EACH(WRAP_PROP_INTERNAL_REFERENCE, CCMPMTCalibration, I3DOMCALPROPS)
#undef I3DOMCALPROPS
        .def("__str__", to_str_CCMSPETemplate)
        .def("__str__", to_str_SPEChargeDistribution)
        .def("__str__", to_str_CCMPMTCalibrationMap)
        .def(dataclass_suite<CCMPMTCalibration>())
        ;

    class_<CCMPMTCalibrationMap, 
        CCMPMTCalibrationMapPtr>("CCMPMTCalibrationMap")
            .def(dataclass_suite<CCMPMTCalibrationMap>())
            ;
}
}
