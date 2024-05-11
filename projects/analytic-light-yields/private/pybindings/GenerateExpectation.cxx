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

#include <analytic-light-yields/CalculateNLLH.h>
#include <analytic-light-yields/GenerateExpectation.h>
#include <icetray/python/copy_suite.hpp>
#include <icetray/python/indexed_property.hpp>
#include <icetray/python/boost_serializable_pickle_suite.hpp>
#include <icetray/python/dataclass_suite.hpp>
#include <icetray/python/stream_to_string.hpp>

using namespace boost::python;

typedef CalculateNLLH::AD AD;
typedef CalculateNLLH::Grad Grad;
constexpr int n_params = CalculateNLLH::n_params;

std::vector<double> LightProfileValue(GenerateExpectation & g, double Rs, double Rt, double tau_s, double tau_t, double tau_rec, double tau_TPB, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<double> const & times) {
    return g.LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_profile_type, times);
}

std::vector<std::tuple<double, Grad>> LightProfile(GenerateExpectation & g, double Rs, double Rt, double tau_s, double tau_t, double tau_rec, double tau_TPB, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<double> const & times) {
    std::vector<AD> grad_times(times.size());
    for (size_t i = 0; i < times.size(); i++) {
        grad_times[i] = AD(times[i], 7);
    }

    std::vector<AD> result = g.LightProfile(AD(Rs, 0), AD(Rt, 1), AD(tau_s, 3), AD(tau_t, 3), AD(tau_rec, 4), AD(tau_TPB, 5), light_profile_type, grad_times);
    std::vector<std::tuple<double, Grad>> ret_val;
    for (size_t i = 0; i < result.size(); i++) {
        double val = result[i].value();
        Grad grad;
        result[i].copyGradient(grad.data());
        ret_val.push_back(std::make_tuple(val, grad));
    }
    return ret_val;
}

std::vector<Grad> LightProfileGrad(GenerateExpectation & g, double Rs, double Rt, double tau_s, double tau_t, double tau_rec, double tau_TPB, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<double> const & times) {
    std::vector<AD> grad_times(times.size());
    for (size_t i = 0; i < times.size(); i++) {
        grad_times[i] = AD(times[i], 7);
    }

    std::vector<AD> result = g.LightProfile(AD(Rs, 0), AD(Rt, 1), AD(tau_s, 3), AD(tau_t, 3), AD(tau_rec, 4), AD(tau_TPB, 5), light_profile_type, grad_times);
    std::vector<Grad> ret_val;
    for (size_t i = 0; i < result.size(); i++) {
        Grad grad;
        result[i].copyGradient(grad.data());
        ret_val.push_back(grad);
    }
    return ret_val;
}

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
GetExpectationGrad(GenerateExpectation & g, CCMPMTKey key, double start_time, double max_time, double peak_time, double Rs, double Rt, double tau_s, double tau_t, double tau_rec, double tau_TPB,
            double normalization, double light_time_offset, double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type) {

    std::tuple<boost::shared_ptr<std::vector<AD>>, boost::shared_ptr<std::vector<AD>>> result = g.GetExpectation<AD>(key, start_time, max_time, peak_time, AD(Rs, 0), AD(Rt, 1), AD(tau_s, 3), AD(tau_t, 3), AD(tau_rec, 4), AD(tau_TPB, 5), AD(normalization, 6), AD(light_time_offset, 7), uv_absorption, z_offset, n_sodium_events, light_profile_type);

    std::vector<std::vector<double>> grads(std::get<0>(result)->size(), std::vector<double>(n_params));
    std::vector<std::vector<double>> grads_squared(std::get<1>(result)->size(), std::vector<double>(n_params));
    for (size_t i = 0; i < std::get<0>(result)->size(); i++) {
        std::get<0>(result)->at(i).copyGradient((grads)[i].data());
        std::get<1>(result)->at(i).copyGradient((grads_squared)[i].data());
    }
    return std::tie(grads, grads_squared);
}

void register_GenerateExpectation() {
{

    class_<GenerateExpectation, boost::noncopyable>("GenerateExpectation")
        .def(init<std::vector<CCMPMTKey>,size_t,I3FramePtr,double,double,double>())
        .def("GetSodiumVertices", &GenerateExpectation::GetSodiumVertices)
        .def("GetYieldsAndOffsets", &GenerateExpectation::GetYieldsAndOffsets)
        .def("LightProfile", LightProfile)
        .def("LightProfileValue", LightProfileValue)
        .def("LightProfileGrad", LightProfileGrad)
        //.def("DLightProfile", &GenerateExpectation::DLightProfile)
        //.def("GetExpectationGrad", GetExpectationGrad)
    ;
}
}


