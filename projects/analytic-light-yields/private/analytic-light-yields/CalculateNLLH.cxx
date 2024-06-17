#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <random>
#include <chrono>
#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <math.h>

#include <icetray/ctpl.h>
#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/I3Double.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include <dataclasses/geometry/CCMGeometry.h>
#include <analytic-light-yields/CalculateNLLH.h>

CalculateNLLH::CalculateNLLH() :
    max_bins(100), n_sodium_events(20), portion_light_reflected_by_tpb(1.0), desired_chunk_width(20.0), desired_chunk_height(20.0), keys_to_fit(I3VectorCCMPMTKey()) {}

CalculateNLLH::CalculateNLLH(I3FramePtr data_frame, I3FramePtr geo_frame, size_t max_bins, size_t n_sodium_events, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height, I3VectorCCMPMTKey keys_to_fit) :
    max_bins(max_bins), n_sodium_events(n_sodium_events), portion_light_reflected_by_tpb(portion_light_reflected_by_tpb), desired_chunk_width(desired_chunk_width), desired_chunk_height(desired_chunk_height), keys_to_fit(keys_to_fit)
{
    SetData(data_frame);
    SetGeo(geo_frame);
}

void CalculateNLLH::SetKeys(I3VectorCCMPMTKey keys){
    keys_to_fit = keys;

    // now let's set up pmt efficiencies that are all 1.0
    for (size_t pmt_it = 0; pmt_it < keys.size(); pmt_it++){
        PMT_eff[keys.at(pmt_it)] = 1.0; 
    }
}

void CalculateNLLH::SetPMTEff(I3MapPMTKeyDouble pmt_eff){
    // clear our defauly pmt efficiencies just to make sure 
    PMT_eff.clear();

    // now use pmt_eff to set efficiencies
    PMT_eff = pmt_eff; 
}

void CalculateNLLH::SetGeo(I3FramePtr geo_frame) {
    this->geo_frame = geo_frame;
    gen_expectation = boost::make_shared<GenerateExpectation>(keys_to_fit, n_sodium_events, geo_frame, portion_light_reflected_by_tpb, desired_chunk_width, desired_chunk_height);
}

void CalculateNLLH::SetData(I3FramePtr data_frame) {
    bool restrict_keys = keys_to_fit.size() > 0;
    // grab our data out of the frame
    I3MapPMTKeyVectorDoubleConstPtr data_map = data_frame->Get<I3MapPMTKeyVectorDoubleConstPtr>("AccumulatedEventsMap");
    for(I3MapPMTKeyVectorDouble::const_iterator i = data_map->begin(); i != data_map->end(); i++) {
        if(restrict_keys and std::find(keys_to_fit.begin(), keys_to_fit.end(), i->first) == keys_to_fit.end()) {
            continue;
        }
        // based on the way I accumulated data, the event starts around the 50th bin
        // so to be safe, let's use the first 40 bins for out pre-event average
        double pre_event_values = 0.0;
        double pre_event_deriv = 0.0;
        double pre_event_bins = 0.0;
        double time_counter = 0.0;
        double prev_data = 0.0;
        for(double const & data_points: i->second) {
            if (time_counter < 80.0){
                pre_event_values += data_points;
                pre_event_bins += 1.0;
                pre_event_deriv += abs(data_points - prev_data);
            } 
            prev_data = data_points;
            time_counter += 2.0;
        }

        SinglePMTInfo pmt_data;
        pmt_data.key = i->first;
        size_t peak_idx = std::distance(i->second.begin(), std::max_element(i->second.begin(), i->second.end()));
        pmt_data.peak_value = i->second.at(peak_idx); 
        size_t start_idx = std::max(size_t(3), peak_idx) - 3;
        size_t min_idx = std::max(size_t(15), peak_idx) - 15;
        size_t max_idx = std::min(min_idx + max_bins, i->second.size());
        pmt_data.data = std::vector(i->second.begin() + min_idx, i->second.begin() + max_idx);
        pmt_data.start_time = (start_idx - min_idx) * 2.0;
        pmt_data.max_time = (std::max(int(min_idx), int(max_idx) - 1) - min_idx) * 2.0 + 30.0;
        pmt_data.peak_time = (peak_idx - min_idx) * 2.0;

        // let's also get the event start bin using derivs
        bool found_start = false;
        double deriv_threshold = (pre_event_deriv / pre_event_bins) * 5.0;
        pmt_data.data_times.push_back(0.0);
        for (size_t data_it = 1; data_it < pmt_data.data.size(); data_it++){
            double deriv = pmt_data.data.at(data_it) - pmt_data.data.at(data_it - 1);
            if (deriv > deriv_threshold and found_start == false){
                pmt_data.event_start_bin = data_it - 1;
                found_start = true;
            }
            pmt_data.data_times.push_back(pmt_data.data_times[data_it - 1] + 2.0);
        }

        // let's subtract off our event start time from pmt_data.data_times
        for (size_t i = 0; i < pmt_data.data_times.size(); i++){
            pmt_data.data_times.at(i) -= (pmt_data.event_start_bin * 2.0);
        }
        pmt_data.pre_event_average = pre_event_values / pre_event_bins;
        data[pmt_data.key] = pmt_data;
    }
    n_data_events = data_frame->Get<I3Double>("TotalEventsPastCuts").value;

}

I3MapPMTKeyVectorDoublePtr CalculateNLLH::GetData() const {
    I3MapPMTKeyVectorDoublePtr data_map = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    for (std::map<CCMPMTKey, SinglePMTInfo>::const_iterator i = data.begin(); i != data.end(); i++) {
        data_map->insert(std::make_pair(i->first, i->second.data));
    }
    return data_map;
}

void CalculateNLLH::SetGenExpectation(boost::shared_ptr<GenerateExpectation> gen_expectation) {
    this->gen_expectation = gen_expectation;
}

CalculateNLLH::AD CalculateNLLH::GetNLLH(CCMPMTKey key, AnalyticLightYieldGenerator const & params, 
                                         const double & late_pulse_mu, const double & late_pulse_sigma, const double & late_pulse_scale, double pmt_efficiency) {
    return this->ComputeNLLH<AD>(key, AD(params.Rs, 0), AD(params.Rt, 1), AD(params.tau_s, 2), AD(params.tau_t, 3), AD(params.tau_rec, 4), AD(params.tau_TPB, 5),
            AD(params.normalization, 6), AD(params.time_offset, 7), AD(params.const_offset, 8), AD(late_pulse_mu, 9), AD(late_pulse_sigma, 10), AD(late_pulse_scale, 11), 
            AD(pmt_efficiency, 12), params.uv_absorption, params.z_offset, params.n_sodium_events, params.light_profile_type);
}

double CalculateNLLH::GetNLLHValue(CCMPMTKey key, AnalyticLightYieldGenerator const & params,
                                   const double & late_pulse_mu, const double & late_pulse_sigma, const double & late_pulse_scale, double pmt_efficiency) {
    return this->ComputeNLLH<double>(key, params.Rs, params.Rt, params.tau_s, params.tau_t, params.tau_rec, params.tau_TPB,
            params.normalization, params.time_offset, params.const_offset, late_pulse_mu, late_pulse_sigma, late_pulse_scale,
            pmt_efficiency, params.uv_absorption, params.z_offset, params.n_sodium_events, params.light_profile_type);
}

std::vector<double> CalculateNLLH::GetNLLHDerivative(CCMPMTKey key, AnalyticLightYieldGenerator const & params,
                                                     const double & late_pulse_mu, const double & late_pulse_sigma, const double & late_pulse_scale, double pmt_efficiency) {
    Grad grad;
    this->GetNLLH(key, params, late_pulse_mu, late_pulse_sigma, late_pulse_scale, pmt_efficiency).copyGradient(grad.data());
    // putting grad into a vector...don't want to deal with pybindings for arrays
    std::vector<double> grad_to_return;
    for (size_t i = 0; i < grad.size(); i++){
        grad_to_return.push_back(grad.at(i));
    }
    return grad_to_return;
}

