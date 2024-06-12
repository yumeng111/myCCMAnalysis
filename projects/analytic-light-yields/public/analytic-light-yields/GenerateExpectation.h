#ifndef GenerateExpectation_H
#define GenerateExpectation_H

#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

#include <map>
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

#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/physics/PhotonYieldSummary.h"
#include "dataclasses/physics/HESodiumEvent.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"
#include <analytic-light-yields/SodiumVertexDistribution.h>
#include <analytic-light-yields/YieldsPerPMT.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>
#include "analytic-light-yields/autodiff.h"

class GenerateExpectation {
    I3FramePtr geo_frame;
    size_t n_sodium_events = 0;
    I3VectorCCMPMTKey keys_to_fit;

    // things we need to do in order to generate expected hits/pmt :
    // 1) get soidum vertex distribution
    // 2) get yields + time offsets per pmt
    // 3) get light profile
    // then we put it all together and bin!
    // so for now, we are going to set sodium vertices ONCE and get yields + offsets ONCE
    // and update light profile with every call to GetExpectation
    std::shared_ptr<SodiumVertexDistribution> sodium_events_constructor = nullptr;
    std::shared_ptr<YieldsPerPMT> yields_and_offset_constructor = nullptr;

    boost::shared_ptr<HESodiumEventSeries> event_vertices = boost::make_shared<HESodiumEventSeries> ();
    std::vector<boost::shared_ptr<PhotonYieldSummarySeriesMap>> yields_per_pmt_per_event;
    // std::map<CCMPMTKey, std::tuple<double, double>> time_extents_per_pmt;
    std::map<CCMPMTKey, std::vector<double>> binned_yields;
    std::map<CCMPMTKey, std::vector<double>> binned_square_yields;

    double uv_absorption = 0.0;
    double z_offset = 0.0;

    double portion_light_reflected_by_tpb;
    double desired_chunk_width;
    double desired_chunk_height;

    void ComputeBinnedYield(CCMPMTKey key, double max_time);

public:
    GenerateExpectation();
    GenerateExpectation(I3VectorCCMPMTKey keys_to_fit, size_t n_sodium_events, I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height);
    void GetSodiumVertices(size_t n_events_to_simulate, double z_position);
    void GetYieldsAndOffsets(double uv_absorption);
    template<typename T>
    std::vector<T> LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                       AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times);

    template<typename T>
    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> GetExpectation(CCMPMTKey key, double start_time, double max_time, double peak_time,
                                            T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                            double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type); 
    
    //static constexpr size_t n_params = 9;
    //typedef double Underlying;

    //typedef phys_tools::autodiff::FD<n_params, Underlying> AD;
    //typedef std::array<Underlying, n_params> Grad;
    
    std::vector<double> light_prof_debug;
    std::vector<double> light_prof_times_debug;
    std::vector<double> light_prof_time_offset_grad_debug;

    std::vector<double> GetLightProfileDebug() {return light_prof_debug;};
    std::vector<double> GetLightProfileTimesDebug() {return light_prof_times_debug;};
    std::vector<double> GetLightProfileTOffsetGradDebug() {return light_prof_time_offset_grad_debug;};

};

template<typename T>
std::vector<T> GenerateExpectation::LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                                    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times) {
    std::vector<T> light_profile(times.size());

    if (light_profile_type == AnalyticLightYieldGenerator::LArLightProfileType::Simplified) {
        get_light_profile_no_recombination(Rs, tau_s, tau_t, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, times, light_profile);
    } else {
        get_total_light_profile (Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, times, light_profile);
    }

    return light_profile;
}

template<typename T>
std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>>
    GenerateExpectation::GetExpectation(CCMPMTKey key, double start_time, double max_time, double peak_time, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale, double uv_absorption, double z_offset, size_t n_sodium_events,
            AnalyticLightYieldGenerator::LArLightProfileType light_profile_type) {

    bool compute_vertices = sodium_events_constructor == nullptr;
    bool compute_yields = yields_and_offset_constructor == nullptr;
    bool bin_yields = binned_yields.find(key) == binned_yields.end();

    // check that we made our sodium event vertices
    if(compute_vertices)
        GetSodiumVertices(n_sodium_events, z_offset);
    // now let's check if we get our yields + time offsets
    if(compute_yields)
        GetYieldsAndOffsets(uv_absorption);
    // now lets check if we have pre-binned the yields
    if(bin_yields)
        ComputeBinnedYield(key, max_time);

    std::vector<double> const & binned_yield = binned_yields.at(key);
    std::vector<double> const & binned_square_yield = binned_square_yields.at(key);

    peak_time += 40.0;
    size_t n_light_bins = size_t(max_time / 2.0) + 1 + 40;
    size_t offset = size_t(start_time / 2.0);
    size_t n_data_bins = n_light_bins - offset;

    std::vector<T> light_times(n_light_bins, 0.0);
    for(size_t i = 0; i < n_light_bins; i++) {
        light_times[i] = 2.0 * i - peak_time;
    }
    
    // now grab our light profile!
    std::vector<T> LAr_light_profile = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, light_profile_type, light_times);
    std::vector<T> LAr_light_profile_squared(n_light_bins);
    for(size_t i = 0; i < n_light_bins; i++) {
        LAr_light_profile_squared[i] = LAr_light_profile[i] * LAr_light_profile[i];
    }

    //// quick debug -- grabbing our light profile again but w very fine time binning
    // empty all our vectors
    //light_prof_debug.clear();
    //light_prof_times_debug.clear();
    //light_prof_time_offset_grad_debug.clear();

    //double fine_grained_bin_spacing = 0.001;
    //size_t n_light_bins_fine_grained = size_t(max_time / fine_grained_bin_spacing) + 1;
    //std::vector<T> light_times_debug(n_light_bins_fine_grained, 0.0);
    //for(size_t i = 0; i < n_light_bins_fine_grained; i++) {
    //    light_times_debug[i] = fine_grained_bin_spacing * i + (light_time_offset - peak_time);
    //    if constexpr (std::is_same<T, double>::value) {
    //        light_prof_times_debug.push_back(light_times_debug[i]);
    //    }
    //}
    //
    //std::vector<T> LAr_light_profile_debug = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_profile_type, light_times_debug);
    //for(size_t i = 0; i < n_light_bins_fine_grained; i++) {
    //    // saving light profile value
    //    if constexpr (std::is_same<T, double>::value) {
    //        light_prof_debug.push_back(LAr_light_profile_debug[i]);
    //    }
    //    // saving time offset gradient
    //    if constexpr (std::is_same<T, AD>::value) {
    //        Grad grad;
    //        LAr_light_profile_debug[i].copyGradient(grad.data());
    //        light_prof_time_offset_grad_debug.push_back(grad.at(7));
    //    }
    //}

    //// done w debug!

    // make a vector of final binned yields to return
    boost::shared_ptr<std::vector<T>> expectation = boost::make_shared<std::vector<T>> (LAr_light_profile.size(), 0.0);
    boost::shared_ptr<std::vector<T>> expectation_squared = boost::make_shared<std::vector<T>> (LAr_light_profile_squared.size(), 0.0);
    boost::shared_ptr<std::vector<T>> extended_light_times = boost::make_shared<std::vector<T>> (LAr_light_profile.size(), 0.0);
    for(size_t i = 0; i < LAr_light_profile.size(); i++) {
        extended_light_times->at(i) = 2.0 * i - peak_time;
    }

    // Convolute the light profile with the binned yields
    for(size_t i = 0; i < n_data_bins; i++) {
        size_t tot_idx = i + offset;
        for(size_t j = 0; j <= tot_idx; j++) {
            size_t k = tot_idx - j;
            expectation->at(i) += binned_yield[k] * LAr_light_profile[j];
            expectation_squared->at(i) += binned_square_yield[k] * LAr_light_profile_squared[j]; 
        }
    }
    
    return std::make_tuple(expectation, expectation_squared, extended_light_times);
}


#endif
