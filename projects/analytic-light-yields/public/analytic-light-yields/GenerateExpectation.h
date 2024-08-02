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

#include "icetray/ctpl.h"
#include "icetray/open.h"
#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3TrayInfo.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/CCMPMTKey.h"
#include "dataio/I3File.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/physics/HESodiumEvent.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"
#include <analytic-light-yields/SodiumVertexDistribution.h>
#include <analytic-light-yields/YieldsPerPMT.h>
#include <analytic-light-yields/G4YieldsPerPMT.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>
#include "analytic-light-yields/autodiff.h"

template<typename T>
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
    std::shared_ptr<G4YieldsPerPMT> g4_yields_and_offset_constructor = nullptr;
    boost::shared_ptr<HESodiumEventSeries> event_vertices = boost::make_shared<HESodiumEventSeries> ();
    double portion_light_reflected_by_tpb;
    double desired_chunk_width;
    double desired_chunk_height;
public:
    GenerateExpectation();
    GenerateExpectation(I3VectorCCMPMTKey keys_to_fit, size_t n_sodium_events, I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height);
    void GetSodiumVertices(size_t n_events_to_simulate, T z_position);
    //void GetYieldsAndOffsets(CCMPMTKey key, T uv_absorption);
    //void ComputeBinnedYield(CCMPMTKey key, double max_time);
    void RunMultiThreadedCode(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, T photons_per_mev, double max_time);
    void GrabG4Yields(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, T photons_per_mev, double max_time, T z_loc);

    std::vector<T> LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                  AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times);

    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> GetExpectation(CCMPMTKey key,
                       double start_time, double max_time, double peak_time,
                       T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                       std::vector<T> uv_absorption, T photons_per_mev, T z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type,
                       bool UseG4Yields); 
    
    std::vector<double> light_prof_debug;
    std::vector<double> light_prof_times_debug;
    std::vector<double> light_prof_time_offset_grad_debug;

    std::vector<double> GetLightProfileDebug() {return light_prof_debug;};
    std::vector<double> GetLightProfileTimesDebug() {return light_prof_times_debug;};
    std::vector<double> GetLightProfileTOffsetGradDebug() {return light_prof_time_offset_grad_debug;};
    
    // now define some class members    
    std::map<CCMPMTKey, std::vector<T>> binned_yields;
    std::map<CCMPMTKey, std::vector<T>> binned_square_yields;
    std::vector<boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>>> yields_per_pmt_per_event;

    std::vector<T> uv_absorption_calculated = {0.0, 0.0};
    T z_offset_calculated = 0.0;
    
    std::vector<std::deque<I3FramePtr>> G4Events;
    std::vector<double> all_z_offsets;
    bool grabbed_g4_events = false;

    void GrabG4Events(size_t n_events, T z_offset);

};

template<typename T> void GenerateExpectation<T>::GrabG4Events(size_t n_events, T z_offset){
    
    std::vector<std::string> g4_fnames;

    if (z_offset == 30.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium32cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium30cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium28cmHEEvents.i3.zst"};
        all_z_offsets = {32.0, 30.0, 28.0};
        //g4_fnames = {"/Users/darcybrewuser/workspaces/CCM/notebooks/G4Sodium30cmHEEvents.i3.zst"};
        //all_z_offsets = {30.0};
    }
    if (z_offset == -30.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium-28cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium-30cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium-32cmHEEvents.i3.zst"};
        all_z_offsets = {-28.0, -30.0, -32.0};
        //g4_fnames = {"/Users/darcybrewuser/workspaces/CCM/notebooks/G4Sodium-30cmHEEvents.i3.zst"};
        //all_z_offsets = {-30.0};
    }
    if (z_offset == 0.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium2cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium0cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium-2cmHEEvents.i3.zst"};
        all_z_offsets = {2.0, 0.0, -2.0};
    }
    if (z_offset == 50.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium52cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium50cmHEEvents.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/G4Sodium48cmHEEvents.i3.zst"};
        all_z_offsets = {52.0, 50.0, 48.0};
    }

    for (size_t f = 0; f < g4_fnames.size(); f++){
        std::deque<I3FramePtr> this_file_events;

        dataio::I3File g4_file(g4_fnames.at(f), dataio::I3File::Mode::read);
        size_t frames_grabbed = 0;
        while (frames_grabbed < n_events or !g4_file.more()){
            I3FramePtr g4_frame = g4_file.pop_frame();
            this_file_events.push_back(g4_frame);
            frames_grabbed += 1;
        }
        std::cout << "grabbed " << this_file_events.size() << " events for sodium at " << all_z_offsets.at(f) << "cm" << std::endl;
        G4Events.push_back(this_file_events);
    }

    grabbed_g4_events = true;
}

template<typename T> std::vector<T> GenerateExpectation<T>::LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                                    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times) {
    std::vector<T> light_profile(times.size());

    if (light_profile_type == AnalyticLightYieldGenerator::LArLightProfileType::Simplified) {
        get_light_profile_no_recombination(Rs, tau_s, tau_t, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, times, light_profile);
    } else {
        get_total_light_profile (Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, times, light_profile);
    }

    return light_profile;
}

template<typename T> std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>>
    GenerateExpectation<T>::GetExpectation(CCMPMTKey key, double start_time, double max_time, double peak_time, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale, std::vector<T> uv_absorption, T photons_per_mev, T z_offset, size_t n_sodium_events,
            AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, bool UseG4Yields) {

    if (grabbed_g4_events == false and UseG4Yields){
        GrabG4Events(n_sodium_events, z_offset);
    }

    // now see if we need to get our yields    
    bool compute_vertices = sodium_events_constructor == nullptr;
    bool compute_yields;
    if (UseG4Yields){
        compute_yields = g4_yields_and_offset_constructor == nullptr;
    } else {
        compute_yields = yields_and_offset_constructor == nullptr;
    }

    // one last check that we've already got the yields for key at this uv absorption and z
    if (compute_yields == false){
        // check if uv absorption has changed!!!
        if (uv_absorption.at(0) != uv_absorption_calculated.at(0) and uv_absorption.at(1) != uv_absorption_calculated.at(1)){
            binned_yields.clear();
            binned_square_yields.clear();
            compute_yields = true;
        }

        // check if z has changed
        if (z_offset != z_offset_calculated){
            binned_yields.clear();
            binned_square_yields.clear();
            compute_yields = true;
        }
    }

    // check that we made our sodium event vertices
    if(compute_vertices and UseG4Yields == false)
        GetSodiumVertices(n_sodium_events, z_offset);
    // now let's check if we got our yields + time offsets
    if(compute_yields){
        //GetYieldsAndOffsets(key, uv_absorption);
        //ComputeBinnedYield(key, max_time);
        uv_absorption_calculated = uv_absorption;
        z_offset_calculated = z_offset;
        if (UseG4Yields){
            GrabG4Yields(keys_to_fit, uv_absorption, photons_per_mev, 2.0 * max_time, z_offset);
        }
        else {
            RunMultiThreadedCode(keys_to_fit, uv_absorption, photons_per_mev, 2.0 * max_time);
        }
    }

    std::vector<T> const & binned_yield = binned_yields.at(key);
    std::vector<T> const & binned_square_yield = binned_square_yields.at(key);
   
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
            expectation->at(i) += binned_yield.at(k) * LAr_light_profile.at(j);
            expectation_squared->at(i) += binned_square_yield.at(k) * LAr_light_profile_squared.at(j); 
        }
    }
    return std::make_tuple(expectation, expectation_squared, extended_light_times);
}

template<typename T> GenerateExpectation<T>::GenerateExpectation() :
    keys_to_fit(I3VectorCCMPMTKey()), geo_frame(geo_frame), n_sodium_events(n_sodium_events), desired_chunk_width(desired_chunk_width), desired_chunk_height(desired_chunk_height) {}

template<typename T> GenerateExpectation<T>::GenerateExpectation(I3VectorCCMPMTKey keys_to_fit, size_t n_sodium_events, I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height) :
    keys_to_fit(keys_to_fit), geo_frame(geo_frame), n_sodium_events(n_sodium_events), desired_chunk_width(desired_chunk_width), desired_chunk_height(desired_chunk_height) {}

template<typename T> void GenerateExpectation<T>::GetSodiumVertices(size_t n_events_to_simulate, T z_position) {
    sodium_events_constructor = std::make_shared<SodiumVertexDistribution> ();
    if constexpr (std::is_same<T, double>::value) {
        event_vertices = sodium_events_constructor->GetEventVertices(n_events_to_simulate, z_position);
    } else {
        event_vertices = sodium_events_constructor->GetEventVertices(n_events_to_simulate, z_position.value());
    }
}

template<typename T> void GenerateExpectation<T>::RunMultiThreadedCode(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, T photons_per_mev, double max_time){
    yields_and_offset_constructor = std::make_shared<YieldsPerPMT>(geo_frame, portion_light_reflected_by_tpb, desired_chunk_width, desired_chunk_height);

    size_t n_threads = 0;
    yields_and_offset_constructor->GetAllYields(n_threads, event_vertices, uv_absorption.at(0), keys_to_fit, max_time, photons_per_mev,
                                                binned_yields, binned_square_yields);
}

template<typename T> void GenerateExpectation<T>::GrabG4Yields(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, T photons_per_mev, double max_time, T z_loc){
    g4_yields_and_offset_constructor = std::make_shared<G4YieldsPerPMT>();

    size_t n_threads = 0;
    
    // let's grab set of g4 events at the z offset above and below z_loc
    double z_loc_double;
    if constexpr (std::is_same<T, double>::value) {
        z_loc_double = z_loc;
    } else {
        z_loc_double = z_loc.value();
    }

    size_t upper_idx;
    size_t lower_idx;
    bool need_to_interpolate = true;

    for (size_t k = 0; k < all_z_offsets.size(); k++){
       // first check if z_loc_double is in all_z_offsets
       if (all_z_offsets.at(k) == z_loc_double){
            need_to_interpolate = false;
            upper_idx = k;
            lower_idx = k;
        }
            
        // now check if this z offset is greater than z_loc_double
        if (all_z_offsets.at(k) > z_loc_double and need_to_interpolate){
             upper_idx = k;
             lower_idx = upper_idx + 1;
        }
        
    }
    if (need_to_interpolate == false){
        std::cout << "for z = " <<  z_loc_double << " exactly in our list at " << all_z_offsets[upper_idx] << std::endl;
    } else {
        std::cout << "for z = " << z_loc_double << ", closest z above = " << all_z_offsets[upper_idx] << " and closest z below = " << all_z_offsets[lower_idx] << std::endl;
    }
    
    g4_yields_and_offset_constructor->GetAllYields(n_threads, keys_to_fit, G4Events.at(upper_idx), G4Events.at(lower_idx), need_to_interpolate, max_time,
                                                   all_z_offsets.at(upper_idx), all_z_offsets.at(lower_idx), uv_absorption, photons_per_mev, z_loc, binned_yields, binned_square_yields);
}

//template<typename T> void GenerateExpectation<T>::GetYieldsAndOffsets(CCMPMTKey key, T uv_absorption) {
//    yields_per_pmt_per_event.clear();
//    binned_yields.clear();
//    binned_square_yields.clear();
//    yields_and_offset_constructor = std::make_shared<YieldsPerPMT>(geo_frame, portion_light_reflected_by_tpb, desired_chunk_width, desired_chunk_height);
//    
//    // now loop over events and get map between CCMPMTKey and std::vector<photon_yield_summary> 
//    for (size_t sodium_it = 0; sodium_it < event_vertices->size(); ++sodium_it) {
//        boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> yields_per_event = yields_and_offset_constructor->GetAllYields(event_vertices->at(sodium_it), uv_absorption, {key});
//        yields_per_pmt_per_event.push_back(yields_per_event);
//        for (typename std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>::const_iterator i = yields_per_event->begin(); i != yields_per_event->end(); i++) {
//            std::vector<photon_yield_summary<T>> const & yields = i->second;
//            if(yields.size() == 0) {
//                continue;
//            }
//        }
//    }
//}
//
//template<typename T> void GenerateExpectation<T>::ComputeBinnedYield(CCMPMTKey key, double max_time) {
//    size_t n_bins = max_time / 2.0;
//    binned_yields[key] = std::vector<T>(n_bins, 0.0);
//    binned_square_yields[key] = std::vector<T>(n_bins, 0.0);
//    std::vector<T> & binned_yields_per_pmt = binned_yields[key];
//    std::vector<T> & binned_square_yields_per_pmt = binned_square_yields[key];
//
//    for (size_t sodium_it = 0; sodium_it < event_vertices->size(); ++sodium_it) {
//        boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> yields_per_event = yields_per_pmt_per_event.at(sodium_it);
//        typename std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>::const_iterator i = yields_per_event->find(key);
//        if (i == yields_per_event->end()) {
//            continue;
//        }
//        std::vector<photon_yield_summary<T>> const & yields = i->second;
//        if(yields.size() == 0) {
//            continue;
//        }
//        for(size_t yield_it = 0; yield_it < yields.size(); ++yield_it) {
//            photon_yield_summary<T> const & yield = yields.at(yield_it);
//            size_t bin_idx;
//            if constexpr (std::is_same<T, double>::value) {
//                bin_idx = yield.time / 2.0;
//            } else {
//                bin_idx = yield.time.value() / 2.0;
//            }
//            if(bin_idx >= n_bins) {
//                continue;
//            }
//            binned_yields_per_pmt.at(bin_idx) += yield.yield;
//            binned_square_yields_per_pmt.at(bin_idx) += yield.yield * yield.yield;
//        }
//    }
//    binned_yields[key] = binned_yields_per_pmt;
//    binned_square_yields[key] = binned_square_yields_per_pmt;
//}


#endif
