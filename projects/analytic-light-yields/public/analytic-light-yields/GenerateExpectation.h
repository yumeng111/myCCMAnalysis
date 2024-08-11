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
#include "simclasses/CCMMCPE.h"
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
    void GrabG4Yields(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, T photons_per_mev, double max_time, T rayl, T z_offset);

    std::vector<T> LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                  AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times);

    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> GetExpectation(CCMPMTKey key,
                       double start_time, double max_time, double peak_time,
                       T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                       std::vector<T> uv_absorption, T rayl, T photons_per_mev, T z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type,
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
    T rayl_calculated = 0.0;
    T z_calculated = 0.0;

    bool grabbed_g4_events = false;
    std::vector<std::pair<float, float>> z_rayl_indexing; // idx of z, rayl used for simulation
    std::vector<std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>>> collated_events; // first idx corresponds to z, rayl in map
                                                                                                     // second idx corresponds to each event
                                                                                                     // map is between CCMPMTKey and photon hits for each event
    void GrabG4Events(size_t n_events, T z_offset);

};

template<typename T> void GenerateExpectation<T>::GrabG4Events(size_t n_events, T z_offset){

    std::vector<std::string> g4_fnames;

    if (z_offset >= 26.0 and z_offset <= 34.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset95.0cmRayleigh2628Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset85.0cmRayleigh2578Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset75.0cmRayleigh2463Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset65.0cmRayleigh2592Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset55.0cmRayleigh2529Events.i3.zst"};
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset45.0cmRayleigh2395Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset35.0cmRayleigh2335Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset25.0cmRayleigh2138Events.i3.zst"};
    }
    if (z_offset >= -34.0 and z_offset <= -26.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset95.0cmRayleigh2732Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset85.0cmRayleigh2713Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset75.0cmRayleigh2816Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset65.0cmRayleigh2622Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset55.0cmRayleigh2596Events.i3.zst"};
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset45.0cmRayleigh2509Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset35.0cmRayleigh2438Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset25.0cmRayleigh2355Events.i3.zst"};
    }
    if (z_offset >= -4.0 and z_offset <= 4.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset95.0cmRayleigh3488Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset85.0cmRayleigh3383Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset75.0cmRayleigh3448Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset65.0cmRayleigh3355Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset55.0cmRayleigh3252Events.i3.zst"};
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset45.0cmRayleigh3288Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset35.0cmRayleigh3114Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset25.0cmRayleigh2906Events.i3.zst"};
    }
    if (z_offset >= 46.0 and z_offset <= 54.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset95.0cmRayleigh1303Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset85.0cmRayleigh1353Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset75.0cmRayleigh1370Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset65.0cmRayleigh1316Events.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset55.0cmRayleigh1359Events.i3.zst"};
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset45.0cmRayleigh1281Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset35.0cmRayleigh1259Events.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset25.0cmRayleigh1292Events.i3.zst"};
    }

    for (size_t f = 0; f < g4_fnames.size(); f++){
        // start with some string parsing to get z and rayl
        std::string path = g4_fnames.at(f);
        size_t n_sodium = 6;
        size_t z_offset_pre = path.find("Sodium");
        size_t z_offset_post = path.find("cmZOffset");
        float this_z_offset = std::stoi(path.substr(z_offset_pre + n_sodium , z_offset_post - (z_offset_pre + n_sodium)));

        size_t n_zoffset = 7;
        size_t rayl_pre = path.find("ZOffset");
        size_t rayl_post = path.find("cmRayleigh");
        float this_rayl = std::stoi(path.substr(rayl_pre + n_zoffset , rayl_post - (rayl_pre + n_zoffset)));

        // now set up object to save each event
        std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> this_simulation_events;

        // now loop over frames
        dataio::I3File g4_file(g4_fnames.at(f), dataio::I3File::Mode::read);
        size_t frames_grabbed = 0;
        while (frames_grabbed < n_events or !g4_file.more()){
            // set up object to hold data
            std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>> this_event_info;

            // grab our frame
            I3FramePtr g4_frame = g4_file.pop_frame();

            // grab + loop through ccmmcpe object
            boost::shared_ptr<CCMMCPESeriesMap const> mcpeseries_source = g4_frame->Get<boost::shared_ptr<CCMMCPESeriesMap const>>("PMTMCHitsMap");

            // Iterate over PMTs in source map
            for (auto it = mcpeseries_source->begin(); it != mcpeseries_source->end(); ++it) {

                // check if this pmt is in this_event_info
                if(this_event_info.find(it->first) == this_event_info.end()) {
                    this_event_info[it->first] = std::vector<ccmmcpe_lightweight>();
                }

                // Iterate over the vector of CCMMCPE in the source map for this PMT
                for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                    CCMMCPE const & pe = *it2;

                    // check the wavelength
                    if (pe.wavelength / I3Units::nanometer < 325.0){
                        continue;
                    }

                    // make a new ccmmcpe_lightweight
                    ccmmcpe_lightweight pe_light;
                    pe_light.g4_time = pe.g4_time;
                    pe_light.g4_distance_uv = pe.g4_distance_uv;

                    // now save
                    this_event_info[it->first].push_back(pe_light);
                }
            }

            // now save this event info!
            this_simulation_events.push_back(this_event_info);
            frames_grabbed += 1;
        }
        // now save this_simulation_events
        collated_events.push_back(this_simulation_events);
        z_rayl_indexing.push_back(std::make_pair(this_z_offset, this_rayl));
        std::cout << "grabbed " << this_simulation_events.size() << " events for sodium at z = " << this_z_offset << "cm and rayl = " << this_rayl << "cm" << std::endl;
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
            T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale, std::vector<T> uv_absorption, T rayl, T photons_per_mev, T z_offset, size_t n_sodium_events,
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

    // one last check that we've already got the yields for key at this uv absorption and rayleigh scattering length
    if (compute_yields == false){
        // check if uv absorption has changed!!!
        if (uv_absorption.at(0) != uv_absorption_calculated.at(0) and uv_absorption.at(1) != uv_absorption_calculated.at(1)){
            binned_yields.clear();
            binned_square_yields.clear();
            compute_yields = true;
        }

        // check if rayl has changed
        if (rayl != rayl_calculated){
            binned_yields.clear();
            binned_square_yields.clear();
            compute_yields = true;
        }

        // check if z has changed
        if (z_offset != z_calculated){
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
        rayl_calculated = rayl;
        z_calculated = z_offset;
        if (UseG4Yields){
            GrabG4Yields(keys_to_fit, uv_absorption, photons_per_mev, 2.0 * max_time, rayl, z_offset);
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

template<typename T> void GenerateExpectation<T>::GrabG4Yields(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, T photons_per_mev, double max_time, T rayl, T z_offset){
    g4_yields_and_offset_constructor = std::make_shared<G4YieldsPerPMT>();

    size_t n_threads = 0;
    
    double z_offset_double;
    double rayl_double;
    if constexpr (std::is_same<T, double>::value) {
        z_offset_double = z_offset;
        rayl_double = rayl;
    } else {
        z_offset_double = z_offset.value();
        rayl_double = rayl.value();
    }

    // let's figure out the z and rayl above / below where we want to evaluate
    size_t Q11; // idx of z below and rayl below
    size_t Q12; // idx of z below and rayl above
    size_t Q21; // idx of z above and rayl below
    size_t Q22; // idx of z above and rayl above

    for (size_t m = 0; m < z_rayl_indexing.size(); m++){
        float this_z = z_rayl_indexing.at(m).first;
        float this_rayl = z_rayl_indexing.at(m).second;

        float delta_z = std::abs(this_z - z_offset_double);
        float delta_rayl = std::abs(this_rayl - rayl_double);

        if (delta_z <= 2.0 and delta_rayl <= 10.0){
            // ok so we are around the z, rayl pair we want
            if (this_z <= z_offset_double and this_rayl <= rayl_double){
                Q11 = m;
            }
            if (this_z <= z_offset_double and this_rayl > rayl_double){
                Q12 = m;
            }
            if (this_z > z_offset_double and this_rayl <= rayl_double){
                Q21 = m;
            }
            if (this_z > z_offset_double and this_rayl > rayl_double){
                Q22 = m;
            }
        }
    }

    std::cout << "Q11{z, rayl} = " << z_rayl_indexing.at(Q11) << ", Q12{z, rayl} = " << z_rayl_indexing.at(Q12)
        << ", Q21{z, rayl} = " << z_rayl_indexing.at(Q21) << ", Q22{z, rayl} = " << z_rayl_indexing.at(Q22) << std::endl;

    g4_yields_and_offset_constructor->GetAllYields(n_threads, keys_to_fit, max_time,
                                                   collated_events.at(Q11), collated_events.at(Q12), collated_events.at(Q21),collated_events.at(Q22),
                                                   z_rayl_indexing.at(Q11), z_rayl_indexing.at(Q12), z_rayl_indexing.at(Q21), z_rayl_indexing.at(Q22),
                                                   uv_absorption, z_offset, rayl, binned_yields, binned_square_yields);
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
