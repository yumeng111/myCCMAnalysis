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
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"
#include <analytic-light-yields/G4YieldsPerPMT.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>
#include "analytic-light-yields/autodiff.h"


template<typename T>
class GenerateExpectation {
    I3VectorCCMPMTKey keys_to_fit;
    std::shared_ptr<G4YieldsPerPMT> g4_yields_and_offset_constructor = nullptr;
public:
    GenerateExpectation();
    GenerateExpectation(I3VectorCCMPMTKey keys_to_fit);

    void GrabG4Events(size_t n_events, T z_offset);
    void GrabG4Yields(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption, double max_time, T rayl, T z_offset,
                      double light_times_per_bin, bool fit_z_rayl);
    std::vector<T> LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                AnalyticLightYieldGenerator::LArLightProfileType light_profile_type,
                                std::vector<T> const & times);

    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> GetExpectation(CCMPMTKey key,
                                                            double start_time, double max_time, double peak_time,
                                                            T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                                            T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                                            std::vector<T> uv_absorption, T rayl, T z_offset,
                                                            size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type,
                                                            bool fit_z_rayl);

    // now define some class members
    std::map<CCMPMTKey, std::vector<T>> binned_yields;
    std::map<CCMPMTKey, std::vector<T>> binned_square_yields;

    // set up some default values
    std::vector<T> uv_absorption_calculated = {0.0, 0.0};
    T rayl_calculated = 0.0;
    T z_calculated = 0.0;

    // set up storage for geant4 simulation data
    bool grabbed_g4_events = false;
    std::vector<std::pair<float, float>> z_rayl_indexing; // idx of z, rayl used for simulation
    std::vector<std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>>> collated_events; // first idx corresponds to z, rayl in map
                                                                                                     // second idx corresponds to each event
                                                                                                     // map is between CCMPMTKey and photon hits for each event
};

template<typename T> void GenerateExpectation<T>::GrabG4Events(size_t n_events, T z_offset){

    std::vector<std::string> g4_fnames;

    if (z_offset >= 26.0 and z_offset <= 34.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium26.0cmZOffset85.0cmRayleigh6952Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium26.0cmZOffset95.0cmRayleigh7168Events_v2.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium26.0cmZOffset75.0cmRayleigh6818Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset85.0cmRayleigh6350Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset95.0cmRayleigh6324Events_v2.i3.zst"};
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium30.0cmZOffset75.0cmRayleigh6417Events_v2.i3.zst"};
    }
    if (z_offset >= -34.0 and z_offset <= -26.0){
        g4_fnames = {//"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-34.0cmZOffset75.0cmRayleigh5451Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-34.0cmZOffset85.0cmRayleigh5517Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-34.0cmZOffset95.0cmRayleigh5640Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset85.0cmRayleigh6583Events_v2.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset75.0cmRayleigh6451Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-30.0cmZOffset95.0cmRayleigh6690Events_v2.i3.zst"};
    }
    if (z_offset >= -4.0 and z_offset <= 4.0){
        g4_fnames = {"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-4.0cmZOffset85.0cmRayleigh8537Events_v2.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-4.0cmZOffset75.0cmRayleigh8350Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium-4.0cmZOffset95.0cmRayleigh8672Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset85.0cmRayleigh8559Events_v2.i3.zst",
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset75.0cmRayleigh8446Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium0.0cmZOffset95.0cmRayleigh8582Events_v2.i3.zst"};
    }
    if (z_offset >= 46.0 and z_offset <= 54.0){
        g4_fnames = {//"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium46.0cmZOffset75.0cmRayleigh6027Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium46.0cmZOffset85.0cmRayleigh6051Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium46.0cmZOffset95.0cmRayleigh6087Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset95.0cmRayleigh4043Events_v2.i3.zst",
                     "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset85.0cmRayleigh4076Events_v2.i3.zst"};
                     //"/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/UpdatedGeometryG4Sodium50.0cmZOffset75.0cmRayleigh4075Events_v2.i3.zst"};
    }

    for (size_t f = 0; f < g4_fnames.size(); f++){
        std::cout << "grabbing events from " << g4_fnames.at(f) << std::endl;
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

                // check if we are fitting to this key
                if (std::find(keys_to_fit.begin(), keys_to_fit.end(), it->first) == keys_to_fit.end()) {
                    continue;
                }

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

template<typename T> std::vector<T> GenerateExpectation<T>::LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                                                         T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                                                         AnalyticLightYieldGenerator::LArLightProfileType light_profile_type,
                                                                         std::vector<T> const & times) {
    std::vector<T> light_profile(times.size());

    if (light_profile_type == AnalyticLightYieldGenerator::LArLightProfileType::Simplified) {
        get_light_profile_no_recombination(Rs, tau_s, tau_t, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, times, light_profile);
    } else {
        get_total_light_profile (Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, times, light_profile);
    }

    return light_profile;
}

template<typename T> std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> GenerateExpectation<T>::GetExpectation(CCMPMTKey key,
                                                                                double start_time, double max_time, double peak_time,
                                                                                T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                                                                T light_time_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale,
                                                                                std::vector<T> uv_absorption, T rayl, T z_offset, size_t n_sodium_events,
                                                                                AnalyticLightYieldGenerator::LArLightProfileType light_profile_type,
                                                                                bool fit_z_rayl) {

    // let's grab out light profile first
    size_t light_times_per_bin = 3;
    size_t n_data_bins = (max_time + 1) * 2.0;
    size_t n_light_bins = n_data_bins * light_times_per_bin;

    std::vector<T> light_times(n_light_bins, 0.0);
    double light_time_delta_t = 2.0 / light_times_per_bin;
    for(size_t i = 0; i < n_light_bins; i++) {
        // light_time_offset is relative to start_time
        light_times[i] = light_time_delta_t * i - (light_time_offset + start_time);
    }

    // now grab our light profile!
    std::vector<T> LAr_light_profile = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, late_pulse_mu, late_pulse_sigma, late_pulse_scale, light_profile_type, light_times);
    std::vector<T> LAr_light_profile_squared(n_light_bins);
    for(size_t i = 0; i < n_light_bins; i++) {
        LAr_light_profile_squared[i] = LAr_light_profile[i] * LAr_light_profile[i];
    }

    if (grabbed_g4_events == false){
        GrabG4Events(n_sodium_events, z_offset);
    }

    // now see if we need to get our yields
    bool compute_yields = g4_yields_and_offset_constructor == nullptr;

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

    // now let's check if we got our yields + time offsets
    if(compute_yields){
        uv_absorption_calculated = uv_absorption;
        rayl_calculated = rayl;
        z_calculated = z_offset;
        GrabG4Yields(keys_to_fit, uv_absorption, max_time, rayl, z_offset, light_times_per_bin, fit_z_rayl);
    }
    std::vector<T> const & binned_yield = binned_yields.at(key);
    std::vector<T> const & binned_square_yield = binned_square_yields.at(key);

    // make a vector of final binned yields to return
    boost::shared_ptr<std::vector<T>> expectation = boost::make_shared<std::vector<T>> (LAr_light_profile.size(), 0.0);
    boost::shared_ptr<std::vector<T>> expectation_squared = boost::make_shared<std::vector<T>> (LAr_light_profile_squared.size(), 0.0);

    // Convolute the light profile with the binned yields
    for(size_t j=0; j<n_light_bins; ++j) {
        double t_light = j * light_time_delta_t;
        T const & light_profile = LAr_light_profile.at(j);
        T const & light_profile_squared = LAr_light_profile_squared.at(j);
        for(size_t k=0; k<binned_yield.size(); ++k) {
            T const & yield = binned_yield.at(k);
            T const & square_yield = binned_square_yield.at(k);
            double const & yield_time = k * light_time_delta_t;
            int data_idx = (t_light + yield_time) / 2.0;
            if(data_idx < 0)
                continue;
            else if(data_idx >= n_data_bins)
                break;
            expectation->at(data_idx) += yield * light_profile;
            expectation_squared->at(data_idx) += square_yield * light_profile_squared;
        }
    }
    return std::make_tuple(expectation, expectation_squared);
}

template<typename T> GenerateExpectation<T>::GenerateExpectation() :
    keys_to_fit(I3VectorCCMPMTKey()) {}

template<typename T> GenerateExpectation<T>::GenerateExpectation(I3VectorCCMPMTKey keys_to_fit) : keys_to_fit(keys_to_fit) {}


template<typename T> void GenerateExpectation<T>::GrabG4Yields(I3VectorCCMPMTKey keys_to_fit, std::vector<T> uv_absorption,
                                                               double max_time, T rayl, T z_offset, double light_times_per_bin,
                                                               bool fit_z_rayl){
    g4_yields_and_offset_constructor = std::make_shared<G4YieldsPerPMT>();

    size_t n_threads = 5;

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
    size_t Q11 = 0; // idx of z below and rayl below
    size_t Q12 = 0; // idx of z below and rayl above
    size_t Q21 = 0; // idx of z above and rayl below
    size_t Q22 = 0; // idx of z above and rayl above

    for (size_t m = 0; m < z_rayl_indexing.size(); m++){
        float this_z = z_rayl_indexing.at(m).first;
        float this_rayl = z_rayl_indexing.at(m).second;
        float delta_z = std::abs(this_z - z_offset_double);
        float delta_rayl = std::abs(this_rayl - rayl_double);
        if (delta_z <= 4.0 and delta_rayl <= 10.0){
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


    float z_below = z_rayl_indexing.at(Q11).first;
    float rayl_below = z_rayl_indexing.at(Q11).second;
    float z_above = z_rayl_indexing.at(Q22).first;
    float rayl_above = z_rayl_indexing.at(Q22).second;

    //std::cout << "For desired z = " << z_offset_double << " and desired rayl = " << rayl_double << ", z below = " << z_below
    //    << ", rayl below = " << rayl_below
    //    << ", z above = " << z_above
    //    << ", and rayl above = " << rayl_above << std::endl;

    g4_yields_and_offset_constructor->GetAllYields(n_threads, keys_to_fit, max_time, light_times_per_bin,
                                                   collated_events.at(Q11), collated_events.at(Q12), collated_events.at(Q21),collated_events.at(Q22),
                                                   z_below, rayl_below, z_above, rayl_above,
                                                   uv_absorption, z_offset, rayl, fit_z_rayl, binned_yields, binned_square_yields);
}



#endif
