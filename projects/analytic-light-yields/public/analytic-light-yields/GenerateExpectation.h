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

#include <tuple>
#include <vector>
#include <random>
#include <cmath>
#include <map>
#include <memory>

#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/physics/PhotonYieldSummary.h"
#include "dataclasses/physics/HESodiumEvent.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"
#include <analytic-light-yields/SodiumVertexDistribution.h>
#include <analytic-light-yields/YieldsPerPMT.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>

class GenerateExpectation {
    I3FramePtr geo_frame;
    size_t n_sodium_events = 0;
    std::vector<CCMPMTKey> keys_to_fit;

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

    // let's also add small noise to every bin
    double noise_photons = 1.0;
    double noise_triggers = 5.0;
    double digitization_time = 16.0 * std::pow(10.0, 3.0); //16 usec in nsec
    double noise_rate = noise_photons / (noise_triggers * digitization_time); // units of photons/nsec
    double noise_rate_per_time_bin = 2.0 * noise_rate; // 2nsec binning

    double uv_absorption = 0.0;
    double z_offset = 0.0;

    double portion_light_reflected_by_tpb;
    double desired_chunk_width;
    double desired_chunk_height;

    void ComputeBinnedYield(CCMPMTKey key, double max_time);

public:
    GenerateExpectation();
    GenerateExpectation(std::vector<CCMPMTKey> keys_to_fit, size_t n_sodium_events, I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height);
    void GetSodiumVertices(size_t n_events_to_simulate, double z_position);
    void GetYieldsAndOffsets(double uv_absorption);
    template<typename T>
    std::vector<T> LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                       AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times);

    template<typename T>
    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> GetExpectation(CCMPMTKey key, double start_time, double max_time, double peak_time, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type);
};

template<typename T>
std::vector<T> GenerateExpectation::LightProfile(T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                                    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, std::vector<T> const & times) {
    std::vector<T> light_profile(times.size());

    if (light_profile_type == AnalyticLightYieldGenerator::LArLightProfileType::Simplified) {
        get_light_profile_no_recombination(Rs, Rt, tau_s, tau_t, tau_TPB, times, light_profile);
    } else {
        get_total_light_profile (Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, times, light_profile);
    }

    return light_profile;
}

template<typename T>
std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>>
    GenerateExpectation::GetExpectation(CCMPMTKey key, double start_time, double max_time, double peak_time, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type) {

    bool compute_vertices = sodium_events_constructor == nullptr or n_sodium_events != this->n_sodium_events or z_offset != this->z_offset;
    bool compute_yields = yields_and_offset_constructor == nullptr or uv_absorption != this->uv_absorption;
    bool bin_yields = binned_yields.find(key) == binned_yields.end();
    compute_yields |= compute_vertices;
    bin_yields |= compute_yields;

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

    size_t n_light_bins = size_t(max_time / 2.0) + 1;
    size_t offset = size_t(start_time / 2.0);
    size_t n_data_bins = n_light_bins - offset;

    std::vector<T> light_times(n_light_bins, 0.0);
    for(size_t i = 0; i < n_light_bins; i++) {
        light_times[i] = 2.0 * i + (light_time_offset - peak_time);
    }

    // now grab our light profile!
    std::vector<T> LAr_light_profile = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_profile_type, light_times);
    std::vector<T> LAr_light_profile_squared(n_light_bins);
    for(size_t i = 0; i < n_light_bins; i++) {
        LAr_light_profile_squared[i] = LAr_light_profile[i] * LAr_light_profile[i];
    }

    // make a vector of final binned yields to return
    boost::shared_ptr<std::vector<T>> final_binned_yields_ = boost::make_shared<std::vector<T>> ();
    boost::shared_ptr<std::vector<T>> final_binned_squared_yields_ = boost::make_shared<std::vector<T>> ();
    std::vector<T> & expectation = *final_binned_yields_;
    std::vector<T> & expectation_squared = *final_binned_squared_yields_;

    // Convolute the light profile with the binned yields
    for(size_t i = 0; i < n_data_bins; i++) {
        size_t tot_idx = i + offset;
        for(size_t j = 0; j <= tot_idx; j++) {
            size_t k = tot_idx - j;
            expectation[i] += binned_yield[k] * LAr_light_profile[j];
            expectation_squared[i] += binned_square_yield[k] * LAr_light_profile_squared[j];
        }
    }

    return std::make_tuple(final_binned_yields_, final_binned_squared_yields_);
}


#endif
