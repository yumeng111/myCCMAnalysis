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

CalculateNLLH::CalculateNLLH(I3FramePtr data_frame, I3FramePtr geo_frame, size_t max_bins, size_t n_sodium_events, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height, std::vector<CCMPMTKey> keys_to_fit) :
    max_bins(max_bins), n_sodium_events(n_sodium_events), portion_light_reflected_by_tpb(portion_light_reflected_by_tpb), desired_chunk_width(desired_chunk_width), desired_chunk_height(desired_chunk_height), keys_to_fit(keys_to_fit)
{
    SetData(data_frame);
    SetGeo(geo_frame);
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
        SinglePMTInfo pmt_data;
        pmt_data.key = i->first;
        size_t peak_idx = std::distance(i->second.begin(), std::max_element(i->second.begin(), i->second.end()));
        size_t start_idx = std::max(size_t(3), peak_idx) - 3;
        size_t min_idx = std::max(size_t(15), peak_idx) - 15;
        size_t max_idx = std::min(min_idx + max_bins, i->second.size());
        pmt_data.data = std::vector(i->second.begin() + min_idx, i->second.begin() + max_idx);
        pmt_data.start_time = (start_idx - min_idx) * 2.0;
        pmt_data.max_time = (std::max(int(min_idx), int(max_idx) - 1) - min_idx) * 2.0;
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

template<typename T>
T CalculateNLLH::ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type) {

    // let's grab our data
    SinglePMTInfo const & pmt_data = data[key];
    double start_time = pmt_data.start_time;
    double max_time = pmt_data.max_time;

    // let's grab our expectation
    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> pred = gen_expectation->GetExpectation(key, start_time, max_time, Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB,
            normalization, light_time_offset, uv_absorption, z_offset, n_sodium_events, light_profile_type);

    // unpack into yields and yields^2
    boost::shared_ptr<std::vector<T>> pred_yields;
    boost::shared_ptr<std::vector<T>> pred_yields_squared;
    pred_yields = std::get<0>(pred);
    pred_yields_squared = std::get<1>(pred);

    T total_nllh(0.0);

    size_t n_bins = pmt_data.data.size();

    // now loop over the time bins in our pred
    for (size_t i=0; i<n_bins; ++i) {
        double k = pmt_data.data.at(i);
        T mu = pred_yields->at(i) * normalization;
        T sigma_squared = pred_yields_squared->at(i) * (normalization * normalization);
        total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
    }

    return total_nllh;
}

CalculateNLLH::AD CalculateNLLH::GetNLLH(CCMPMTKey key, AnalyticLightYieldGenerator const & params) {
    return this->ComputeNLLH<AD>(key, AD(params.Rs, 0), AD(params.Rt, 1), AD(params.tau_s, 2), AD(params.tau_t, 3), AD(params.tau_rec, 4), AD(params.tau_TPB, 5),
            AD(params.normalization, 6), AD(params.time_offset, 7), params.uv_absorption, params.z_offset, params.n_sodium_events, params.light_profile_type);
}

double CalculateNLLH::GetNLLHValue(CCMPMTKey key, AnalyticLightYieldGenerator const & params) {
    return this->ComputeNLLH<double>(key, params.Rs, params.Rt, params.tau_s, params.tau_t, params.tau_rec, params.tau_TPB,
            params.normalization, params.time_offset, params.uv_absorption, params.z_offset, params.n_sodium_events, params.light_profile_type);
}

CalculateNLLH::Grad CalculateNLLH::GetNLLHDerivative(CCMPMTKey key, AnalyticLightYieldGenerator const & params) {
    Grad grad;
    this->GetNLLH(key, params).copyGradient(grad.data());
    return grad;
}


