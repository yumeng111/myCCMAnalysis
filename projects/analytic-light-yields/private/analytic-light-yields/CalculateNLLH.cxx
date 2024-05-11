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

CalculateNLLH::CalculateNLLH(I3FramePtr data_frame, I3FramePtr geo_frame, std::vector<CCMPMTKey> keys_to_fit) :
    keys_to_fit(keys_to_fit)
{
    SetData(data_frame);
    SetGeo(geo_frame);
    gen_expectation = std::make_shared<GenerateExpectation> ();
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
        size_t max_idx = std::distance(i->second.begin(), std::max_element(i->second.begin(), i->second.end()));
        max_idx = std::min(max_idx, max_bins);
        size_t start_idx = std::max(size_t(3), max_idx) - 3;
        size_t min_idx = std::max(size_t(15), max_idx) - 15;
        pmt_data.data = std::vector(i->second.begin() + min_idx, i->second.begin() + max_idx);
        pmt_data.start_time = (start_idx - min_idx) * 2.0;
        pmt_data.peak_time = (peak_idx - min_idx) * 2.0;
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

void CalculateNLLH::SetGenExpectation(std::shared_ptr<GenerateExpectation> gen_expectation) {
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

    // for every pmt, we want to align data and pred on the peak bins
    // let's also grab the smallest pred peak bin to make sure we are not going too far to the left of the peaks
    size_t closest_pred_peak_bin = 3;
    for (I3MapPMTKeyDouble::const_iterator i = pred_peak_bins.begin(); i != pred_peak_bins.end(); i++) {
        if (i->second < closest_pred_peak_bin){
            closest_pred_peak_bin = i->second;
        }
    }

    // ok so we have data and we have prediction
    // before we can calculate our LLH, we need to figure our the timing

    size_t llh_bins_from_peak = closest_pred_peak_bin;
    size_t total_llh_time_bins = 100;

    double total_nllh = 0.0;

    for (I3MapPMTKeyVectorDouble::const_iterator i = pred_yields->begin(); i != pred_yields->end(); i++) {
        // looping over each pmt key in our pred map between CCMPMTKey and std::vector<double>
        // now loop over the time bins in our pred
        for (size_t time_bin_it = (pred_peak_bins[i->first] - llh_bins_from_peak); time_bin_it < (pred_peak_bins[i->first] - llh_bins_from_peak + total_llh_time_bins); time_bin_it++) {
            size_t corresponding_data_bin = data_peak_bins[i->first] - llh_bins_from_peak + (time_bin_it - (pred_peak_bins[i->first] - llh_bins_from_peak)); 
            double k = data.at(i->first).at(corresponding_data_bin);
            double mu = (i->second).at(time_bin_it) * pinned_nuissance[i->first];  
            double sigma_squared = pred_yields_squared->at(i->first).at(time_bin_it) * std::pow(pinned_nuissance[i->first], 2.0);
            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
            
            // save some things for debugging
            // add this key to our maps between data and pred for debugging
            if (debug_data.find(i->first) == debug_data.end()) {
                debug_data[i->first] = std::vector<double> (total_llh_time_bins, 0.0);
            }
            if (debug_pred.find(i->first) == debug_pred.end()) {
                debug_pred[i->first] = std::vector<double> (total_llh_time_bins, 0.0);
            }
            if (debug_sigma2.find(i->first) == debug_sigma2.end()) {
                debug_sigma2[i->first] = std::vector<double> (total_llh_time_bins, 0.0);
            }

            // now save our data and pred
            debug_data.at(i->first).at(time_bin_it + llh_bins_from_peak - pred_peak_bins[i->first]) = k;
            debug_pred.at(i->first).at(time_bin_it + llh_bins_from_peak - pred_peak_bins[i->first]) = mu;
            debug_sigma2.at(i->first).at(time_bin_it + llh_bins_from_peak - pred_peak_bins[i->first]) = sigma_squared;
        }
    }

    return total_nllh;
}


