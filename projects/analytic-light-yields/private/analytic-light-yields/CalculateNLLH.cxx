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

CalculateNLLH::CalculateNLLH() {}

void CalculateNLLH::GrabData(I3FramePtr data_frame){

    // grab our data out of the frame
    data = data_frame->Get<I3MapPMTKeyVectorDouble>("AccumulatedEventsMap");
    n_data_events = data_frame->Get<I3Double>("TotalEventsPastCuts").value;
    //std::cout << "n_data_events = " << n_data_events << std::endl;
    grabbed_data = true;
}

I3MapPMTKeyDouble GetPeakBin(boost::shared_ptr<I3MapPMTKeyVectorDouble> data){

    // ok let's get peak bin in data
    I3MapPMTKeyDouble data_peak_bins;
    for (I3MapPMTKeyVectorDouble::const_iterator i = data->begin(); i != data->end(); i++) {
        size_t peak_bin;
        double peak_val = 0.0;
        std::vector<double> this_pmt_data = i->second;
        
        // now loop over wf in this pmt
        for (size_t j = 0; j < this_pmt_data.size(); j++){
            if (this_pmt_data.at(j) > peak_val){
                peak_val = this_pmt_data.at(j);
                peak_bin = j;
            }
        }

        // now save
        data_peak_bins[i->first] = peak_bin;
    }

    return data_peak_bins;
}

I3MapPMTKeyDouble GetPeakBin(I3MapPMTKeyVectorDouble data){

    // ok let's get peak bin in data
    I3MapPMTKeyDouble data_peak_bins;
    for (I3MapPMTKeyVectorDouble::const_iterator i = data.begin(); i != data.end(); i++) {
        size_t peak_bin;
        double peak_val = 0.0;
        std::vector<double> this_pmt_data = i->second;
        
        // now loop over wf in this pmt
        for (size_t j = 0; j < this_pmt_data.size(); j++){
            if (this_pmt_data.at(j) > peak_val){
                peak_val = this_pmt_data.at(j);
                peak_bin = j;
            }
        }

        // now save
        data_peak_bins[i->first] = peak_bin;
    }

    return data_peak_bins;
}

I3MapPMTKeyDouble PinNuisance(boost::shared_ptr<I3MapPMTKeyVectorDouble> pred, I3MapPMTKeyVectorDouble data){

    // ok let's get peak height in data and pred for each pmt, then figure out scaling factor to make them match 
    I3MapPMTKeyDouble data_peak_height;
    I3MapPMTKeyDouble pred_peak_height;

    for (I3MapPMTKeyVectorDouble::const_iterator i = pred->begin(); i != pred->end(); i++) {
        double pred_peak_val = 0.0;
        std::vector<double> this_pmt_pred = i->second;
        
        // now loop over wf in this pmt
        for (size_t j = 0; j < this_pmt_pred.size(); j++){
            if (this_pmt_pred.at(j) > pred_peak_val){
                pred_peak_val = this_pmt_pred.at(j);
            }
        }

        // now save
        pred_peak_height[i->first] = pred_peak_val;
    }
    
    for (I3MapPMTKeyVectorDouble::const_iterator i = data.begin(); i != data.end(); i++) {
        double data_peak_val = 0.0;
        std::vector<double> this_pmt_data = i->second;
        
        // now loop over wf in this pmt
        for (size_t j = 0; j < this_pmt_data.size(); j++){
            if (this_pmt_data.at(j) > data_peak_val){
                data_peak_val = this_pmt_data.at(j);
            }
        }

        // now save
        data_peak_height[i->first] = data_peak_val;
    }

    // now make our scaling map
    I3MapPMTKeyDouble scaling;
    for (I3MapPMTKeyDouble::const_iterator j = data_peak_height.begin(); j != data_peak_height.end(); j++){
        double this_scaling_factor = (j->second / pred_peak_height[j->first]);
        scaling[j->first] = this_scaling_factor;
    }


    return scaling;
}

std::vector<I3MapPMTKeyVectorDouble> CalculateNLLH::DebugDatavsPred(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, I3FramePtr data_frame,
                                                                    std::vector<CCMPMTKey> keys_to_fit, double single_pmt_norm, double single_pmt_offset, double light_time_offset){
    GetNLLHDerivative(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
    
    std::vector<I3MapPMTKeyVectorDouble> final_debug_info;
    final_debug_info.push_back(debug_data);
    final_debug_info.push_back(debug_pred);
    final_debug_info.push_back(debug_sigma2);
    return final_debug_info;

}

I3Vector<double> CalculateNLLH::GetNLLHDerivative(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, I3FramePtr data_frame,
                                                  std::vector<CCMPMTKey> keys_to_fit, double single_pmt_norm, double single_pmt_offset, double light_time_offset){
    // check to see if we've grabbed our data
    if (grabbed_data == false){
        GrabData(data_frame);
    }

    // let's check to see if we've initialized our GenerateExpectation constructor
    if (gen_expectation == nullptr){
        gen_expectation = std::make_shared<GenerateExpectation> ();
    }
    
    // let's grab our expectation
    std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>> pred = gen_expectation->GetExpectationWithDerivs(analytic_light_yield_setup, geo_frame, 
                                                                                                                single_pmt_norm, single_pmt_offset, keys_to_fit, light_time_offset);
    
    // unpack into yields and yields^2
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    // and unpack derivs of yields wrt light profile variables
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_Rs = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared_Rs = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_Rt = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared_Rt = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_tau_s = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared_tau_s = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_tau_t = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared_tau_t = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_tau_TPB = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared_tau_TPB = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_time_offset = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared_time_offset = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    std::tie(pred_yields, pred_yields_squared,
             pred_yields_Rs, pred_yields_squared_Rs,
             pred_yields_Rt, pred_yields_squared_Rt, 
             pred_yields_tau_s, pred_yields_squared_tau_s,
             pred_yields_tau_t, pred_yields_squared_tau_t, 
             pred_yields_tau_TPB, pred_yields_squared_tau_TPB,
             pred_yields_time_offset, pred_yields_squared_time_offset) = pred;

    // and grab number of simulation events
    double n_simulation_events = analytic_light_yield_setup.n_sodium_events;

    // let grab our peak bins in data and in pred
    I3MapPMTKeyDouble data_peak_bins = GetPeakBin(data);
    I3MapPMTKeyDouble pred_peak_bins = GetPeakBin(pred_yields);

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
    double total_dRs = 0.0;
    double total_dRt = 0.0;
    double total_dtau_s = 0.0;
    double total_dtau_t = 0.0;
    double total_dtau_TPB = 0.0;
    double total_dScaling = 0.0;
    double total_dOffset = 0.0;
    double total_dtime_offset = 0.0;

    for (I3MapPMTKeyVectorDouble::const_iterator i = pred_yields->begin(); i != pred_yields->end(); i++) {
        // looping over each pmt key in our pred map between CCMPMTKey and std::vector<double>
        // now loop over the time bins in our pred
        for (size_t time_bin_it = (pred_peak_bins[i->first] - llh_bins_from_peak); time_bin_it < (pred_peak_bins[i->first] - llh_bins_from_peak + total_llh_time_bins); time_bin_it++){
        //for (size_t time_bin_it = 0; time_bin_it < total_llh_time_bins; time_bin_it++){
            size_t corresponding_data_bin = data_peak_bins[i->first] - llh_bins_from_peak + (time_bin_it - (pred_peak_bins[i->first] - llh_bins_from_peak)); 
            double k = data.at(i->first).at(corresponding_data_bin);
            double mu = (i->second).at(time_bin_it); 
            double sigma_squared = pred_yields_squared->at(i->first).at(time_bin_it);
            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);

            // now let's get our derivs
            double dw_dRs = pred_yields_Rs->at(i->first).at(time_bin_it);
            double dw2_dRs = pred_yields_squared_Rs->at(i->first).at(time_bin_it);
            total_dRs += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dRs, dw2_dRs);

            double dw_dRt = pred_yields_Rt->at(i->first).at(time_bin_it);
            double dw2_dRt = pred_yields_squared_Rt->at(i->first).at(time_bin_it);
            total_dRt += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dRt, dw2_dRt);

            double dw_dtau_s = pred_yields_tau_s->at(i->first).at(time_bin_it);
            double dw2_dtau_s = pred_yields_squared_tau_s->at(i->first).at(time_bin_it);
            total_dtau_s += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dtau_s, dw2_dtau_s);

            double dw_dtau_t = pred_yields_tau_t->at(i->first).at(time_bin_it);
            double dw2_dtau_t = pred_yields_squared_tau_t->at(i->first).at(time_bin_it);
            total_dtau_t += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dtau_t, dw2_dtau_t);

            double dw_dtau_TPB = pred_yields_tau_TPB->at(i->first).at(time_bin_it);
            double dw2_dtau_TPB = pred_yields_squared_tau_TPB->at(i->first).at(time_bin_it);
            total_dtau_TPB += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dtau_TPB, dw2_dtau_TPB);

            double dw_dScaling = (mu - (n_simulation_events * single_pmt_offset)) / single_pmt_norm;
            double dw2_dScaling = 2 * mu * dw_dScaling;
            total_dScaling += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dScaling, dw2_dScaling);

            double dw_dOffset = n_simulation_events;
            double dw2_dOffset = 2.0 * mu;
            total_dOffset += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dOffset, dw2_dOffset);
            
            double dw_dtime_offset = pred_yields_time_offset->at(i->first).at(time_bin_it);
            double dw2_dtime_offset = pred_yields_squared_time_offset->at(i->first).at(time_bin_it);
            total_dtime_offset += MCLLH::LEffDeriv()(k, mu, sigma_squared, dw_dtime_offset, dw2_dtime_offset);

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
            debug_data.at(i->first).at(time_bin_it - (pred_peak_bins[i->first] - llh_bins_from_peak)) = k;
            debug_pred.at(i->first).at(time_bin_it - (pred_peak_bins[i->first] - llh_bins_from_peak)) = mu;
            debug_sigma2.at(i->first).at(time_bin_it - (pred_peak_bins[i->first] - llh_bins_from_peak)) = sigma_squared;
        }
    }

    // now make our vector to return everybody 
    I3Vector<double> nllh_and_derivs;
    nllh_and_derivs.push_back(total_nllh);
    nllh_and_derivs.push_back(total_dRs);
    nllh_and_derivs.push_back(total_dRt);
    nllh_and_derivs.push_back(total_dtau_s);
    nllh_and_derivs.push_back(total_dtau_t);
    nllh_and_derivs.push_back(total_dtau_TPB);
    nllh_and_derivs.push_back(total_dScaling);
    //nllh_and_derivs.push_back(total_dOffset);
    nllh_and_derivs.push_back(total_dtime_offset);

    return nllh_and_derivs;

}



double CalculateNLLH::GetNLLH(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, std::vector<CCMPMTKey> keys_to_fit, I3FramePtr data_frame, double light_time_offset){

    // check to see if we've grabbed our data
    if (grabbed_data == false){
        GrabData(data_frame);
    }

    // let's check to see if we've initialized our GenerateExpectation constructor
    if (gen_expectation == nullptr){
        gen_expectation = std::make_shared<GenerateExpectation> ();
    }
    
    // let's grab our expectation
    std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>> pred = gen_expectation->GetExpectation(analytic_light_yield_setup, geo_frame,
                                                                                                            keys_to_fit, light_time_offset);

    // unpack into yields and yields^2
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    std::tie(pred_yields, pred_yields_squared) = pred;
    

    // let grab our peak bins in data and in pred
    I3MapPMTKeyDouble data_peak_bins = GetPeakBin(data);
    I3MapPMTKeyDouble pred_peak_bins = GetPeakBin(pred_yields);

    // let's also grab the nuissance params...pin them so they scale data and pred peaks to be at the same height
    I3MapPMTKeyDouble pinned_nuissance = PinNuisance(pred_yields, data);

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
        for (size_t time_bin_it = (pred_peak_bins[i->first] - llh_bins_from_peak); time_bin_it < (pred_peak_bins[i->first] - llh_bins_from_peak + total_llh_time_bins); time_bin_it++){
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







