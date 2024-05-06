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
    std::cout << "n_data_events = " << n_data_events << std::endl;
    grabbed_data = true;
}

I3Vector<double> CalculateNLLH::GetNLLHDerivative(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
                                                  std::vector<double> nuisance_params, I3FramePtr data_frame, size_t time_bin_offset){

    // check to see if we've grabbed our data
    if (grabbed_data == false){
        GrabData(data_frame);
    }

    // let's check to see if we've initialized our GenerateExpectation constructor
    if (gen_expectation == nullptr){
        gen_expectation = new GenerateExpectation();
    }
    // let's grab our expectation
    std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
               boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>> pred = gen_expectation->GetExpectationWithDerivs(analytic_light_yield_setup, geo_frame);
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
    std::tie(pred_yields, pred_yields_squared,
             pred_yields_Rs, pred_yields_squared_Rs,
             pred_yields_Rt, pred_yields_squared_Rt, 
             pred_yields_tau_s, pred_yields_squared_tau_s,
             pred_yields_tau_t, pred_yields_squared_tau_t, 
             pred_yields_tau_TPB, pred_yields_squared_tau_TPB) = pred;

    // and grab number of simulation events
    double n_simulation_events = analytic_light_yield_setup.n_sodium_events;

    // ok so we have data and we have prediction
    // before we can calculate our LLH, we need to figure our the timing
    // let's try aligning on event start time but allowing a time offset to float
    // then we will calculate our llh from 6nsec after event start until 70nsec after event start

    size_t llh_bins_from_start = 5;
    size_t total_llh_time_bins = 35;

    size_t corresponding_data_bin;
    double mu;
    double sigma_squared;
    double k;
    double DmuDtheta;
    double Dsigma_squaredDtheta;

    double total_nllh = 0.0;
    double total_dRs = 0.0;
    double total_dRt = 0.0;
    double total_dtau_s = 0.0;
    double total_dtau_t = 0.0;
    double total_dtau_TPB = 0.0;

    // let's also set up our vector for derivs of nuisance_params
    std::vector<double> nuisance_param_derivs(nuisance_params.size() , 0.0);

    size_t nuisance_counter = 0;
    for (I3MapPMTKeyVectorDouble::const_iterator i = pred_yields->begin(); i != pred_yields->end(); i++) {
        // looping over each pmt key in our pred map between CCMPMTKey and std::vector<double>
        // now loop over the time bins in our pred
        for (size_t time_bin_it = llh_bins_from_start; time_bin_it < (total_llh_time_bins + llh_bins_from_start); time_bin_it++){
            corresponding_data_bin = time_bin_it + time_bin_offset;
            k = data.at(i->first).at(corresponding_data_bin);
            mu = (i->second).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            sigma_squared = pred_yields_squared->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
            // now let's get our derivs
            DmuDtheta = pred_yields_Rs->at(i->first).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            Dsigma_squaredDtheta = pred_yields_squared_Rs->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_dRs += MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
            DmuDtheta = pred_yields_Rt->at(i->first).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            Dsigma_squaredDtheta = pred_yields_squared_Rt->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_dRt += MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
            DmuDtheta = pred_yields_tau_s->at(i->first).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            Dsigma_squaredDtheta = pred_yields_squared_tau_s->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_dtau_s += MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
            DmuDtheta = pred_yields_tau_t->at(i->first).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            Dsigma_squaredDtheta = pred_yields_squared_tau_t->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_dtau_t += MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
            DmuDtheta = pred_yields_tau_TPB->at(i->first).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            Dsigma_squaredDtheta = pred_yields_squared_tau_TPB->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_dtau_TPB += MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
            // and can't forget about the nuissance params! luckily the derivs are easy
            DmuDtheta = mu / nuisance_params[nuisance_counter];
            Dsigma_squaredDtheta = 2 * nuisance_params[nuisance_counter] * pred_yields_squared->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events), 2.0);
            nuisance_param_derivs[nuisance_counter] += MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
        }
        nuisance_counter += 1;
    }

    // now make our vector to return everybody 
    I3Vector<double> nllh_and_derivs;
    nllh_and_derivs.push_back(total_nllh);
    nllh_and_derivs.push_back(total_dRs);
    nllh_and_derivs.push_back(total_dRt);
    nllh_and_derivs.push_back(total_dtau_s);
    nllh_and_derivs.push_back(total_dtau_t);
    nllh_and_derivs.push_back(total_dtau_TPB);
    for (size_t k = 0; k < nuisance_param_derivs.size(); k++){
        nllh_and_derivs.push_back(nuisance_param_derivs.at(k));
    }
    
    return nllh_and_derivs;

}


double CalculateNLLH::GetNLLH(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, std::vector<double> nuisance_params, I3FramePtr data_frame, size_t time_bin_offset){

    // check to see if we've grabbed our data
    if (grabbed_data == false){
        GrabData(data_frame);
    }
    
    // let's check to see if we've initialized our GenerateExpectation constructor
    if (gen_expectation == nullptr){
        gen_expectation = new GenerateExpectation();
    }
    // let's grab our expectation
    std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>> pred = gen_expectation->GetExpectation(analytic_light_yield_setup, geo_frame);
    // unpack into yields and yields^2
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred_yields_squared = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    std::tie(pred_yields, pred_yields_squared) = pred;

    // and grab number of simulation events
    double n_simulation_events = analytic_light_yield_setup.n_sodium_events;

    // ok so we have data and we have prediction
    // before we can calculate our LLH, we need to figure our the timing
    // let's try aligning on event start time but allowing a time offset to float
    // then we will calculate our llh from 6nsec after event start until 70nsec after event start

    size_t llh_bins_from_start = 5;
    size_t total_llh_time_bins = 35;

    size_t corresponding_data_bin;
    double mu;
    double sigma_squared;
    double k;
    double total_nllh = 0.0;

    size_t nuisance_counter = 0;
    for (I3MapPMTKeyVectorDouble::const_iterator i = pred_yields->begin(); i != pred_yields->end(); i++) {
        // looping over each pmt key in our pred map between CCMPMTKey and std::vector<double>
        // now loop over the time bins in our pred
        for (size_t time_bin_it = llh_bins_from_start; time_bin_it < (total_llh_time_bins + llh_bins_from_start); time_bin_it++){
            corresponding_data_bin = time_bin_it + time_bin_offset;
            k = data.at(i->first).at(corresponding_data_bin);
            mu = (i->second).at(time_bin_it) * (n_data_events / n_simulation_events) * nuisance_params[nuisance_counter];
            sigma_squared = pred_yields_squared->at(i->first).at(time_bin_it) * std::pow((n_data_events / n_simulation_events) * nuisance_params[nuisance_counter], 2.0);
            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
        }
        nuisance_counter += 1;
    }

    return total_nllh;

}







