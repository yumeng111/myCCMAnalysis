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
#include <analytic-light-yields/GenerateExpectation.h>


GenerateExpectation::GenerateExpectation() {}

void GenerateExpectation::GetSodiumVertices(size_t const & n_events_to_simulate, double const & z_position){
    sodium_events_constructor = std::make_shared<SodiumVertexDistribution> ();
    event_vertices = sodium_events_constructor->GetEventVertices(n_events_to_simulate, z_position);
}

void GenerateExpectation::GetYieldsAndOffsets(I3FramePtr geo_frame, double const & uv_absorption, std::vector<CCMPMTKey> keys_to_fit){
    yields_and_offset_constructor = std::make_shared<YieldsPerPMT> ();
    // now loop over events and get map between CCMPMTKey and PhotonYieldSummarySeries
    for (size_t sodium_it = 0; sodium_it < event_vertices->size(); sodium_it++){
        boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_per_event = yields_and_offset_constructor->GetAllYields(event_vertices->at(sodium_it), geo_frame, uv_absorption, keys_to_fit);
        yields_per_pmt_per_event.push_back(yields_per_event);
    }
}

I3Vector<double> GenerateExpectation::LightProfile(double const & Rs, double const & Rt, double const & tau_s, double const & tau_t, double const & tau_rec, double const & tau_TPB,
                                                    AnalyticLightYieldGenerator::LArLightProfileType const & light_profile_type, double light_time_offset){

    if (LAr_scintillation_light_constructor == nullptr){
        LAr_scintillation_light_constructor = std::make_shared<LArScintillationLightProfile> ();
    }

    I3Vector<double> light_profile;

    if (light_profile_type == AnalyticLightYieldGenerator::LArLightProfileType::Simplified){
        light_profile = LAr_scintillation_light_constructor->GetSimplifiedLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, light_time_offset);
    }
    else {
        light_profile = LAr_scintillation_light_constructor->GetFullLightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_time_offset);
    }

    return light_profile;
}

I3Vector<double> GenerateExpectation::DLightProfile(double const & Rs, double const & Rt, double const & tau_s, double const & tau_t, double const & tau_TPB, 
                                                    std::string deriv_variable, double light_time_offset){

    if (LAr_scintillation_light_constructor == nullptr){
        LAr_scintillation_light_constructor = std::make_shared<LArScintillationLightProfile> ();
    }

    I3Vector<double> light_profile;
    light_profile = LAr_scintillation_light_constructor->GetSimplifiedLightProfileDeriv(Rs, Rt, tau_s, tau_t, tau_TPB, deriv_variable, light_time_offset);

    return light_profile;
}

std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>> GenerateExpectation::GetExpectation(AnalyticLightYieldGenerator analytic_light_yield_setup,
                                                                                                                                       I3FramePtr geo_frame, std::vector<CCMPMTKey> keys_to_fit,
                                                                                                                                       double light_time_offset){

    // let's parse the info necessary to get our expectation
    double Rs = analytic_light_yield_setup.Rs;
    double Rt = analytic_light_yield_setup.Rt;
    double tau_s = analytic_light_yield_setup.tau_s;
    double tau_t = analytic_light_yield_setup.tau_t;
    double tau_rec = analytic_light_yield_setup.tau_rec;
    double tau_TPB = analytic_light_yield_setup.tau_TPB;
    double uv_absorption = analytic_light_yield_setup.uv_absorption;
    double normalization = analytic_light_yield_setup.normalization;
    double n_sodium_events = analytic_light_yield_setup.n_sodium_events;
    double z_offset = analytic_light_yield_setup.z_offset;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = analytic_light_yield_setup.light_profile_type;

    // check that we made our sodium event vertices
    if (sodium_events_constructor == nullptr){
        GetSodiumVertices(n_sodium_events, z_offset);
    }

    // now let's check if we get our yields + time offsets
    if (yields_and_offset_constructor == nullptr){
        GetYieldsAndOffsets(geo_frame, uv_absorption, keys_to_fit);
    }

    // now grab our light profile!
    I3Vector<double> LAr_light_profile = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_profile_type, light_time_offset);

    // so let's recap where we are
    // 1) we generated sodium event vertices
    // 2) we took those vertices and calculated yields and time offsets in PMTs
    // 3) we calculated our light profile
    // so all that's left to do is take our yields and apply to LAr light profile then bin charge
    // after binning, we are going to be ready to compare to data!

    // make map of final binned yields to return
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();

    double light_profile_start_time = 0.0;
    size_t n_light_profile_time_bins = LAr_light_profile.size();
    I3Vector<double> light_profile_times;
    for (size_t i = 0; i < n_light_profile_time_bins; i++){
        light_profile_times.push_back(light_profile_start_time + 2.0 * (double)i);
    }

    // let's also make the time binning for our expectation
    double expectation_start_time = 0.0;
    size_t n_expectation_time_bins = n_light_profile_time_bins + 10; // making it a bit longer than light profile to account for time shifts
    I3Vector<double> expectation_times;
    for (size_t i = 0; i < n_expectation_time_bins; i++){
        expectation_times.push_back(expectation_start_time + 2.0 * (double)i);
    }

    double time_offset;
    double relative_yields;
    double light_time;
    double light_val;
    size_t bin_idx;
    for (size_t event_it = 0; event_it < yields_per_pmt_per_event.size(); event_it ++){
        boost::shared_ptr<PhotonYieldSummarySeriesMap> this_event_yields_per_pmt = yields_per_pmt_per_event.at(event_it); 
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields = boost::make_shared<I3MapPMTKeyVectorDouble> ();

        for (PhotonYieldSummarySeriesMap::const_iterator i = this_event_yields_per_pmt->begin(); i != this_event_yields_per_pmt->end(); i++) {
            // check if this pmt is in this_event_binned_yields -- if not, add it
            if (this_event_binned_yields->find(i->first) == this_event_binned_yields->end()) {
                (*this_event_binned_yields)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }

            // now let's loop over all the PhotonYieldSummary for this PMT
            for(PhotonYieldSummary const & this_yield: i->second) {
                time_offset = this_yield.time;
                relative_yields = this_yield.yield;

                // now apply time offsets and yields to our light profile
                for (size_t light_time_it = 0; light_time_it < n_light_profile_time_bins; light_time_it++){
                    light_time = light_profile_times.at(light_time_it) + time_offset;
                    light_val = LAr_light_profile.at(light_time_it) * relative_yields;

                    // now binning!
                    bin_idx = (size_t) (light_time / 2.0);
                    (*this_event_binned_yields)[i->first].at(bin_idx) += light_val;
                }
            }
        }

        // ok so for this event we have finished binning the expected light yield
        // let's add this event yields and yields^2 to our final lists
        for (I3MapPMTKeyVectorDouble::const_iterator j = this_event_binned_yields->begin(); j != this_event_binned_yields->end(); j++) {
            // check if this pmt key is in final_binned_yields_
            if (final_binned_yields_->find(j->first) == final_binned_yields_->end()) {
                (*final_binned_yields_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_
            if (final_binned_squared_yields_->find(j->first) == final_binned_squared_yields_->end()) {
                (*final_binned_squared_yields_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }

            // now loop over the yields in every time bin and square
            for (size_t time_bin_it = 0; time_bin_it < this_event_binned_yields->at(j->first).size(); time_bin_it++){
                final_binned_yields_->at(j->first).at(time_bin_it) += (this_event_binned_yields->at(j->first).at(time_bin_it));
                final_binned_squared_yields_->at(j->first).at(time_bin_it) += std::pow(this_event_binned_yields->at(j->first).at(time_bin_it), 2.0);
            }
        }
    }

    return std::make_tuple(final_binned_yields_, final_binned_squared_yields_);

}


std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>, 
           boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
           boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
           boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
           boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
           boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>,
           boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>>
           GenerateExpectation::GetExpectationWithDerivs(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
                                                         double const & single_pmt_norm, double const & single_pmt_offset, std::vector<CCMPMTKey> keys_to_fit, double light_time_offset){
    // let's parse the info necessary to get our expectation
    double Rs = analytic_light_yield_setup.Rs;
    double Rt = analytic_light_yield_setup.Rt;
    double tau_s = analytic_light_yield_setup.tau_s;
    double tau_t = analytic_light_yield_setup.tau_t;
    double tau_rec = analytic_light_yield_setup.tau_rec;
    double tau_TPB = analytic_light_yield_setup.tau_TPB;
    double uv_absorption = analytic_light_yield_setup.uv_absorption;
    double normalization = analytic_light_yield_setup.normalization;
    double n_sodium_events = analytic_light_yield_setup.n_sodium_events;
    double z_offset = analytic_light_yield_setup.z_offset;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = analytic_light_yield_setup.light_profile_type;

    // check that we made our sodium event vertices
    if (sodium_events_constructor == nullptr){
        GetSodiumVertices(n_sodium_events, z_offset);
    }

    // now let's check if we get our yields + time offsets
    if (yields_and_offset_constructor == nullptr){
        GetYieldsAndOffsets(geo_frame, uv_absorption, keys_to_fit);
    }

    // now grab our light profile!
    I3Vector<double> LAr_light_profile = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_profile_type, light_time_offset);
    I3Vector<double> DLAr_light_profile_dRs = DLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, "Rs", light_time_offset);
    I3Vector<double> DLAr_light_profile_dRt = DLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, "Rt", light_time_offset);
    I3Vector<double> DLAr_light_profile_dtau_s = DLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, "tau_s", light_time_offset);
    I3Vector<double> DLAr_light_profile_dtau_t = DLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, "tau_t", light_time_offset);
    I3Vector<double> DLAr_light_profile_dtau_TPB = DLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, "tau_TPB", light_time_offset);
    I3Vector<double> DLAr_light_profile_dtime_offset = DLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, "time_offset", light_time_offset);

    // so let's recap where we are
    // 1) we generated sodium event vertices
    // 2) we took those vertices and calculated yields and time offsets in PMTs
    // 3) we calculated our light profile
    // so all that's left to do is take our yields and apply to LAr light profile then bin charge
    // after binning, we are going to be ready to compare to data!

    // make map of final binned yields to return
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_deriv_Rs_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_deriv_Rs_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_deriv_Rt_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_deriv_Rt_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_deriv_tau_s_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_deriv_tau_s_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_deriv_tau_t_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_deriv_tau_t_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_deriv_tau_TPB_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_deriv_tau_TPB_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_deriv_time_offset_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_squared_yields_deriv_time_offset_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();

    double light_profile_start_time = 0.0;
    size_t n_light_profile_time_bins = LAr_light_profile.size();
    I3Vector<double> light_profile_times;
    for (size_t i = 0; i < n_light_profile_time_bins; i++){
        light_profile_times.push_back(light_profile_start_time + 2.0 * (double)i);
    }

    // let's also make the time binning for our expectation
    double expectation_start_time = 0.0;
    size_t n_expectation_time_bins = n_light_profile_time_bins + 10; // making it a bit longer than light profile to account for time shifts
    I3Vector<double> expectation_times;
    for (size_t i = 0; i < n_expectation_time_bins; i++){
        expectation_times.push_back(expectation_start_time + 2.0 * (double)i);
    }

    for (size_t event_it = 0; event_it < yields_per_pmt_per_event.size(); event_it ++){
        boost::shared_ptr<PhotonYieldSummarySeriesMap> this_event_yields_per_pmt = yields_per_pmt_per_event.at(event_it); 
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields = boost::make_shared<I3MapPMTKeyVectorDouble> ();
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields_deriv_Rs = boost::make_shared<I3MapPMTKeyVectorDouble> ();
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields_deriv_Rt = boost::make_shared<I3MapPMTKeyVectorDouble> ();
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields_deriv_tau_s = boost::make_shared<I3MapPMTKeyVectorDouble> ();
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields_deriv_tau_t = boost::make_shared<I3MapPMTKeyVectorDouble> ();
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields_deriv_tau_TPB = boost::make_shared<I3MapPMTKeyVectorDouble> ();
        boost::shared_ptr<I3MapPMTKeyVectorDouble> this_event_binned_yields_deriv_time_offset = boost::make_shared<I3MapPMTKeyVectorDouble> ();

        for (PhotonYieldSummarySeriesMap::const_iterator i = this_event_yields_per_pmt->begin(); i != this_event_yields_per_pmt->end(); i++) {
            // check if this pmt is in this_event_binned_yields -- if not, add it
            if (this_event_binned_yields->find(i->first) == this_event_binned_yields->end()) {
                (*this_event_binned_yields)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt is in this_event_binned_yields_deriv_Rs -- if not, add it
            if (this_event_binned_yields_deriv_Rs->find(i->first) == this_event_binned_yields_deriv_Rs->end()) {
                (*this_event_binned_yields_deriv_Rs)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt is in this_event_binned_yields_deriv_Rt -- if not, add it
            if (this_event_binned_yields_deriv_Rt->find(i->first) == this_event_binned_yields_deriv_Rt->end()) {
                (*this_event_binned_yields_deriv_Rt)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt is in this_event_binned_yields_deriv_tau_s -- if not, add it
            if (this_event_binned_yields_deriv_tau_s->find(i->first) == this_event_binned_yields_deriv_tau_s->end()) {
                (*this_event_binned_yields_deriv_tau_s)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt is in this_event_binned_yields_deriv_tau_t -- if not, add it
            if (this_event_binned_yields_deriv_tau_t->find(i->first) == this_event_binned_yields_deriv_tau_t->end()) {
                (*this_event_binned_yields_deriv_tau_t)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt is in this_event_binned_yields_deriv_tau_TPB -- if not, add it
            if (this_event_binned_yields_deriv_tau_TPB->find(i->first) == this_event_binned_yields_deriv_tau_TPB->end()) {
                (*this_event_binned_yields_deriv_tau_TPB)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt is in this_event_binned_yields_deriv_time_offset -- if not, add it
            if (this_event_binned_yields_deriv_time_offset->find(i->first) == this_event_binned_yields_deriv_time_offset->end()) {
                (*this_event_binned_yields_deriv_time_offset)[i->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }

            // now let's loop over all the PhotonYieldSummary for this PMT
            for(PhotonYieldSummary const & this_yield: i->second) {
                double time_offset = this_yield.time;
                double relative_yields = this_yield.yield;

                // now apply time offsets and yields to our light profile
                for (size_t light_time_it = 0; light_time_it < n_light_profile_time_bins; light_time_it++){
                    // now binning!
                    // let's get time bin
                    double light_time = light_profile_times.at(light_time_it) + time_offset;
                    size_t bin_idx = (size_t) (light_time / 2.0);
                    
                    // now light profile
                    double light_val = LAr_light_profile.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields)[i->first].at(bin_idx) += light_val;
                
                    // now do again for derivs
                    double light_val_dRs = DLAr_light_profile_dRs.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields_deriv_Rs)[i->first].at(bin_idx) += light_val_dRs;

                    double light_val_dRt = DLAr_light_profile_dRt.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields_deriv_Rt)[i->first].at(bin_idx) += light_val_dRt;

                    double light_val_dtau_s = DLAr_light_profile_dtau_s.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields_deriv_tau_s)[i->first].at(bin_idx) += light_val_dtau_s;

                    double light_val_dtau_t = DLAr_light_profile_dtau_t.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields_deriv_tau_t)[i->first].at(bin_idx) += light_val_dtau_t;

                    double light_val_dtau_TPB = DLAr_light_profile_dtau_TPB.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields_deriv_tau_TPB)[i->first].at(bin_idx) += light_val_dtau_TPB;
                    
                    double light_val_dtime_offset = DLAr_light_profile_dtime_offset.at(light_time_it) * relative_yields;
                    (*this_event_binned_yields_deriv_time_offset)[i->first].at(bin_idx) += light_val_dtime_offset;

                }
            }
        }

        // ok so for this event we have finished binning the expected light yield
        // let's add this event yields and yields^2 to our final lists
        for (I3MapPMTKeyVectorDouble::const_iterator j = this_event_binned_yields->begin(); j != this_event_binned_yields->end(); j++) {
            // check if this pmt key is in final_binned_yields_
            if (final_binned_yields_->find(j->first) == final_binned_yields_->end()) {
                (*final_binned_yields_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_
            if (final_binned_squared_yields_->find(j->first) == final_binned_squared_yields_->end()) {
                (*final_binned_squared_yields_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_yields_deriv_Rs_
            if (final_binned_yields_deriv_Rs_->find(j->first) == final_binned_yields_deriv_Rs_->end()) {
                (*final_binned_yields_deriv_Rs_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_deriv_Rs_
            if (final_binned_squared_yields_deriv_Rs_->find(j->first) == final_binned_squared_yields_deriv_Rs_->end()) {
                (*final_binned_squared_yields_deriv_Rs_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_yields_deriv_Rt_
            if (final_binned_yields_deriv_Rt_->find(j->first) == final_binned_yields_deriv_Rt_->end()) {
                (*final_binned_yields_deriv_Rt_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_deriv_Rt_
            if (final_binned_squared_yields_deriv_Rt_->find(j->first) == final_binned_squared_yields_deriv_Rt_->end()) {
                (*final_binned_squared_yields_deriv_Rt_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_yields_deriv_tau_s_
            if (final_binned_yields_deriv_tau_s_->find(j->first) == final_binned_yields_deriv_tau_s_->end()) {
                (*final_binned_yields_deriv_tau_s_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_deriv_tau_s_
            if (final_binned_squared_yields_deriv_tau_s_->find(j->first) == final_binned_squared_yields_deriv_tau_s_->end()) {
                (*final_binned_squared_yields_deriv_tau_s_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_yields_deriv_tau_t_
            if (final_binned_yields_deriv_tau_t_->find(j->first) == final_binned_yields_deriv_tau_t_->end()) {
                (*final_binned_yields_deriv_tau_t_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_deriv_tau_t_
            if (final_binned_squared_yields_deriv_tau_t_->find(j->first) == final_binned_squared_yields_deriv_tau_t_->end()) {
                (*final_binned_squared_yields_deriv_tau_t_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_yields_deriv_tau_TPB_
            if (final_binned_yields_deriv_tau_TPB_->find(j->first) == final_binned_yields_deriv_tau_TPB_->end()) {
                (*final_binned_yields_deriv_tau_TPB_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_deriv_tau_TPB_
            if (final_binned_squared_yields_deriv_tau_TPB_->find(j->first) == final_binned_squared_yields_deriv_tau_TPB_->end()) {
                (*final_binned_squared_yields_deriv_tau_TPB_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_yields_deriv_time_offset_
            if (final_binned_yields_deriv_time_offset_->find(j->first) == final_binned_yields_deriv_time_offset_->end()) {
                (*final_binned_yields_deriv_time_offset_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }
            // check if this pmt key is in final_binned_squared_yields_deriv_time_offset_
            if (final_binned_squared_yields_deriv_time_offset_->find(j->first) == final_binned_squared_yields_deriv_time_offset_->end()) {
                (*final_binned_squared_yields_deriv_time_offset_)[j->first] = std::vector<double> (n_expectation_time_bins, 0.0);
            }


            // now loop over the yields in every time bin and square
            for (size_t time_bin_it = 0; time_bin_it < this_event_binned_yields->at(j->first).size(); time_bin_it++){
                double w = this_event_binned_yields->at(j->first).at(time_bin_it) * single_pmt_norm + single_pmt_offset;
                double w2 = w * w;
                (*final_binned_yields_)[j->first].at(time_bin_it) += w;
                (*final_binned_squared_yields_)[j->first].at(time_bin_it) += w2; 
                
                // and for the derivs
                // note -- dwdtheta = derivs we binned * scaling
                // but dw2dtheta = 2 * w * dwdtheta
                double dw_dRs = this_event_binned_yields_deriv_Rs->at(j->first).at(time_bin_it) * single_pmt_norm;
                (*final_binned_yields_deriv_Rs_)[j->first].at(time_bin_it) += dw_dRs; 
                (*final_binned_squared_yields_deriv_Rs_)[j->first].at(time_bin_it) += 2.0 * w * dw_dRs;
                
                double dw_dRt = this_event_binned_yields_deriv_Rt->at(j->first).at(time_bin_it) * single_pmt_norm;
                (*final_binned_yields_deriv_Rt_)[j->first].at(time_bin_it) += dw_dRt; 
                (*final_binned_squared_yields_deriv_Rt_)[j->first].at(time_bin_it) += 2.0 * w * dw_dRt;
                
                double dw_dtau_s = this_event_binned_yields_deriv_tau_s->at(j->first).at(time_bin_it) * single_pmt_norm;
                (*final_binned_yields_deriv_tau_s_)[j->first].at(time_bin_it) += dw_dtau_s; 
                (*final_binned_squared_yields_deriv_tau_s_)[j->first].at(time_bin_it) += 2.0 * w * dw_dtau_s; 
                
                double dw_dtau_t = this_event_binned_yields_deriv_tau_t->at(j->first).at(time_bin_it) * single_pmt_norm;
                (*final_binned_yields_deriv_tau_t_)[j->first].at(time_bin_it) += dw_dtau_t;
                (*final_binned_squared_yields_deriv_tau_t_)[j->first].at(time_bin_it) += 2.0 * w * dw_dtau_t; 
                
                double dw_dtau_TPB = this_event_binned_yields_deriv_tau_TPB->at(j->first).at(time_bin_it) * single_pmt_norm;
                (*final_binned_yields_deriv_tau_TPB_)[j->first].at(time_bin_it) += dw_dtau_TPB; 
                (*final_binned_squared_yields_deriv_tau_TPB_)[j->first].at(time_bin_it) += 2.0 * w * dw_dtau_TPB; 
                
                double dw_dtime_offset = this_event_binned_yields_deriv_time_offset->at(j->first).at(time_bin_it) * single_pmt_norm;
                (*final_binned_yields_deriv_time_offset_)[j->first].at(time_bin_it) += dw_dtime_offset; 
                (*final_binned_squared_yields_deriv_time_offset_)[j->first].at(time_bin_it) += 2.0 * w * dw_dtime_offset; 
            }
        }
    }


    return std::make_tuple(final_binned_yields_, final_binned_squared_yields_,
                           final_binned_yields_deriv_Rs_, final_binned_squared_yields_deriv_Rs_,
                           final_binned_yields_deriv_Rt_, final_binned_squared_yields_deriv_Rt_,
                           final_binned_yields_deriv_tau_s_, final_binned_squared_yields_deriv_tau_s_,
                           final_binned_yields_deriv_tau_t_, final_binned_squared_yields_deriv_tau_t_,
                           final_binned_yields_deriv_tau_TPB_, final_binned_squared_yields_deriv_tau_TPB_,
                           final_binned_yields_deriv_time_offset_, final_binned_squared_yields_deriv_time_offset_);

}





