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
    sodium_events_constructor = new SodiumVertexDistribution();
    event_vertices = sodium_events_constructor->GetEventVertices(n_events_to_simulate, z_position);
    std::cout << "in GenerateExpectation::GetSodiumVertices" << std::endl;
}

void GenerateExpectation::GetYieldsAndOffsets(I3FramePtr geo_frame, double const & uv_absorption){
    yields_and_offset_constructor = new YieldsPerPMT();
    // now loop over events and get map between CCMPMTKey and PhotonYieldSummarySeries
    for (size_t sodium_it = 0; sodium_it < event_vertices->size(); sodium_it++){
        boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_per_event = yields_and_offset_constructor->GetAllYields(event_vertices->at(sodium_it), geo_frame, uv_absorption);
        yields_per_pmt_per_event.push_back(yields_per_event);
    }
    std::cout << "in GenerateExpectation::GetYieldsAndOffsets" << std::endl;
}

I3Vector<double> GenerateExpectation::LightProfile(double const & Rs, double const & Rt, double const & tau_s, double const & tau_t, double const & tau_rec, double const & tau_TPB,
                                       AnalyticLightYieldGenerator::LArLightProfileType const & light_profile_type){

    if (LAr_scintillation_light_constructor == nullptr){
        LAr_scintillation_light_constructor = new LArScintillationLightProfile();
    }

    I3Vector<double> light_profile;

    if (light_profile_type == AnalyticLightYieldGenerator::LArLightProfileType::Simplified){
        light_profile = LAr_scintillation_light_constructor->GetSimplifiedLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB);
    }
    else {
        light_profile = LAr_scintillation_light_constructor->GetFullLightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB);
    }
    std::cout << "in GenerateExpectation::LightProfile" << std::endl;

    return light_profile;
}

std::tuple<boost::shared_ptr<I3MapPMTKeyVectorDouble>, boost::shared_ptr<I3MapPMTKeyVectorDouble>> GenerateExpectation::GetExpectation(AnalyticLightYieldGenerator analytic_light_yield_setup,
                                                                                                                                       I3FramePtr geo_frame){
//boost::shared_ptr<I3MapPMTKeyVectorDouble> GenerateExpectation::GetExpectation(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame){

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
        GetYieldsAndOffsets(geo_frame, uv_absorption);
    }

    // now grab our light profile!
    I3Vector<double> LAr_light_profile = LightProfile(Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_profile_type);

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
                    light_val = LAr_light_profile.at(light_time_it) * relative_yields * normalization;

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
                final_binned_yields_->at(j->first).at(time_bin_it) += this_event_binned_yields->at(j->first).at(time_bin_it);
                final_binned_squared_yields_->at(j->first).at(time_bin_it) += std::pow(this_event_binned_yields->at(j->first).at(time_bin_it), 2.0);
            }
        }
    }


    return std::make_tuple(final_binned_yields_, final_binned_squared_yields_);
    //return final_binned_yields_;

}







