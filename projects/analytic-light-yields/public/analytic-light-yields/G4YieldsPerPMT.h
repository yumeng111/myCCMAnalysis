#ifndef G4YieldsPerPMT_H
#define G4YieldsPerPMT_H

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

#include "dataio/I3File.h"
#include "icetray/ctpl.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "simclasses/CCMMCPE.h"
#include "simclasses/PhotonSummary.h"

struct ccmmcpe_lightweight {
    float g4_time;
    float g4_distance_uv;
};


template <typename T>
struct G4PhotonPropagationJob {
    std::atomic<bool> running = false;
    size_t event_idx = 0;
    std::pair<size_t, size_t> event_ranges;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q11_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q11_binned_yields_squared = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q12_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q12_binned_yields_squared = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q21_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q21_binned_yields_squared = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q22_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q22_binned_yields_squared = nullptr;
};

template <typename T>
struct G4PhotonPropagationResult {
    size_t event_idx = 0;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q11_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q11_binned_yields_squared = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q12_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q12_binned_yields_squared = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q21_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q21_binned_yields_squared = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q22_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> Q22_binned_yields_squared = nullptr;
    bool done = false;
};

class G4YieldsPerPMT {
    ctpl::thread_pool pool;
    size_t max_cached_vertices = (size_t) 3000;
    std::exception_ptr teptr = nullptr;

public:
    G4YieldsPerPMT();
    template<typename T> void GetAllYields(size_t n_threads, std::vector<CCMPMTKey> const & keys_to_fit, double max_time, double light_times_per_bin,
                                           std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q11_events,
                                           std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q12_events,
                                           std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q21_events,
                                           std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q22_events,
                                           float z_below,
                                           float rayl_below,
                                           float z_above,
                                           float rayl_above,
                                           std::vector<T> UV_absorption_length,
                                           T desired_z_offset,
                                           T desired_rayl,
                                           bool fit_z_rayl,
                                           std::map<CCMPMTKey, std::vector<T>> & binned_yields,
                                           std::map<CCMPMTKey, std::vector<T>> & binned_yields_squared);

    std::vector<double> GetPlottingInformation(CCMPMTKey key, double uv_absorption, double scaling);
    template<typename T> std::tuple<std::map<CCMPMTKey, T>, std::map<CCMPMTKey, T>> GetSummedYieldsMap(std::vector<CCMPMTKey> keys_to_fit, T uv_abs, T scaling, double z_offset);
    void TimeComparison(std::vector<CCMPMTKey> keys_to_fit, double z_offset);

    std::map<double, std::deque<I3FramePtr>> G4Events_;
    bool grabbed_g4_events = false;
    std::vector<double> G4Times_;
    std::vector<double> CalculatedTimes_;

    std::vector<double> GetG4Times() { return G4Times_; };
    std::vector<double> GetCalculatedTimes() { return CalculatedTimes_; }
};

template<typename T> void ScaleAndBinG4Yields(std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & all_sodium_events,
                                              std::vector<CCMPMTKey> const & keys_to_fit,
                                              double & max_time,
                                              double & light_times_per_bin,
                                              size_t & event_start,
                                              size_t & event_end,
                                              std::vector<T> UV_absorption_length,
                                              std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields,
                                              std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields_squared){

    size_t n_data_bins = (max_time + 1) * 2.0;
    size_t n_bins = n_data_bins * light_times_per_bin;
    double bin_width = 2.0 / light_times_per_bin;
    // now let's loop over each event in our vector sodium_events
    for (size_t event_it = event_start; event_it < event_end; event_it ++){

        std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>> const & parsed_event = all_sodium_events.at(event_it);

        // now loop over all pmts we are fitting
        for(size_t k = 0; k < keys_to_fit.size(); k++){
            CCMPMTKey this_key = keys_to_fit.at(k);
 
            // let's check if this pmt is in our map to save binned events to, if not we add it!
            if (this_event_binned_yields->find(this_key) == this_event_binned_yields->end()) {
                (*this_event_binned_yields)[this_key] = std::vector<T>(n_bins, 0.0);
            }
            if (this_event_binned_yields_squared->find(this_key) == this_event_binned_yields_squared->end()) {
                (*this_event_binned_yields_squared)[this_key] = std::vector<T>(n_bins, 0.0);
            }

            // grab the list of lightweight CCMMCPEs for this key
            std::vector<ccmmcpe_lightweight> const & this_key_ccmmcpe = parsed_event.at(this_key);

            // loop over each lightweight CCMMCPE
            for (size_t m = 0; m < this_key_ccmmcpe.size(); m++){
                ccmmcpe_lightweight const & this_ccmmcpe = this_key_ccmmcpe.at(m);

                double distance_travelled_uv = this_ccmmcpe.g4_distance_uv / I3Units::cm;
                // now figure out scaling due to uv absorption
                T uv_scaling = exp(- distance_travelled_uv / UV_absorption_length.at(0));
                if (UV_absorption_length.at(1) > 0.0){
                    uv_scaling += exp(- distance_travelled_uv / UV_absorption_length.at(1));
                }
                double time = this_ccmmcpe.g4_time / I3Units::nanosecond;

                // now grab time bin and save!
                size_t bin_idx = time / bin_width;
                if(bin_idx >= n_bins) {
                    continue;
                }
                this_event_binned_yields->at(this_key).at(bin_idx) += uv_scaling;
                this_event_binned_yields_squared->at(this_key).at(bin_idx) += (uv_scaling * uv_scaling);
            }
        }
    }
}

template<typename T>
void FrameThread(std::atomic<bool> & running,
                 std::vector<CCMPMTKey> const & keys_to_fit,
                 std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q11_events,
                 std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q12_events,
                 std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q21_events,
                 std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q22_events,
                 size_t & event_start,
                 size_t & event_end,
                 std::vector<T> UV_absorption_length,
                 double max_time,
                 double light_times_per_bin,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q11_binned_yields,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q11_binned_yields_squared,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q12_binned_yields,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q12_binned_yields_squared,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q21_binned_yields,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q21_binned_yields_squared,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q22_binned_yields,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & Q22_binned_yields_squared,
                 bool fit_z_rayl){

    // just need to bin!
    ScaleAndBinG4Yields<T>(Q11_events, keys_to_fit, max_time, light_times_per_bin, event_start, event_end, UV_absorption_length, Q11_binned_yields, Q11_binned_yields_squared);

    if (fit_z_rayl){
        ScaleAndBinG4Yields<T>(Q12_events, keys_to_fit, max_time, light_times_per_bin, event_start, event_end, UV_absorption_length, Q12_binned_yields, Q12_binned_yields_squared);
        ScaleAndBinG4Yields<T>(Q21_events, keys_to_fit, max_time, light_times_per_bin, event_start, event_end, UV_absorption_length, Q21_binned_yields, Q21_binned_yields_squared);
        ScaleAndBinG4Yields<T>(Q22_events, keys_to_fit, max_time, light_times_per_bin, event_start, event_end, UV_absorption_length, Q22_binned_yields, Q22_binned_yields_squared);
    }


    running.store(false);
}

template<typename T>
void RunFrameThread(ctpl::thread_pool & pool,
                    G4PhotonPropagationJob<T> * job,
                    std::vector<CCMPMTKey> const & keys_to_fit,
                    std::vector<T> UV_absorption_length,
                    double max_time,
                    double light_times_per_bin,
                    std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q11_events,
                    std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q12_events,
                    std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q21_events,
                    std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q22_events,
                    bool fit_z_rayl){

    job->running.store(true);
    pool.push([ &running = job->running,
                &keys_to_fit,
                & Q11_events,
                & Q12_events,
                & Q21_events,
                & Q22_events,
                &event_start = job->event_ranges.first,
                &event_end = job->event_ranges.second,
                UV_absorption_length,
                max_time,
                light_times_per_bin,
                job,
                fit_z_rayl
    ] (int id) {
    FrameThread(running,
                keys_to_fit,
                Q11_events,
                Q12_events,
                Q21_events,
                Q22_events,
                event_start,
                event_end,
                UV_absorption_length,
                max_time,
                light_times_per_bin,
                job->Q11_binned_yields,
                job->Q11_binned_yields_squared,
                job->Q12_binned_yields,
                job->Q12_binned_yields_squared,
                job->Q21_binned_yields,
                job->Q21_binned_yields_squared,
                job->Q22_binned_yields,
                job->Q22_binned_yields_squared,
                fit_z_rayl);
    });

}

template<typename T> void GeneralInterpolation(T desired_x, float x_below, float x_above,
                                               std::vector<T> & Q_below, std::vector<T> & Q_above,
                                               std::vector<T> interpolation_results) {

    T coeff_below = (x_above - desired_x) / (x_above - x_below);
    T coeff_above = (desired_x - x_below) / (x_above - x_below);
    std::cout << "interpolating at desired x = " << desired_x << std::endl;
    std::cout << "x below = " << x_below << " and x above = " << x_above << std::endl;
    // now loop over each time
    for (size_t i = 0; i < Q_below.size(); i++){
        T this_Q_below = Q_below.at(i);
        T this_Q_above = Q_above.at(i);
        std::cout << "Q below = " << this_Q_below << " and Q above = " << this_Q_above << std::endl;
        T interp_result = (coeff_below * this_Q_below) + (coeff_above * this_Q_above);
        std::cout << "interp result at " << i << " = " << interp_result << std::endl;
        interpolation_results.at(i) += interp_result; 
    }
}

template<typename T> void G4YieldsPerPMT::GetAllYields(size_t n_threads, std::vector<CCMPMTKey> const & keys_to_fit, double max_time, double light_times_per_bin,
                                                       std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q11_events,
                                                       std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q12_events,
                                                       std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q21_events,
                                                       std::vector<std::map<CCMPMTKey, std::vector<ccmmcpe_lightweight>>> const & Q22_events,
                                                       float z_below,
                                                       float rayl_below,
                                                       float z_above,
                                                       float rayl_above,
                                                       std::vector<T> UV_absorption_length,
                                                       T desired_z_offset,
                                                       T desired_rayl,
                                                       bool fit_z_rayl,
                                                       std::map<CCMPMTKey, std::vector<T>> & binned_yields,
                                                       std::map<CCMPMTKey, std::vector<T>> & binned_yields_squared){

    // set up maps to save binned yields from each quadrant
    std::map<CCMPMTKey, std::vector<T>> Q11_binned_yields;
    std::map<CCMPMTKey, std::vector<T>> Q11_binned_yields_squared;
    std::map<CCMPMTKey, std::vector<T>> Q12_binned_yields;
    std::map<CCMPMTKey, std::vector<T>> Q12_binned_yields_squared;
    std::map<CCMPMTKey, std::vector<T>> Q21_binned_yields;
    std::map<CCMPMTKey, std::vector<T>> Q21_binned_yields_squared;
    std::map<CCMPMTKey, std::vector<T>> Q22_binned_yields;
    std::map<CCMPMTKey, std::vector<T>> Q22_binned_yields_squared;

    // make way to keep track of if our final results are initialized or not per pmt
    std::map<CCMPMTKey, bool> final_results_uninitialized;
    for (size_t k = 0; k < keys_to_fit.size(); k++){
        final_results_uninitialized[keys_to_fit.at(k)] = true;
    }

    // will be used for multi-threading our simulation jobs
    std::deque<G4PhotonPropagationJob<T> *> free_jobs;
    std::deque<G4PhotonPropagationJob<T> *> running_jobs;
    std::deque<G4PhotonPropagationResult<T>> results;
    size_t min_vertex_idx = 0;

    // set up our num threads
    size_t num_threads;
    if (n_threads == 0){
        num_threads = std::thread::hardware_concurrency();
    } else{
        num_threads = n_threads;
    }
    pool.resize(num_threads);

    // now let's chunk up our event to give each thread a range of events

    size_t n_events_per_thread = Q11_events.size() / num_threads;
    size_t left_over_events = Q11_events.size() - (num_threads * n_events_per_thread);

    std::vector<std::pair<size_t, size_t>> thread_event_ranges;

    size_t event_range_start = 0;
    size_t event_range_end = event_range_start;
    size_t accounted_for_left_over_events = 0;
    for (size_t thread_it = 0; thread_it < num_threads; thread_it++){

        event_range_start = event_range_end;
        event_range_end = event_range_start + n_events_per_thread; // each thread gets mandatory n_events_per_thread

        // check if we've dolled out all of the left over events yet
        if (accounted_for_left_over_events < left_over_events){
            event_range_end += 1;
            accounted_for_left_over_events += 1;
        }

        // save idxs
        thread_event_ranges.push_back(std::make_pair(event_range_start, event_range_end));
    }

    // now let's loop over our pre-made lists of sodium events for each thread
    for (size_t thread_it = 0; thread_it < thread_event_ranges.size(); ++thread_it) {

        while(true) {
            // Check if any jobs have finished
            for(int i=int(running_jobs.size())-1; i>=0; --i) {
			    if (teptr) {
				    try{
					    std::rethrow_exception(teptr);
				    }
				    catch(const std::exception &ex)
				    {
					    std::cerr << "Thread exited with exception: " << ex.what() << "\n";
				    }
			    }
			    if(not running_jobs.at(i)->running.load()) {
				    G4PhotonPropagationJob<T> * job = running_jobs.at(i);
                    running_jobs.erase(running_jobs.begin() + i);
                    free_jobs.push_back(job);
                    results.at(job->event_idx - min_vertex_idx).done = true;
                } else {
                    G4PhotonPropagationJob<T> * job = running_jobs.at(i);
                }
            }

            // Check for any done results and push the corresponding frames
            size_t results_done = 0;
            for(size_t i=0; i<results.size(); ++i) {
                if(results.at(i).done) {
                    // let's save our results binned yields to appropriate quadrant
                    for (auto e = results.at(i).Q11_binned_yields->begin(); e != results.at(i).Q11_binned_yields->end(); ++e) {
                        if (final_results_uninitialized[e->first]){
                            Q11_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q11_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q12_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q12_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q21_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q21_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q22_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                            Q22_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                            final_results_uninitialized[e->first] = false;
                        }
                        for (size_t b = 0; b < e->second.size(); b++){
                            Q11_binned_yields[e->first].at(b) += e->second.at(b);
                            Q11_binned_yields_squared[e->first].at(b) += (results.at(i).Q11_binned_yields_squared)->at(e->first).at(b);
                            if (fit_z_rayl){
                                Q12_binned_yields[e->first].at(b) += e->second.at(b);
                                Q12_binned_yields_squared[e->first].at(b) += (results.at(i).Q12_binned_yields_squared)->at(e->first).at(b);
                                Q21_binned_yields[e->first].at(b) += e->second.at(b);
                                Q21_binned_yields_squared[e->first].at(b) += (results.at(i).Q21_binned_yields_squared)->at(e->first).at(b);
                                Q22_binned_yields[e->first].at(b) += e->second.at(b);
                                Q22_binned_yields_squared[e->first].at(b) += (results.at(i).Q22_binned_yields_squared)->at(e->first).at(b);
                            }
                        }
                    }
                    // now reset our results object
                    results.at(i).Q11_binned_yields = nullptr;
                    results.at(i).Q11_binned_yields_squared = nullptr;
                    results.at(i).Q12_binned_yields = nullptr;
                    results.at(i).Q12_binned_yields_squared = nullptr;
                    results.at(i).Q21_binned_yields = nullptr;
                    results.at(i).Q21_binned_yields_squared = nullptr;
                    results.at(i).Q22_binned_yields = nullptr;
                    results.at(i).Q22_binned_yields_squared = nullptr;
                    results.at(i).event_idx = 0;
                    results_done += 1;
                } else {
                    break;
                }
            }
            if(results_done > 0) {
                results.erase(results.begin(), results.begin() + results_done);
                min_vertex_idx += results_done;
            }

            // Attempt to queue up a new job for the frame
            G4PhotonPropagationJob<T> * job = nullptr;

            if(free_jobs.size() > 0) {
                job = free_jobs.front();
                job->running.store(false);
                free_jobs.pop_front();
            } else if(running_jobs.size() < num_threads) {
                job = new G4PhotonPropagationJob<T>();
                job->running.store(false);
            }

            if(job != nullptr and results.size() < max_cached_vertices) {
                job->running.store(true);
                running_jobs.push_back(job);
                job->event_idx = thread_it;
                // add event ranges for this job
                job->event_ranges = thread_event_ranges.at(thread_it);
                // add emtpy maps to store binned yields for each quadrant
                job->Q11_binned_yields = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q11_binned_yields_squared = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q12_binned_yields = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q12_binned_yields_squared = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q21_binned_yields = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q21_binned_yields_squared = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q22_binned_yields = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                job->Q22_binned_yields_squared = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
                results.emplace_back();
                results.back().event_idx = job->event_idx;
                results.back().Q11_binned_yields = job->Q11_binned_yields;
                results.back().Q11_binned_yields_squared = job->Q11_binned_yields_squared;
                results.back().Q12_binned_yields = job->Q12_binned_yields;
                results.back().Q12_binned_yields_squared = job->Q12_binned_yields_squared;
                results.back().Q21_binned_yields = job->Q21_binned_yields;
                results.back().Q21_binned_yields_squared = job->Q21_binned_yields_squared;
                results.back().Q22_binned_yields = job->Q22_binned_yields;
                results.back().Q22_binned_yields_squared = job->Q22_binned_yields_squared;
                results.back().done = false;
                RunFrameThread<T>(pool, job, keys_to_fit, UV_absorption_length, max_time, light_times_per_bin,
                                  Q11_events, Q12_events, Q21_events, Q22_events, fit_z_rayl);
                break;
            } else if(job != nullptr) {
                free_jobs.push_back(job);
            }
        }
    }

    // final check for any running jobs
    while(running_jobs.size() > 0) {
        // Check if any jobs have finished
        for(int i=int(running_jobs.size())-1; i>=0; --i) {
            if(not running_jobs.at(i)->running.load()) {
                G4PhotonPropagationJob<T> * job = running_jobs.at(i);
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                results.at(job->event_idx - min_vertex_idx).done = true;
            }
        }
        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results.at(i).done) {
                for (auto e = results.at(i).Q11_binned_yields->begin(); e != results.at(i).Q11_binned_yields->end(); ++e) {
                    if (final_results_uninitialized[e->first]){
                        Q11_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q11_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q12_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q12_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q21_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q21_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q22_binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                        Q22_binned_yields_squared[e->first] = std::vector<T>(e->second.size(), 0.0);
                        final_results_uninitialized[e->first] = false;
                    }
                    for (size_t b = 0; b < e->second.size(); b++){
                        Q11_binned_yields[e->first].at(b) += e->second.at(b);
                        Q11_binned_yields_squared[e->first].at(b) += (results.at(i).Q11_binned_yields_squared)->at(e->first).at(b);
                            if (fit_z_rayl){
                                Q12_binned_yields[e->first].at(b) += e->second.at(b);
                                Q12_binned_yields_squared[e->first].at(b) += (results.at(i).Q12_binned_yields_squared)->at(e->first).at(b);
                                Q21_binned_yields[e->first].at(b) += e->second.at(b);
                                Q21_binned_yields_squared[e->first].at(b) += (results.at(i).Q21_binned_yields_squared)->at(e->first).at(b);
                                Q22_binned_yields[e->first].at(b) += e->second.at(b);
                                Q22_binned_yields_squared[e->first].at(b) += (results.at(i).Q22_binned_yields_squared)->at(e->first).at(b);
                            }
                    }
                }
                // now reset our results object
                results.at(i).Q11_binned_yields = nullptr;
                results.at(i).Q11_binned_yields_squared = nullptr;
                results.at(i).Q12_binned_yields = nullptr;
                results.at(i).Q12_binned_yields_squared = nullptr;
                results.at(i).Q21_binned_yields = nullptr;
                results.at(i).Q21_binned_yields_squared = nullptr;
                results.at(i).Q22_binned_yields = nullptr;
                results.at(i).Q22_binned_yields_squared = nullptr;
                results.at(i).event_idx = 0;
                results_done += 1;
            } else {
                break;
            }
        }
        if(results_done > 0) {
            results.erase(results.begin(), results.begin() + results_done);
            min_vertex_idx += results_done;
        }
    }
    // we also need to delete free_jobs now that we're done threading
    if (free_jobs.size() > 0){
        for(G4PhotonPropagationJob<T> * obj : free_jobs){
            delete obj;
        }
    }
    // ok we finished binning!!! time to interpolate!!!!
    // iterate through all the keys in our binned yields

    for (auto it = Q11_binned_yields.begin(); it != Q11_binned_yields.end(); ++it) {
        std::vector<T> & this_key_Q11 = it->second;
        std::vector<T> & this_key_Q11_squared = Q11_binned_yields_squared.at(it->first);

        if (fit_z_rayl){
            std::vector<T> & this_key_Q12 = Q12_binned_yields.at(it->first);
            std::vector<T> & this_key_Q21 = Q21_binned_yields.at(it->first);
            std::vector<T> & this_key_Q22 = Q22_binned_yields.at(it->first);
            std::vector<T> & this_key_Q12_squared = Q12_binned_yields_squared.at(it->first);
            std::vector<T> & this_key_Q21_squared = Q21_binned_yields_squared.at(it->first);
            std::vector<T> & this_key_Q22_squared = Q22_binned_yields_squared.at(it->first);


            // now interpolate in the x direction
            std::vector<T> interpolate_x_y_below (this_key_Q11.size(), 0.0);
            std::vector<T> interpolate_x_y_above (this_key_Q11.size(), 0.0);
            std::vector<T> interpolate_x_y_below_squared (this_key_Q11.size(), 0.0);
            std::vector<T> interpolate_x_y_above_squared (this_key_Q11.size(), 0.0);
            GeneralInterpolation<T>(desired_z_offset, z_below, z_above, this_key_Q11, this_key_Q21, interpolate_x_y_below);
            GeneralInterpolation<T>(desired_z_offset, z_below, z_above, this_key_Q12, this_key_Q22, interpolate_x_y_above);
            GeneralInterpolation<T>(desired_z_offset, z_below, z_above, this_key_Q11_squared, this_key_Q21_squared, interpolate_x_y_below_squared);
            GeneralInterpolation<T>(desired_z_offset, z_below, z_above, this_key_Q12_squared, this_key_Q22_squared, interpolate_x_y_above_squared);

            // now interpolate in the y direction
            std::vector<T> final_interpolation(this_key_Q11.size(), 0.0);
            std::vector<T> final_interpolation_squared(this_key_Q11_squared.size(), 0.0);
            GeneralInterpolation<T>(desired_rayl, rayl_below, rayl_above, interpolate_x_y_below, interpolate_x_y_above, final_interpolation);
            GeneralInterpolation<T>(desired_rayl, rayl_below, rayl_above, interpolate_x_y_below_squared, interpolate_x_y_above_squared, final_interpolation_squared);

            // all done! now save
            binned_yields[it->first] = final_interpolation;
            binned_yields_squared[it->first] = final_interpolation_squared;
        } else {
            binned_yields[it->first] = this_key_Q11;
            binned_yields_squared[it->first] = this_key_Q11_squared;
        }
    }

}

template<typename T> std::tuple<std::map<CCMPMTKey, T>, std::map<CCMPMTKey, T>> G4YieldsPerPMT::GetSummedYieldsMap(std::vector<CCMPMTKey> keys_to_fit, T uv_abs, T scaling, double z_offset){

    //// make sure we have filled G4Events_
    //if (grabbed_g4_events == false){
    //    std::vector<std::string> g4_fnames = {"/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumCenterHEEvents.i3.zst", "/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumPlus50HEEvents.i3.zst"};
    //    std::vector<double> z_locs = {0.0, 50.0};

    //    for (size_t f = 0; f < g4_fnames.size(); f++){
    //        std::deque<I3FramePtr> this_file_events;

    //        dataio::I3File g4_file(g4_fnames.at(f), dataio::I3File::Mode::read);
    //        while (g4_file.more()){
    //            I3FramePtr g4_frame = g4_file.pop_frame();
    //            this_file_events.push_back(g4_frame);
    //        }

    //        G4Events_[z_locs.at(f)] = this_file_events;
    //    }

    //}

    //// now get all yields
    //size_t n_threads = 0;
    //double max_time = 50.0;
    //std::map<CCMPMTKey, std::vector<T>> binned_yields;
    //std::map<CCMPMTKey, std::vector<T>> binned_yields_squared;

    //GetAllYields<T>(n_threads, keys_to_fit, G4Events_[z_offset], G4Events_[z_offset], false, max_time, 50.0, 50.0, {uv_abs, 0.0}, z_offset, scaling, binned_yields, binned_yields_squared);

    //// now sum binned_yields
    //std::map<CCMPMTKey, T> summed_yields;
    //std::map<CCMPMTKey, T> summed_yields_squared;

    //for (const auto& pair : binned_yields) {
    //    T total = 0;
    //    for (size_t i = 0; i < pair.second.size(); i++){
    //        total += pair.second.at(i);
    //    }
    //    summed_yields.insert(std::make_pair(pair.first, total));
    //}

    //for (const auto& pair : binned_yields_squared) {
    //    T total = 0;
    //    for (size_t i = 0; i < pair.second.size(); i++){
    //        total += pair.second.at(i);
    //    }
    //    summed_yields_squared.insert(std::make_pair(pair.first, total));
    //}

    //return std::make_tuple(summed_yields, summed_yields_squared);
    std::map<CCMPMTKey, T> summed_yields;
    std::map<CCMPMTKey, T> summed_yields_squared;
    return std::make_tuple(summed_yields, summed_yields_squared);
}

#endif
