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

template <typename T>
struct g4_photon_yield_summary {
    double time; // photon hit time
    T yield; // number of photon hits
};


template <typename T>
struct G4PhotonPropagationJob {
    std::atomic<bool> running = false;
    size_t event_idx = 0;
    std::vector<I3FramePtr> this_thread_events_above;
    std::vector<I3FramePtr> this_thread_events_below;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_squared = nullptr;
};

template <typename T>
struct G4PhotonPropagationResult {
    size_t event_idx = 0;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_squared = nullptr;
    bool done = false;
};

class G4YieldsPerPMT {
    ctpl::thread_pool pool;
    size_t max_cached_vertices = (size_t) 3000;
    std::exception_ptr teptr = nullptr;
    
public:
    G4YieldsPerPMT();
    template<typename T> void GetAllYields(size_t n_threads, std::vector<CCMPMTKey> const & keys_to_fit, std::deque<I3FramePtr> G4Events_above, std::deque<I3FramePtr> G4Events_below,
                                           bool need_to_interpolate, double max_time, double z_above, double z_below, std::vector<T> UV_absorption_length, T scaling, T z_offset, std::map<CCMPMTKey,
                                           std::vector<T>> & binned_yields, std::map<CCMPMTKey, std::vector<T>> & binned_square_yields);
    
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


template<typename T> void GrabG4Yields(I3FramePtr frame,
                                       std::vector<T> UV_absorption_length,
                                       std::shared_ptr<std::map<CCMPMTKey, std::vector<g4_photon_yield_summary<T>>>> & all_pmt_yields_map,
                                       std::vector<CCMPMTKey> const & keys_to_fit) {
    // some constants
    double c = 2.998 * std::pow(10.0, 8.0); // speed of light in m/s
    double c_cm_per_nsec = c * std::pow(10.0, -7.0); // speed of light in cm/nsec
    double uv_index_of_refraction = 1.358;
    double vis_index_of_refraction = 1.23;
    
    // lets grab necessary things from each frame
    boost::shared_ptr<CCMMCPESeriesMap const> CCMMCPEMap = frame->Get<boost::shared_ptr<CCMMCPESeriesMap const>>("PMTMCHitsMap");

    // now loop over all pmts we are fitting
    for(size_t k = 0; k < keys_to_fit.size(); k++){
        // let's check if this pmt is in our map to save all photon yield summary info, if not we add it!
        if (all_pmt_yields_map->find(keys_to_fit.at(k)) == all_pmt_yields_map->end()) {
            (*all_pmt_yields_map)[keys_to_fit.at(k)] = std::vector<g4_photon_yield_summary<T>>{};
        }
        
        // grab the list of CCMMCPEs for this key
        CCMMCPESeries this_key_ccmmcpe = CCMMCPEMap->at(keys_to_fit.at(k));
        
        // loop over each CCMMCPE 
        for (size_t m = 0; m < this_key_ccmmcpe.size(); m++){
            CCMMCPE this_ccmmcpe = this_key_ccmmcpe.at(m);
            // first check the wavelength!! only continue if visible
            double wavelength = this_ccmmcpe.wavelength * 1e9; // units of nm
            if (wavelength > 325.0){
                // check to make sure it had photon summary information
                if ((this_ccmmcpe.g4_distance_uv + this_ccmmcpe.g4_distance_visible) > 0.0){
                    
                    double distance_travelled_uv = this_ccmmcpe.g4_distance_uv / I3Units::cm;
                    // now figure out scaling due to uv absorption
                    T uv_scaling = exp(- distance_travelled_uv / UV_absorption_length.at(0));
                    if (UV_absorption_length.at(1) > 0.0){
                        uv_scaling += exp(- distance_travelled_uv / UV_absorption_length.at(1));
                    }
                    double time = this_ccmmcpe.g4_time; // already in units of nsec 
                    
                    // now make a new yield to save
                    g4_photon_yield_summary<T> this_yield_summary;
                    this_yield_summary.time = time; 
                    this_yield_summary.yield = uv_scaling; 
                    all_pmt_yields_map->at(keys_to_fit.at(k)).push_back(this_yield_summary);
                }
            }
        }
    }
}

template<typename T> void BinEvents(std::vector<I3FramePtr> sodium_events, std::vector<T> UV_absorption_length, double max_time, std::vector<CCMPMTKey> const & keys_to_fit,
                                    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields, std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields_squared){

    // now let's loop over each event in our vector sodium_events
    for (size_t event_it = 0; event_it < sodium_events.size(); event_it ++){

        // call simulation code
        std::shared_ptr<std::map<CCMPMTKey, std::vector<g4_photon_yield_summary<T>>>> this_event_pmt_yields_map = std::make_shared<std::map<CCMPMTKey, std::vector<g4_photon_yield_summary<T>>>>();
        GrabG4Yields(sodium_events.at(event_it), UV_absorption_length, this_event_pmt_yields_map, keys_to_fit);
        
        // now let's bin this_event_pmt_yields_map 
        size_t n_bins = max_time / 2.0;
        
        // loop over each key are we calculating yields + offset for
        for (size_t k = 0; k < keys_to_fit.size(); k++){ 
            CCMPMTKey key = keys_to_fit.at(k);
            (*this_event_binned_yields)[key] = std::vector<T>(n_bins, 0.0);
            (*this_event_binned_yields_squared)[key] = std::vector<T>(n_bins, 0.0);
            
            auto i = this_event_pmt_yields_map->find(key);
            if (i == this_event_pmt_yields_map->end()) {
                continue;
            }
            std::vector<g4_photon_yield_summary<T>> const & yields = i->second;
            if(yields.size() == 0) {
                continue;
            }
            for(size_t yield_it = 0; yield_it < yields.size(); ++yield_it) {
                g4_photon_yield_summary<T> const & yield = yields.at(yield_it);
                size_t bin_idx = yield.time / 2.0;
                if(bin_idx >= n_bins) {
                    continue;
                }
                this_event_binned_yields->at(key).at(bin_idx) += yield.yield;
                this_event_binned_yields_squared->at(key).at(bin_idx) += (yield.yield * yield.yield);
            }
        }
    }

}

template<typename T>
void FrameThread(std::atomic<bool> & running,
                 std::vector<CCMPMTKey> const & keys_to_fit,
                 std::vector<I3FramePtr> sodium_events_above,
                 std::vector<I3FramePtr> sodium_events_below,
                 std::vector<T> UV_absorption_length,
                 T scaling,
                 T z_offset,
                 double max_time,
                 double z_above,
                 double z_below,
                 bool need_to_interpolate,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields_squared) {

    //if (need_to_interpolate){
    // make some maps to hold binned yields for above and below idx
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_above = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_squared_above = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
    
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_below = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_squared_below = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>();

    // now grab our g4 yields + bin both of these guys
    BinEvents<T>(sodium_events_above, UV_absorption_length, max_time, keys_to_fit, this_event_binned_yields_above, this_event_binned_yields_squared_above);
    BinEvents<T>(sodium_events_below, UV_absorption_length, max_time, keys_to_fit, this_event_binned_yields_below, this_event_binned_yields_squared_below);

    size_t n_bins = max_time / 2.0;
    
    // now go through each binned yields + linear interpolate
    for (auto it = this_event_binned_yields_above->begin(); it != this_event_binned_yields_above->end(); it++) {
        // add this key to this_event_binned_yields
        (*this_event_binned_yields)[it->first] = std::vector<T>(n_bins, 0.0);
        (*this_event_binned_yields_squared)[it->first] = std::vector<T>(n_bins, 0.0);

        // grab our binned yields for each key
        std::vector<T>& this_pmt_yields_above = it->second;
        std::vector<T>& this_pmt_yields_below = this_event_binned_yields_below->at(it->first);
        
        std::vector<T>& this_pmt_yields_squared_above = this_event_binned_yields_squared_above->at(it->first);
        std::vector<T>& this_pmt_yields_squared_below = this_event_binned_yields_squared_below->at(it->first);
        
        // iterate through binned yields + linear interpolate + scale
        for (size_t y = 0; y < this_pmt_yields_above.size(); y++){
            T yields_above = this_pmt_yields_above.at(y);
            T yields_below = this_pmt_yields_below.at(y);

            T interpolated_yields = yields_below + (z_offset - z_below) * ((yields_above - yields_below) / (z_above - z_below));    

            // now do the same for yields squared 
            T yields_squared_above = this_pmt_yields_squared_above.at(y);
            T yields_squared_below = this_pmt_yields_squared_below.at(y);

            T interpolated_yields_squared = yields_squared_below + (z_offset - z_below) * ((yields_squared_above - yields_squared_below) / (z_above - z_below));    
            
            // now save
            this_event_binned_yields->at(it->first).at(y) += interpolated_yields * scaling;
            this_event_binned_yields_squared->at(it->first).at(y) += interpolated_yields_squared * scaling * scaling;
        }
    }
    //}
    //else {
    //    BinEvents<T>(sodium_events_above, UV_absorption_length, max_time, keys_to_fit, this_event_binned_yields, this_event_binned_yields_squared);
    //    
    //    // apply scaling
    //    for (auto it = this_event_binned_yields->begin(); it != this_event_binned_yields->end(); it++) {

    //        // grab our binned yields for each key
    //        std::vector<T>& this_pmt_yields = it->second;
    //        
    //        // iterate through binned yields + scale
    //        for (size_t y = 0; y < this_pmt_yields.size(); y++){
    //            // now save
    //            this_event_binned_yields->at(it->first).at(y) *= scaling;
    //            this_event_binned_yields_squared->at(it->first).at(y) *= (scaling * scaling);
    //        }
    //    }
    //}
     

    running.store(false);
}

template<typename T>
void RunFrameThread(ctpl::thread_pool & pool,
                    G4PhotonPropagationJob<T> * job,
                    std::vector<CCMPMTKey> const & keys_to_fit,
                    std::vector<T> UV_absorption_length,
                    T scaling,
                    T z_offset,
                    double max_time,
                    double z_above,
                    double z_below,
                    bool need_to_interpolate){

    job->running.store(true);
    pool.push([ &running = job->running, &keys_to_fit, sodium_events_above = job->this_thread_events_above, sodium_events_below = job->this_thread_events_below,
                UV_absorption_length, scaling, z_offset, max_time, z_above, z_below, need_to_interpolate, job 
    ] (int id) {
    FrameThread(running,
                keys_to_fit,
                sodium_events_above,
                sodium_events_below,
                UV_absorption_length, 
                scaling,
                z_offset,
                max_time,
                z_above,
                z_below,
                need_to_interpolate,
                job->this_event_binned_yields,
                job->this_event_binned_yields_squared);
    });

}

template<typename T> void G4YieldsPerPMT::GetAllYields(size_t n_threads, std::vector<CCMPMTKey> const & keys_to_fit, std::deque<I3FramePtr> G4Events_above, std::deque<I3FramePtr> G4Events_below, 
                                                       bool need_to_interpolate, double max_time, double z_above, double z_below, std::vector<T> UV_absorption_length, T scaling, T z_offset, 
                                                       std::map<CCMPMTKey, std::vector<T>> & binned_yields, std::map<CCMPMTKey, std::vector<T>> & binned_square_yields){
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

    // now let's chunk up our event_vertices so we can give 1 vector of event_vertices to each thread
    std::vector<std::vector<I3FramePtr>> events_per_thread_above;
    std::vector<std::vector<I3FramePtr>> events_per_thread_below;
    
    size_t n_events_per_thread = G4Events_above.size() / num_threads;
    size_t left_over_events = G4Events_above.size() - (num_threads * n_events_per_thread);
    
    size_t sodium_event_idx = 0;
    size_t accounted_for_left_over_events = 0;
    for (size_t thread_it = 0; thread_it < num_threads; thread_it++){
        // make vector to hold events
        std::vector<I3FramePtr> this_thread_vector_of_events_above;
        std::vector<I3FramePtr> this_thread_vector_of_events_below;

        // things to keep track of
        size_t n_events_on_this_thread = 0;

        // now time to actaully add the sodium events to this_thread_vector_of_events
        while (n_events_on_this_thread < n_events_per_thread){
            this_thread_vector_of_events_above.push_back(G4Events_above.at(sodium_event_idx));
            this_thread_vector_of_events_below.push_back(G4Events_below.at(sodium_event_idx));
            n_events_on_this_thread += 1;
            sodium_event_idx += 1;
        }

        // check if we've dolled out all of the left over events yet
        if (accounted_for_left_over_events < left_over_events){
            this_thread_vector_of_events_above.push_back(G4Events_above.at(sodium_event_idx));
            this_thread_vector_of_events_below.push_back(G4Events_below.at(sodium_event_idx));
            sodium_event_idx += 1;
            accounted_for_left_over_events += 1;
        }

        events_per_thread_above.push_back(this_thread_vector_of_events_above);
        events_per_thread_below.push_back(this_thread_vector_of_events_below);
    }

    // now let's loop over our pre-made lists of sodium events for each thread
    for (size_t sodium_it = 0; sodium_it < events_per_thread_above.size(); ++sodium_it) {
        std::vector<I3FramePtr> this_thread_events_above = events_per_thread_above.at(sodium_it);
        std::vector<I3FramePtr> this_thread_events_below = events_per_thread_below.at(sodium_it);

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
                    // let's save this_event_binned_yields to our binned_yields map
                    for (auto e = results.at(i).this_event_binned_yields->begin(); e != results.at(i).this_event_binned_yields->end(); ++e) {
                        // now let's take our binned_yields at this pmt and add vector appropraitly
                        // and same for the binned_yields_squared
                        if (binned_yields[e->first].size() == 0){
                            binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                            binned_square_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                        }
                        for (size_t b = 0; b < e->second.size(); b++){
                            binned_yields[e->first].at(b) += e->second.at(b);
                            binned_square_yields[e->first].at(b) += (results.at(i).this_event_binned_yields_squared)->at(e->first).at(b);
                        }
                    }
                    // now reset our results object
                    results.at(i).this_event_binned_yields = nullptr;
                    results.at(i).this_event_binned_yields_squared = nullptr;
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
                job->event_idx = sodium_it;
                job->this_thread_events_above = this_thread_events_above; 
                job->this_thread_events_below = this_thread_events_below; 
                job->this_event_binned_yields = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>(); 
                job->this_event_binned_yields_squared = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>(); 
                results.emplace_back();
                results.back().event_idx = job->event_idx;
                results.back().this_event_binned_yields = job->this_event_binned_yields;
                results.back().this_event_binned_yields_squared = job->this_event_binned_yields_squared;
                results.back().done = false;
                RunFrameThread<T>(pool, job, keys_to_fit, UV_absorption_length, scaling, z_offset, max_time, z_above, z_below, need_to_interpolate);
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
                // let's save this_event_binned_yields to our binned_yields map
                for (auto e = results.at(i).this_event_binned_yields->begin(); e != results.at(i).this_event_binned_yields->end(); ++e) {
                    // now let's take our binned_yields at this pmt and add vector appropraitly
                    // and same for the binned_yields_squared
                    if (binned_yields[e->first].size() == 0){
                        binned_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                        binned_square_yields[e->first] = std::vector<T>(e->second.size(), 0.0);
                    }
                    for (size_t b = 0; b < e->second.size(); b++){
                        binned_yields[e->first].at(b) += e->second.at(b);
                        binned_square_yields[e->first].at(b) += results.at(i).this_event_binned_yields_squared->at(e->first).at(b);
                    }
                }
                // now reset our results object
                results.at(i).this_event_binned_yields = nullptr;
                results.at(i).this_event_binned_yields_squared = nullptr;
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

}

template<typename T> std::tuple<std::map<CCMPMTKey, T>, std::map<CCMPMTKey, T>> G4YieldsPerPMT::GetSummedYieldsMap(std::vector<CCMPMTKey> keys_to_fit, T uv_abs, T scaling, double z_offset){
    
    // make sure we have filled G4Events_
    if (grabbed_g4_events == false){
        std::vector<std::string> g4_fnames = {"/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumCenterHEEvents.i3.zst", "/Users/darcybrewuser/workspaces/CCM/notebooks/G4SodiumPlus50HEEvents.i3.zst"};
        std::vector<double> z_locs = {0.0, 50.0};

        for (size_t f = 0; f < g4_fnames.size(); f++){
            std::deque<I3FramePtr> this_file_events;

            dataio::I3File g4_file(g4_fnames.at(f), dataio::I3File::Mode::read);
            while (g4_file.more()){
                I3FramePtr g4_frame = g4_file.pop_frame();
                this_file_events.push_back(g4_frame);
            }

            G4Events_[z_locs.at(f)] = this_file_events;
        }

    }
    
    // now get all yields
    size_t n_threads = 0;
    double max_time = 50.0;
    std::map<CCMPMTKey, std::vector<T>> binned_yields;
    std::map<CCMPMTKey, std::vector<T>> binned_square_yields;

    GetAllYields<T>(n_threads, keys_to_fit, G4Events_[z_offset], G4Events_[z_offset], false, max_time, 50.0, 50.0, {uv_abs, 0.0}, z_offset, scaling, binned_yields, binned_square_yields);

    // now sum binned_yields
    std::map<CCMPMTKey, T> summed_yields;
    std::map<CCMPMTKey, T> summed_yields_squared;

    for (const auto& pair : binned_yields) {
        T total = 0;
        for (size_t i = 0; i < pair.second.size(); i++){
            total += pair.second.at(i);
        }
        summed_yields.insert(std::make_pair(pair.first, total));
    }
    
    for (const auto& pair : binned_square_yields) {
        T total = 0;
        for (size_t i = 0; i < pair.second.size(); i++){
            total += pair.second.at(i);
        }
        summed_yields_squared.insert(std::make_pair(pair.first, total));
    }

    return std::make_tuple(summed_yields, summed_yields_squared);
}

#endif
