#ifndef YieldsPerPMT_H
#define YieldsPerPMT_H

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

#include "icetray/ctpl.h"
#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/physics/HESodiumEvent.h"

template <typename T>
struct photon_yield_summary {
    enum class PhotonSource : int8_t {
        Unknown = 0,
        Vertex = 1,
        TPBFoil = 2,
    };
    T time; // photon hit time
    T yield; // number of photon hits
    PhotonSource photon_source;
};

struct yields_pmt_info {
    CCMPMTKey key;
    double pmt_x = 0.0;
    double pmt_y = 0.0;
    double pmt_z = 0.0;
    double facing_dir_x = 0.0;
    double facing_dir_y = 0.0;
    double facing_dir_z = 0.0;
    double coating_flag = 0.0;
    double pmt_facing_area = 0.0;
    double pmt_side_area = 0.0;
};

struct secondary_loc_info {
    double loc_x = 0.0;
    double loc_y = 0.0;
    double loc_z = 0.0;
    double facing_dir_x = 0.0;
    double facing_dir_y = 0.0;
    double facing_dir_z = 0.0;
    double pmt_portion = 0.0;
    double loc_facing_area = 0.0;
    double loc_side_area = 0.0;
};

template <typename T>
struct PhotonPropagationJob {
    std::atomic<bool> running = false;
    size_t event_idx = 0;
    std::vector<HESodiumEvent> this_event;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_squared = nullptr;
};

template <typename T>
struct PhotonPropagationResult {
    size_t event_idx = 0;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields = nullptr;
    std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> this_event_binned_yields_squared = nullptr;
    bool done = false;
};

class YieldsPerPMT {
    std::string geometry_name_ = std::string("CCMGeometry");
    I3FramePtr geo_frame;

    double portion_light_reflected_by_tpb_ = 1.0;
    double desired_chunk_width_ = 20.0; // use 5 for finer binning
    double desired_chunk_height_ = 20.0; // use 5 for finer binning
    std::vector<yields_pmt_info> pmt_parsed_information_;
    std::vector<secondary_loc_info> locations_to_check_information_;

    unsigned int coated_omtype = (unsigned int)10;
    unsigned int uncoated_omtype = (unsigned int)20;
    std::vector<std::vector<double>> uncoated_pmt_locs_top_;
    std::vector<std::vector<double>> coated_pmt_locs_top_;
    std::vector<std::vector<double>> uncoated_pmt_locs_bottom_;
    std::vector<std::vector<double>> coated_pmt_locs_bottom_;
    std::vector<std::vector<double>> uncoated_pmt_locs_side_;
    std::vector<std::vector<double>> coated_pmt_locs_side_;

    // defining some geomtry things for modelling the detector...not the most elegant
    double pmt_radius = 10.16; //radius in cm^2
    double pmt_facing_area = M_PI * std::pow(pmt_radius, 2.0); // units of cm^2
    double pmt_side_area_factor = 0.217549; // unitless
    double pmt_side_area = pmt_facing_area * pmt_side_area_factor;
    double cylinder_height = 48.4 * 2.54; // units of cm!
    double cylinder_radius = 40.0 * 2.54; // units of cm!
    double cylinder_circumference = M_PI * 2.0 * cylinder_radius;
    double cylinder_max_x = cylinder_radius;
    double cylinder_min_x = - cylinder_max_x;
    double cylinder_max_y = cylinder_radius;
    double cylinder_min_y = - cylinder_max_y;
    double cylinder_max_z = cylinder_height / 2.0;
    double cylinder_min_z = - cylinder_max_z;
    double chunk_side_area_factor = 0.1; // this number is a guess..maybe model it one day

    // some constants we use for the simulation
    double c = 2.998 * std::pow(10.0, 8.0); // speed of light in m/s
    double c_cm_per_nsec_ = c * std::pow(10.0, -7.0); // speed of light in cm/nsec
    double uv_index_of_refraction_ = 1.358;
    double vis_index_of_refraction_ = 1.23;

    bool set_secondary_locs_ = false;

    ctpl::thread_pool pool;
    size_t max_cached_vertices = (size_t) 2000;
    std::exception_ptr teptr = nullptr;

public:
    YieldsPerPMT() = default;
    YieldsPerPMT(I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height);

    void SetPMTInformation(I3FramePtr frame);
    void SetGeoFrame(I3FramePtr geo_frame);
    void SetChunks(double chunk_width, double chunk_height);

    template<typename T> void GetAllYields(size_t n_threads, boost::shared_ptr<HESodiumEventSeries> event_vertices, T UV_absorption_length, std::vector<CCMPMTKey> const & keys_to_fit,
                                                     double max_time, std::map<CCMPMTKey, std::vector<T>> & binned_yields, std::map<CCMPMTKey, std::vector<T>> & binned_square_yields);

};

inline void get_solid_angle_and_distance_vertex_to_location(double const & vertex_x,
                                                     double const & vertex_y,
                                                     double const & vertex_z,
                                                     double const & loc_x,
                                                     double const & loc_y,
                                                     double const & loc_z,
                                                     double const & facing_dir_x,
                                                     double const & facing_dir_y,
                                                     double const & facing_dir_z,
                                                     double const & facing_area,
                                                     double const & side_area,
                                                     double & omega,
                                                     double & D) {
    double v_to_loc_dir_x = loc_x - vertex_x;
    double v_to_loc_dir_y = loc_y - vertex_y;
    double v_to_loc_dir_z = loc_z - vertex_z;

    double v_to_loc_dir_norm = std::sqrt(std::pow(v_to_loc_dir_x, 2) + std::pow(v_to_loc_dir_y, 2) + std::pow(v_to_loc_dir_z, 2));
    v_to_loc_dir_x /= v_to_loc_dir_norm;
    v_to_loc_dir_y /= v_to_loc_dir_norm;
    v_to_loc_dir_z /= v_to_loc_dir_norm;

    double dot = std::abs((v_to_loc_dir_x * facing_dir_x) + (v_to_loc_dir_y * facing_dir_y) + (v_to_loc_dir_z * facing_dir_z));
    double effective_area = facing_area * dot + (1.0 - dot) * side_area;
    double effective_radius = std::sqrt(effective_area / M_PI);
    double effective_diameter = 2.0 * effective_radius;

    D = std::sqrt(std::pow(loc_x - vertex_x, 2) + std::pow(loc_y - vertex_y, 2) + std::pow(loc_z - vertex_z, 2));
    double d = effective_diameter;

    if (D == 0.0){
        omega = 2.0 * M_PI;
    }
    else {
        double delta = 2.0 * std::atan(d / (2.0 * D)); // opening angle
        omega = 2.0 * M_PI * (1 - std::cos(delta / 2.0));
    }
}

template<typename T> T get_light_yield(double const & omega,
                                       double const & distance_travelled,
                                       T absorption_length) {

    T light_yield = omega * exp(- distance_travelled / absorption_length);
    return light_yield;
}

// let's make a function to calculate how much light from a vertex goes into each PMT

template<typename T> void vertex_to_pmt_propagation(double const & c_cm_per_nsec,
                                                    double const & uv_index_of_refraction,
                                                    double const & vis_index_of_refraction,
                                                    T n_photons_produced,
                                                    T time_offset,
                                                    T absorption_length,
                                                    std::vector<yields_pmt_info> const & pmt_parsed_information_,
                                                    bool const & is_visible,
                                                    I3Position const & vertex,
                                                    std::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> & all_pmt_yields_map,
                                                    std::vector<CCMPMTKey> const & keys_to_fit) {

    // let's start by looping over our pmt parsed information
    for (size_t pmt_it = 0; pmt_it < pmt_parsed_information_.size(); ++pmt_it) {

        double const & pmt_x_loc = pmt_parsed_information_.at(pmt_it).pmt_x;
        double const & pmt_y_loc = pmt_parsed_information_.at(pmt_it).pmt_y;
        double const & pmt_z_loc = pmt_parsed_information_.at(pmt_it).pmt_z;
        double const & facing_dir_x = pmt_parsed_information_.at(pmt_it).facing_dir_x;
        double const & facing_dir_y = pmt_parsed_information_.at(pmt_it).facing_dir_y;
        double const & facing_dir_z = pmt_parsed_information_.at(pmt_it).facing_dir_z;
        double const & coating_flag = pmt_parsed_information_.at(pmt_it).coating_flag; // coating flag == 1.0 for coated pmts and 0.0 for uncoated!
        double const & pmt_facing_area = pmt_parsed_information_.at(pmt_it).pmt_facing_area;
        double const & pmt_side_area = pmt_parsed_information_.at(pmt_it).pmt_side_area;
        CCMPMTKey const & pmt_key = pmt_parsed_information_.at(pmt_it).key;

        // check if this is a pmt key we should fit
        bool fit = keys_to_fit.size() == 0 or std::find(keys_to_fit.begin(), keys_to_fit.end(), pmt_key) != keys_to_fit.end();
        if (!fit) {
            continue;
        }
        // ok, this must be a pmt we want to fit!
        if (is_visible == false and coating_flag == 0.0) {
            // this is the case where we are propagating UV light and we have an uncoated pmt...so we won't see anything
            continue;
        }

        // let's check if this pmt is in our map to save all photon yield summary info, if not we add it!
        if (all_pmt_yields_map->find(pmt_key) == all_pmt_yields_map->end()) {
            (*all_pmt_yields_map)[pmt_key] = std::vector<photon_yield_summary<T>>{};
        }

        else {
            double omega;
            double distance_travelled;
            // let's get solid angle and distance from vertex to this pmt
            get_solid_angle_and_distance_vertex_to_location(vertex.GetX() / I3Units::cm, vertex.GetY() / I3Units::cm, vertex.GetZ() / I3Units::cm,
                                                            pmt_x_loc, pmt_y_loc, pmt_z_loc,
                                                            facing_dir_x, facing_dir_y, facing_dir_z,
                                                            pmt_facing_area, pmt_side_area, omega, distance_travelled);

            // now call function to get light yields
            T photons_in_this_pmt = get_light_yield(omega, distance_travelled, absorption_length);
            photons_in_this_pmt *= n_photons_produced;

            T travel_time;
            // ok so we have the photons seen by this pmt, let's get the propagation times
            if (is_visible) {
                travel_time = distance_travelled / (c_cm_per_nsec / vis_index_of_refraction); // units of nsec
            }
            else{
                travel_time = distance_travelled / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec
            }

            // now all that's left to do is save!
            photon_yield_summary<T> this_pmt_yield_summary;
            this_pmt_yield_summary.time = travel_time + time_offset;
            this_pmt_yield_summary.yield = photons_in_this_pmt;
            if (is_visible){
                this_pmt_yield_summary.photon_source = photon_yield_summary<T>::PhotonSource::TPBFoil;
            }
            else {
                this_pmt_yield_summary.photon_source = photon_yield_summary<T>::PhotonSource::Vertex;
            }
            for (auto it = all_pmt_yields_map->begin(); it != all_pmt_yields_map->end(); ++it) {
                CCMPMTKey map_key = it->first;
                if (map_key == pmt_key){
                    all_pmt_yields_map->at(map_key).push_back(this_pmt_yield_summary);
                }
            }
        }
    }
}

// now let's do vertex --> TPB --> PMT propagation

template<typename T> void vertex_to_TPB_to_PMT_propagation(double const & c_cm_per_nsec,
                                                           double const & uv_index_of_refraction,
                                                           double const & vis_index_of_refraction,
                                                           T n_photons_produced,
                                                           T UV_absorption_length,
                                                           T vis_absorption_length,
                                                           I3Position const & vertex,
                                                           std::vector<yields_pmt_info> const & pmt_parsed_information_,
                                                           std::vector<secondary_loc_info> const & locations_to_check_information_,
                                                           std::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> & all_pmt_yields_map,
                                                           std::vector<CCMPMTKey> const & keys_to_fit) {
    // define some things
    double omega;
    double distance_travelled;
    double loc_x;
    double loc_y;
    double loc_z;
    double facing_dir_x;
    double facing_dir_y;
    double facing_dir_z;
    double pmt_portion;
    double facing_area;
    double side_area;

    // let's start by looping our our secondary locations parsed information
    for (size_t loc_it = 0; loc_it < locations_to_check_information_.size(); ++loc_it) {

        loc_x = locations_to_check_information_.at(loc_it).loc_x;
        loc_y = locations_to_check_information_.at(loc_it).loc_y;
        loc_z = locations_to_check_information_.at(loc_it).loc_z;
        facing_dir_x = locations_to_check_information_.at(loc_it).facing_dir_x;
        facing_dir_y = locations_to_check_information_.at(loc_it).facing_dir_y;
        facing_dir_z = locations_to_check_information_.at(loc_it).facing_dir_z;
        pmt_portion = locations_to_check_information_.at(loc_it).pmt_portion; // portion_light_reflected_by_tpb_ if we are NOT on a pmt, otherwise 0.5 * portion_light_reflected_by_tpb_
        facing_area = locations_to_check_information_.at(loc_it).loc_facing_area;
        side_area = locations_to_check_information_.at(loc_it).loc_side_area;

        // let's get solid angle and distance from vertex to this loc
        get_solid_angle_and_distance_vertex_to_location(vertex.GetX() / I3Units::cm, vertex.GetY() / I3Units::cm, vertex.GetZ() / I3Units::cm,
                                                        loc_x, loc_y, loc_z,
                                                        facing_dir_x, facing_dir_y, facing_dir_z,
                                                        facing_area, side_area, omega, distance_travelled);
        // now call function to get light yields
        // we are propagating scint light from vertex to TPB locs, so want to use UV absorption length
        T photons_at_secondary_location = get_light_yield(omega, distance_travelled, UV_absorption_length);
        photons_at_secondary_location *= n_photons_produced;
        photons_at_secondary_location *= pmt_portion;

        // now let's get the travel time as well
        T travel_time = distance_travelled / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec

        // ok so now we have time offset and yields at this loc, let's calculate what gets into our PMTs
        vertex_to_pmt_propagation<T>(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction,
                                  photons_at_secondary_location, travel_time, vis_absorption_length,
                                  pmt_parsed_information_, true, I3Position(loc_x * I3Units::cm, loc_y * I3Units::cm, loc_z * I3Units::cm),
                                  all_pmt_yields_map, keys_to_fit);

        // since vertex_to_pmt_propagation saves info -- we're done!
    }
}

template<typename T> void PutSimulationStepsTogether(HESodiumEvent const & soidum_event,
                                                     double const & c_cm_per_nsec,
                                                     double const & uv_index_of_refraction,
                                                     double const & vis_index_of_refraction,
                                                     T UV_absorption_length,
                                                     T vis_absorption_length,
                                                     std::vector<yields_pmt_info> const & pmt_parsed_information_,
                                                     std::vector<secondary_loc_info> const & locations_to_check_information_,
                                                     std::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> & all_pmt_yields_map,
                                                     std::vector<CCMPMTKey> const & keys_to_fit) {

    // ok so let's grab our vertex locations from our soidum_event
    I3Position photon_vertex = soidum_event.photon_vertex;
    I3Position electron_vertex = soidum_event.electron_vertex;
    I3Position positron_vertex = soidum_event.positron_vertex;

    T n_photons_1275 = 1.275;
    T n_photons_511 = 0.511;
    T default_time_offset = 0.0;

    // now let's call our functions to get direct vertex --> PMT yields (UV scint light)
    vertex_to_pmt_propagation<T>(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_1275,
                              default_time_offset, UV_absorption_length, pmt_parsed_information_,
                              false, photon_vertex, all_pmt_yields_map, keys_to_fit);
    vertex_to_pmt_propagation<T>(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_511,
                              default_time_offset, UV_absorption_length, pmt_parsed_information_,
                              false, electron_vertex, all_pmt_yields_map, keys_to_fit);
    vertex_to_pmt_propagation<T>(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_511,
                              default_time_offset, UV_absorption_length, pmt_parsed_information_,
                              false, positron_vertex, all_pmt_yields_map, keys_to_fit);
    
    // now let's call our functions to get indirect vertex --> TPB --> PMT yields
    vertex_to_TPB_to_PMT_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_1275,
                                     UV_absorption_length, vis_absorption_length, photon_vertex,
                                     pmt_parsed_information_, locations_to_check_information_, all_pmt_yields_map, keys_to_fit);
    vertex_to_TPB_to_PMT_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_511,
                                     UV_absorption_length, vis_absorption_length, electron_vertex,
                                     pmt_parsed_information_, locations_to_check_information_, all_pmt_yields_map, keys_to_fit);
    vertex_to_TPB_to_PMT_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_511,
                                     UV_absorption_length, vis_absorption_length, positron_vertex,
                                     pmt_parsed_information_, locations_to_check_information_, all_pmt_yields_map, keys_to_fit);
    // ok done!
}


template<typename T>
void FrameThread(std::atomic<bool> & running,
                 std::vector<CCMPMTKey> const & keys_to_fit,
                 std::vector<HESodiumEvent> const & sodium_events,
                 double const & c_cm_per_nsec,
                 double const & uv_index_of_refraction,
                 double const & vis_index_of_refraction,
                 T UV_absorption_length,
                 T vis_absorption_length,
                 std::vector<yields_pmt_info> const & pmt_parsed_information_,
                 std::vector<secondary_loc_info> const & locations_to_check_information_,
                 double max_time,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields,
                 std::shared_ptr<std::map<CCMPMTKey, std::vector<T>>> & this_event_binned_yields_squared) {

    // now let's loop over each event in our vector sodium_events
    for (size_t event_it = 0; event_it < sodium_events.size(); event_it ++){

        // call simulation code
        std::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> this_event_pmt_yields_map = std::make_shared<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>>();
        PutSimulationStepsTogether(sodium_events.at(event_it), c_cm_per_nsec, uv_index_of_refraction,
                                   vis_index_of_refraction, UV_absorption_length, vis_absorption_length,
                                   pmt_parsed_information_, locations_to_check_information_, this_event_pmt_yields_map, keys_to_fit);
        
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
            std::vector<photon_yield_summary<T>> const & yields = i->second;
            if(yields.size() == 0) {
                continue;
            }
            for(size_t yield_it = 0; yield_it < yields.size(); ++yield_it) {
                photon_yield_summary<T> const & yield = yields.at(yield_it);
                size_t bin_idx;
                if constexpr (std::is_same<T, double>::value) {
                    bin_idx = yield.time / 2.0;
                } else {
                    bin_idx = yield.time.value() / 2.0;
                }
                if(bin_idx >= n_bins) {
                    continue;
                }
                this_event_binned_yields->at(key).at(bin_idx) += yield.yield;
                this_event_binned_yields_squared->at(key).at(bin_idx) += (yield.yield * yield.yield);
            }
        }
    }

    running.store(false);
}

template<typename T>
void RunFrameThread(ctpl::thread_pool & pool,
                    PhotonPropagationJob<T> * job,
                    std::vector<CCMPMTKey> const & keys_to_fit,
                    double const & c_cm_per_nsec,
                    double const & uv_index_of_refraction,
                    double const & vis_index_of_refraction,
                    T UV_absorption_length,
                    T vis_absorption_length,
                    double max_time,
                    std::vector<yields_pmt_info> const & pmt_parsed_information_,
                    std::vector<secondary_loc_info> const & locations_to_check_information_) {

    job->running.store(true);
    pool.push([ &running = job->running, &keys_to_fit, &sodium_event = job->this_event,
                &c_cm_per_nsec, &uv_index_of_refraction, &vis_index_of_refraction,
                UV_absorption_length, vis_absorption_length, &pmt_parsed_information_, &locations_to_check_information_,
                max_time, job 
    ] (int id) {
    FrameThread(running,
                keys_to_fit,
                sodium_event,
                c_cm_per_nsec,
                uv_index_of_refraction,
                vis_index_of_refraction,
                UV_absorption_length, 
                vis_absorption_length,
                pmt_parsed_information_,
                locations_to_check_information_,
                max_time,
                job->this_event_binned_yields,
                job->this_event_binned_yields_squared);
    });

}

template<typename T> void YieldsPerPMT::GetAllYields(size_t n_threads, boost::shared_ptr<HESodiumEventSeries> event_vertices, T UV_absorption_length, std::vector<CCMPMTKey> const & keys_to_fit,
                                                     double max_time, std::map<CCMPMTKey, std::vector<T>> & binned_yields, std::map<CCMPMTKey, std::vector<T>> & binned_square_yields){

    // this function takes the list of event vertices, uv absorption length (which we will be fitting for), keys to fit (which will usually be a list of just one key),
    // and finally binned_yields and binned_square_yields as references which we will be updating
    // we are going to multi thread this code -- chunk up event_vertices into vector for each thread
    
    // set our absorption length for visible light (in cm -- shouldnt really affect things)
    T vis_absorption_length_ = 2000.0;

    // will be used for multi-threading our simulation jobs
    std::deque<PhotonPropagationJob<T> *> free_jobs;
    std::deque<PhotonPropagationJob<T> *> running_jobs;
    std::deque<PhotonPropagationResult<T>> results;
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
    std::vector<std::vector<HESodiumEvent>> events_per_thread;
    size_t n_events_per_thread = event_vertices->size() / num_threads;
    size_t left_over_events = event_vertices->size() - (num_threads * n_events_per_thread);
    
    size_t sodium_event_idx = 0;
    size_t accounted_for_left_over_events = 0;
    for (size_t thread_it = 0; thread_it < num_threads; thread_it++){
        // make vector to hold events
        std::vector<HESodiumEvent> this_thread_vector_of_events;

        // things to keep track of
        size_t n_events_on_this_thread = 0;

        // now time to actaully add the sodium events to this_thread_vector_of_events
        while (n_events_on_this_thread < n_events_per_thread){
            this_thread_vector_of_events.push_back(event_vertices->at(sodium_event_idx));
            n_events_on_this_thread += 1;
            sodium_event_idx += 1;
        }

        // check if we've dolled out all of the left over events yet
        if (accounted_for_left_over_events < left_over_events){
            this_thread_vector_of_events.push_back(event_vertices->at(sodium_event_idx));
            sodium_event_idx += 1;
            accounted_for_left_over_events += 1;
        }

        events_per_thread.push_back(this_thread_vector_of_events);
    }

    // now let's loop over our pre-made lists of sodium events for each thread
    for (size_t sodium_it = 0; sodium_it < events_per_thread.size(); ++sodium_it) {
        std::vector<HESodiumEvent> this_event = events_per_thread.at(sodium_it);

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
				    PhotonPropagationJob<T> * job = running_jobs.at(i);
                    running_jobs.erase(running_jobs.begin() + i);
                    free_jobs.push_back(job);
                    results.at(job->event_idx - min_vertex_idx).done = true;
                } else {
                    PhotonPropagationJob<T> * job = running_jobs.at(i);
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
            PhotonPropagationJob<T> * job = nullptr;

            if(free_jobs.size() > 0) {
                job = free_jobs.front();
                job->running.store(false);
                free_jobs.pop_front();
            } else if(running_jobs.size() < num_threads) {
                job = new PhotonPropagationJob<T>();
                job->running.store(false);
            }

            if(job != nullptr and results.size() < max_cached_vertices) {
                job->running.store(true);
                running_jobs.push_back(job);
                job->event_idx = sodium_it;
                job->this_event = this_event; 
                job->this_event_binned_yields = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>(); 
                job->this_event_binned_yields_squared = std::make_shared<std::map<CCMPMTKey, std::vector<T>>>(); 
                results.emplace_back();
                results.back().event_idx = job->event_idx;
                results.back().this_event_binned_yields = job->this_event_binned_yields;
                results.back().this_event_binned_yields_squared = job->this_event_binned_yields_squared;
                results.back().done = false;
                RunFrameThread<T>(pool, job, keys_to_fit, c_cm_per_nsec_, uv_index_of_refraction_, vis_index_of_refraction_, UV_absorption_length, vis_absorption_length_,
                               max_time, pmt_parsed_information_, locations_to_check_information_);
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
                PhotonPropagationJob<T> * job = running_jobs.at(i);
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
        for(PhotonPropagationJob<T> * obj : free_jobs){
            delete obj;
        }
    }

}

#endif
