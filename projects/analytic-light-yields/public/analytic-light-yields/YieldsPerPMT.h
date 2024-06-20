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

public:
    YieldsPerPMT() = default;
    YieldsPerPMT(I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height);

    void SetPMTInformation(I3FramePtr frame);
    void SetGeoFrame(I3FramePtr geo_frame);
    void SetChunks(double chunk_width, double chunk_height);

    template<typename T> boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> GetAllYields(HESodiumEvent const & event_vertex, T UV_absorption_length, std::vector<CCMPMTKey> const & keys_to_fit);

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
                                                    boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> & all_pmt_yields_map,
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
                                                           boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> & all_pmt_yields_map,
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
                                                     boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> & all_pmt_yields_map,
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

template<typename T> boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> YieldsPerPMT::GetAllYields(HESodiumEvent const & event_vertex, T UV_absorption_length, std::vector<CCMPMTKey> const & keys_to_fit) {

    // so we have an event we want to simulate
    // probably want to to multi-thread at some point, but for now we will just loop over all high energy sodium events

    // let's also make our photon_yield_summary map to save to
    boost::shared_ptr<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> all_pmt_yields_map_ = boost::make_shared<std::map<CCMPMTKey, std::vector<photon_yield_summary<T>>>> ();

    // ok we've set up the geometry stuff, now we can take the vertex and calculate light yields in each pmt
    T vis_absorption_length_ = 2000.0;
    PutSimulationStepsTogether(event_vertex, c_cm_per_nsec_, uv_index_of_refraction_,
                               vis_index_of_refraction_, UV_absorption_length, vis_absorption_length_,
                               pmt_parsed_information_, locations_to_check_information_, all_pmt_yields_map_, keys_to_fit);

    return all_pmt_yields_map_;
}

#endif
