#ifndef PhotonPropagation_H
#define PhotonPropagation_H

#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

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

#include <icetray/ctpl.h>
#include <icetray/open.h>
#include <icetray/I3Frame.h>
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
#include "icetray/robust_statistics.h"
#include <dataclasses/geometry/CCMGeometry.h>

struct PhotonPropagationJob {
    std::atomic<bool> running = false;
    std::thread thread;
    size_t thread_index = 0;
    std::vector<double>* vertex = nullptr;
    size_t vertex_index = 0;
    std::vector<std::vector<double>>* binned_charges = nullptr; // this is where we save the binned charges for each pmt for each event
};

struct PhotonPropagationResult {
    size_t vertex_index = 0;
    std::vector<std::vector<double>>* binned_charges = nullptr; // this is where we save the binned charges for each pmt for each event
    bool done = false;
};

class PhotonPropagation {
    std::exception_ptr teptr = nullptr;
    bool geo_seen = false;
    std::string geometry_name_ = std::string("CCMGeometry");
    double smearing_mu_ = 1.0;
    double smearing_sigma_ = 1.0;
    double smearing_xi_ = 0.5;
    size_t n_convolution_chunks_ = 200.0;
    double desired_chunk_width_ = 5.0;
    double desired_chunk_height_ = 5.0;
    double n_chunks_top_ = 50.0;
    double portion_light_reflected_by_tpb_ = 1.0;
    size_t n_events_to_simulate_ = (size_t)1000;
    double visible_absorption_length_ = 2000.0;

    // place to store relevant information about our pmts!!!
    // pmt_x_loc, pmt_y_loc, pmt_z_loc, facing direction_x, facing_direction_y, facing_direction_z, coating flag, pmt facing area, pmt side area
    std::vector<std::vector<double>> pmt_parsed_information_;

    // similar place to store relevant information about our secondary locations
    // for our secondary locations, I think we can probably also pre-compute the yield and travel time for 1 photon from each location to each pmt
    // x_loc, y_loc, z_loc, facing_direction_x, facing_direction_y, facing_direction_z, TPB portion, facing area, side area
    std::vector<std::vector<double>> locations_to_check_information_;
    std::vector<std::vector<double>> locations_to_check_to_pmt_yield_;
    std::vector<std::vector<double>> locations_to_check_to_pmt_travel_time_;

    // place to store list of vertices to simulate
    std::vector<std::vector<double>> verticies_to_simuate_;
    unsigned int coated_omtype = (unsigned int)10;
    unsigned int uncoated_omtype = (unsigned int)20;

    // defining some geomtry things for modelling the detector...not the most elegant
    double pmt_radius = 10.16; //radius in cm^2
    double pmt_facing_area = M_PI * std::pow(pmt_radius, 2);
    double pmt_side_area_factor = 0.217549;
    double pmt_side_area = pmt_facing_area * pmt_side_area_factor;
    double cylinder_max_x = 96.0;
    double cylinder_min_x = - cylinder_max_x;
    double cylinder_max_y = 96.0;
    double cylinder_min_y = - cylinder_max_y;
    double cylinder_max_z = 58.0;
    double cylinder_min_z = - cylinder_max_z;
    double cylinder_radius = cylinder_max_x;
    double cylinder_circumference = M_PI * 2 * cylinder_max_x;
    double cylinder_height = cylinder_max_z * 2;
    double chunk_side_area_factor = 0.1; // this number is a guess..maybe model it one day

    // now some geometry things for throwing source events
    double rod_diameter = 1.0;
    double source_diameter = 0.8;
    double rod_width = (rod_diameter - source_diameter)/2;
    double source_inset = -0.25;
    double decay_constant = 10.0;
    double source_rod_lower_end_cap = - source_inset;
    double detector_lower_end_cap = cylinder_min_z;
    double detector_radius = cylinder_radius;
    double pos_rad = source_diameter / 2;
    size_t total_events_that_escaped = 0;

    // pmt noise rate
    double noise_photons = 1.0;
    double noise_triggers = 5.0;
    double digitization_time = 16 * std::pow(10, 3); //16 usec in nsec
    double noise_rate = noise_photons / (noise_triggers * digitization_time); // units of photons/nsec
    double noise_rate_per_time_bin = 2.0 * noise_rate; // 2nsec binning

    // some constants we use for the simulation
    double c = 2.998 * std::pow(10, 8); // speed of light in m/s
    double c_cm_per_nsec = c * std::pow(10, -7); // speed of light in cm/nsec
    double uv_index_of_refraction = 1.358;
    double vis_index_of_refraction = 1.23;
    double pmt_quantum_efficiency = 0.25;
    double full_acceptance = 4.0 * M_PI;
    size_t n_pmts_to_simulate = (size_t) 0;

    size_t num_threads = std::thread::hardware_concurrency();
    size_t max_cached_vertices = (size_t) 2000;

    std::deque<PhotonPropagationJob *> free_jobs;
    std::deque<PhotonPropagationJob *> running_jobs;
    std::deque<PhotonPropagationResult> results;

public:
    PhotonPropagation();
    void Geometry(I3FramePtr frame);
    I3Vector<I3Vector<double>> GetSimulation(double const & singlet_ratio_,
                                             double const & triplet_ratio_,
                                                   double const & singlet_tau_,
                                                   double const & triplet_tau_,
                                                   double const & recombination_tau_,
                                                   double const & TPB_ratio_,
                                                   double const & TPB_tau_,
                                                   double const & UV_absorption_length_,
                                                   double const & n_photons_produced_);
};
#endif
