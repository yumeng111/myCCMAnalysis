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
#include "icetray/robust_statistics.h"
#include <dataclasses/geometry/CCMGeometry.h>
#include <daqtools/PhotonPropagation.h>

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

class PhotonPropagation: public I3Module {
    std::exception_ptr teptr = nullptr;
    bool geo_seen;
    std::string geometry_name_;
    std::string output_name_;
    double singlet_ratio_;
    double triplet_ratio_;
    double singlet_tau_;
    double triplet_tau_;
    double recombination_tau_;
    double TPB_ratio_;
    double TPB_tau_;
    double smearing_mu_;
    double smearing_sigma_;
    double smearing_xi_;
    size_t n_convolution_chunks_;
    double desired_chunk_width_;
    double desired_chunk_height_;
    double n_chunks_top_;
    double portion_light_reflected_by_tpb_;
    size_t n_events_to_simulate_;
    double UV_absorption_length_;
    double visible_absorption_length_;
    double n_photons_produced_;

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

    size_t num_threads;
    size_t max_cached_vertices;

    std::deque<PhotonPropagationJob *> free_jobs;
    std::deque<PhotonPropagationJob *> running_jobs;
    std::deque<PhotonPropagationResult> results;
public:
    PhotonPropagation(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void Physics(I3FramePtr frame);
};

I3_MODULE(PhotonPropagation);

PhotonPropagation::PhotonPropagation(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string("CCMGeometry"));
        AddParameter("R_s", "ratio of singlet light", (double)0.23);
        AddParameter("R_t", "ratio of triplet light", (double)0.71);
        AddParameter("tau_s", "singlet time constant", (double)8.2);
        AddParameter("tau_t", "triplet time constant", (double)1445.0);
        AddParameter("tau_rec", "recombination time constant", (double)75.5);
        AddParameter("R_TPB", "ratio of TPB light", (double)0.1);
        AddParameter("tau_TPB", "TPB time constant", (double)2.0);
        AddParameter("mu", "mu for GEV function to smear light", (double)1.0);
        AddParameter("sigma", "sigma for GEV function to smear light", (double)1.0);
        AddParameter("xi", "xi for GEV function to smear light", (double)0.5);
        AddParameter("NConvolutionChunks", "number of chunks to convolve over", (size_t)200);
        AddParameter("ChunkWidth", "width to chunk up the sides of the detector for light propagation", 5.0); // used to be 5
        AddParameter("ChunkHeight", "height to chunk up the sides of the detector for light propagation", 5.0); // used to be 5
        AddParameter("NChunksTop", "number of chunks to split top and bottom faces of detector into", 50.0); // used to be 50
        AddParameter("TPBPortion", "portion of light that tpb reflects", 1.0);
        AddParameter("NEventsToSimulate", "number of events we want to simulate for sodium source like ensemble", (size_t)1000);
        AddParameter("UVAbsorptionLength", "how far UV light travels in cm before 1/e is absorped", 40.0);
        AddParameter("VisAbsorptionLength", "how far visible light travels in cm before 1/e is absorped", 2000.0);
        AddParameter("NPhotonsProduced", "how many photons are produced in this event (~40,000 photons/MeV is the hope)", 40000.0);
        AddParameter("OutputName", "Key to save output CCMWaveformDoubleSeries to", std::string("CCMCalibratedWaveforms"));
        AddParameter("NumThreads", "Number of worker threads to use for baseline estimation", (size_t)0);
        AddParameter("MaxCachedVertices", "The maximum number of vertices this module is allowed to have cached", (size_t)(2000));
}


void  PhotonPropagation::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("R_s", singlet_ratio_);
    GetParameter("R_t", triplet_ratio_);
    GetParameter("tau_s", singlet_tau_);
    GetParameter("tau_t", triplet_tau_);
    GetParameter("tau_rec", recombination_tau_);
    GetParameter("R_TPB", TPB_ratio_);
    GetParameter("tau_TPB", TPB_tau_);
    GetParameter("mu", smearing_mu_);
    GetParameter("sigma", smearing_sigma_);
    GetParameter("xi", smearing_xi_);
    GetParameter("NConvolutionChunks", n_convolution_chunks_);
    GetParameter("ChunkWidth", desired_chunk_width_);
    GetParameter("ChunkHeight", desired_chunk_height_);
    GetParameter("NChunksTop", n_chunks_top_);
    GetParameter("TPBPortion", portion_light_reflected_by_tpb_);
    GetParameter("NEventsToSimulate", n_events_to_simulate_);
    GetParameter("UVAbsorptionLength", UV_absorption_length_);
    GetParameter("VisAbsorptionLength", visible_absorption_length_);
    GetParameter("NPhotonsProduced", n_photons_produced_);
    GetParameter("OutputName", output_name_);
    GetParameter("NumThreads", num_threads);
    GetParameter("MaxCachedVertices", max_cached_vertices);
    if(num_threads == 0) {
        size_t const processor_count = std::thread::hardware_concurrency();
        num_threads = processor_count;
    }
}

void get_ray_intersections(double const & x,
                           double const & y,
                           double const & z,
                           double const & nx,
                           double const & ny,
                           double const & nz,
                           double const & cz2,
                           double & return_x,
                           double & return_y,
                           double & return_z){
    double nx2 = std::pow(nx, 2);
    double ny2 = std::pow(ny, 2);
    double nr2 = nx2 + ny2;
    double nr = std::sqrt(nr2);
    double r0_2 = std::pow(x, 2) + std::pow(y, 2);
    double r0 = std::sqrt(r0_2);

    // now let's check for intersection with lower z cap
    double t1 = (cz2 - z) / nz;
    double xx = x + nx * t1;
    double yy = y + ny * t1;
    double zz = cz2;
    return_x = xx;
    return_y = yy;
    return_z = zz;

}

void check_escape(double const & x,
                  double const & y,
                  double const & z,
                  double const & nx,
                  double const & ny,
                  double const & nz,
                  double const & source_rod_lower_end_cap,
                  double const & detector_radius,
                  double const & decay_constant,
                  double const & pos_rad,
                  double const & detector_lower_end_cap,
                  std::mt19937 & gen,
                  std::uniform_real_distribution<double> & dis_0_1,
                  bool & escaped,
                  double & final_x,
                  double & final_y,
                  double & final_z){

    // first let's get our intersection
    double x1;
    double y1;
    double z1;
    get_ray_intersections(x, y, z, nx, ny, nz, source_rod_lower_end_cap, x1, y1, z1);

    // need to check the radius of intersection
    double intersection_radius = std::sqrt(std::pow(x1, 2) + std::pow(y1, 2));

    // Generate a random number
    double cdf_prob = dis_0_1(gen);
    // distance travelled
    double dist = -decay_constant * std::log(cdf_prob + 1/decay_constant);

    if (intersection_radius < pos_rad){

        // now let's check if inside the detector
        final_x = x + dist * nx;
        final_y = y + dist * ny;
        final_z = z + dist * nz;

        if (-detector_radius < final_x  and final_x < detector_radius and -detector_radius < final_y  and final_y < detector_radius and detector_lower_end_cap < final_z){
            escaped = true;
            if (final_z > 0){
                escaped = false;
            }
        }
    }

    if (intersection_radius >= pos_rad){
        // ok so we are a sodium photon going up through the rod... about 1/3 chance of survival (this is pretty ad hoc..)
        double survival = dis_0_1(gen);
        if (survival > 1.0/3.0){
            escaped = false;
        }
        else{
            // now let's check if inside the detector
            final_x = x + dist * nx;
            final_y = y + dist * ny;
            final_z = z + dist * nz;

            if (-detector_radius < final_x  and final_x < detector_radius and -detector_radius < final_y  and final_y < detector_radius and detector_lower_end_cap < final_z){
                escaped = true;
            }
        }
    }

}

void get_1275kev_photon(double const & source_diameter,
                        double const & source_rod_lower_end_cap,
                        double const & pos_rad,
                        double const & decay_constant,
                        double const & detector_radius,
                        double const & detector_lower_end_cap,
                        std::mt19937 & gen,
                        std::uniform_real_distribution<double> & dis_angle,
                        std::uniform_real_distribution<double> & dis_0_1,
                        std::uniform_real_distribution<double> & dis_neg_1_1,
                        bool & escaped,
                        double & final_x,
                        double & final_y,
                        double & final_z){
    // let's get an initial position
    double theta_pos = dis_angle(gen);
    double r = std::sqrt(dis_0_1(gen)) * source_diameter/2;
    double x = r * std::cos(theta_pos);
    double y = r * std::sin(theta_pos);
    double z = 0;

    // now let's get an initial direction
    double phi = dis_angle(gen);
    double cos_theta = dis_neg_1_1(gen);
    double theta = std::acos(cos_theta);
    double nx = std::cos(phi) * std::sin(theta);
    double ny = std::sin(phi) * std::sin(theta);
    double nz = std::cos(theta);

    // now let's see if the event escapes the rod
    check_escape(x, y, z, nx, ny, nz, source_rod_lower_end_cap, detector_radius, decay_constant, pos_rad, detector_lower_end_cap, gen, dis_0_1, escaped, final_x, final_y, final_z);
}

void get_solid_angle_and_distance_vertex_to_location(double const & vertex_x,
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
                                                     double & D){
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

void get_light_yield(double const & omega,
                     double const & total_omega,
                     double const & distance_travelled,
                     double const & normalization,
                     double const & absorption_length,
                     double const & efficiency,
                     double & light_yield){

    light_yield = efficiency * (omega / total_omega) * normalization * std::exp(- distance_travelled / absorption_length);
}

void secondary_loc_to_pmt_propagation(double const & full_acceptance,
                                      double const & c_cm_per_nsec,
                                      double const & vis_index_of_refraction,
                                      double const & quantum_efficiency,
                                      double const & vis_absorption_length,
                                      std::vector<std::vector<double>> const & pmt_parsed_information_,
                                      double const & loc_x,
                                      double const & loc_y,
                                      double const & loc_z,
                                      std::vector<double> & pmt_photon_yields,
                                      std::vector<double> & pmt_photon_propagation_times){

    // we pass this function two vectors to store yields and time offsets for charge in each ptm
    // let's define some things
    bool is_visible = true; // we are propagating wavelengthshifted visible light from TPB locs to PMTs...always visible!
    double n_photons_produced = 1.0; // we are only propagating 1 photons! so we can mulitply by actually number of photons later on!

    double efficiency;
    double omega;
    double distance_travelled;
    double photons_in_this_pmt;
    double travel_time;
    double pmt_x_loc;
    double pmt_y_loc;
    double pmt_z_loc;
    double facing_dir_x;
    double facing_dir_y;
    double facing_dir_z;
    double coating_flag;
    double pmt_facing_area;
    double pmt_side_area;

    // let's start by looping over our pmt parsed information
    for (size_t pmt_it = 0; pmt_it < pmt_parsed_information_.size(); pmt_it ++){
        efficiency = 1.0; // this describes light that goes into pmts
                          // for vis light on coated, it's 0.5
                          // for vis light on unocated, it's 1.0

        pmt_x_loc = pmt_parsed_information_[pmt_it][0];
        pmt_y_loc = pmt_parsed_information_[pmt_it][1];
        pmt_z_loc = pmt_parsed_information_[pmt_it][2];
        facing_dir_x = pmt_parsed_information_[pmt_it][3];
        facing_dir_y = pmt_parsed_information_[pmt_it][4];
        facing_dir_z = pmt_parsed_information_[pmt_it][5];
        coating_flag = pmt_parsed_information_[pmt_it][6]; // coating flag == 1.0 for coated pmts and 0.0 for uncoated!
        pmt_facing_area = pmt_parsed_information_[pmt_it][7];
        pmt_side_area = pmt_parsed_information_[pmt_it][8];

        // let's get solid angle and distance from loc to this pmt
        get_solid_angle_and_distance_vertex_to_location(loc_x, loc_y, loc_z,
                                                     pmt_x_loc, pmt_y_loc, pmt_z_loc,
                                                     facing_dir_x, facing_dir_y, facing_dir_z,
                                                     pmt_facing_area, pmt_side_area, omega, distance_travelled);

       // now let's do a check for the efficiency
        if (is_visible and coating_flag == 1.0){
            // visible light on coated pmts
            efficiency = 0.5;
        }

        // now call function to get light yields
        get_light_yield(omega, full_acceptance, distance_travelled, n_photons_produced, vis_absorption_length, efficiency, photons_in_this_pmt);
        photons_in_this_pmt *= quantum_efficiency;

        // ok so we have the photons seen by this pmt, let's get the propagation times
        travel_time = distance_travelled / (c_cm_per_nsec / vis_index_of_refraction); // units of nsec

        // now all that's left to do is save!
        pmt_photon_yields.push_back(photons_in_this_pmt);
        pmt_photon_propagation_times.push_back(travel_time);
    }

}
PhotonPropagation::PhotonPropagation(){}

void PhotonPropagation::SetData(I3Vector<I3Vector<double>> data_series){
    data_series_ = data_series;
    // data_series_ is index by n_pmts x n_time_bins
    // while we are setting the data, let's also set the time
    // we want the peak of the data time to be at 0 and 2 nsec binning

    times_of_data_points_.clear();
    for (size_t second_dim_it = 0; second_dim_it < data_series_[0].size(); second_dim_it ++ ){
        times_of_data_points_.push_back((double)second_dim_it * 2.0);
    }

    max_data_value_ = 0;
    time_of_max_data_value_ = 0;

    // let's first find the time of max, then set that equal to zero
    double total_charge_per_time_bin;
    for (size_t time_bin_it = 0; time_bin_it < times_of_data_points_.size(); time_bin_it ++){
        total_charge_per_time_bin = 0;
        for (size_t pmt_it = 0; pmt_it < data_series_.size(); pmt_it ++){
            total_charge_per_time_bin += data_series[pmt_it][time_bin_it];
        }
        if (total_charge_per_time_bin > max_data_value_){
            max_data_value_ = total_charge_per_time_bin;
            time_of_max_data_value_ = times_of_data_points_[time_bin_it];
        }
    }

    // ok now we can subtract time_of_max_data_value_ from times_of_data_points_s
    for (size_t time_bin_it = 0; time_bin_it < times_of_data_points_.size(); time_bin_it ++){
        times_of_data_points_[time_bin_it] -= time_of_max_data_value_;
    }


}

void PhotonPropagation::SetDataSampleSize(size_t n_data_samples){
    n_data_samples_ = n_data_samples;
}

void PhotonPropagation::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    I3Map<CCMPMTKey, CCMOMGeo> const & pmt_geo = geo.pmt_geo;
    std::vector<double> this_pmt_info (9);
    std::vector<double> this_loc_info (9);
    std::vector<double> loc_to_pmt_photon_yields;
    std::vector<double> loc_to_pmt_photon_propagation_times;

    for(std::pair<CCMPMTKey const, CCMOMGeo> const & it : pmt_geo) {
        CCMPMTType omtype = it.second.omtype;
        double coating_flag;
        if (omtype == coated_omtype){
            coating_flag = 1.0; // coating flag = 1 ==> pmt is coated!
        }
        else if (omtype == uncoated_omtype){
            coating_flag = 0.0; // uncoated pmt!
        }
        else{
            continue; // this pmt is not an 8in! we dont need it!
        }
        // now let's compute our facing direction
        I3Position position = it.second.position;
        double pos_x = position.GetX();
        double pos_y = position.GetY();
        double pos_z = position.GetZ();

        double facing_dir_x;
        double facing_dir_y;
        double facing_dir_z;

        if (pos_z == 58.0){
            // region 0 pmts!
            facing_dir_x = 0.0;
            facing_dir_y = 0.0;
            facing_dir_z = -1.0;
        }
        else if (pos_z == -58.0){
            // region 6 pmts!
            facing_dir_x = 0.0;
            facing_dir_y = 0.0;
            facing_dir_z = 1.0;
        }
        else{
            double facing_radius = std::pow(pos_x, 2) + std::pow(pos_y, 2);
            double facing_r_x = - pos_x / facing_radius;
            double facing_r_y = - pos_y / facing_radius;
            double facing_dir_norm_factor = std::sqrt(std::pow(facing_r_x, 2) + std::pow(facing_r_y, 2));
            facing_dir_x = facing_r_x / facing_dir_norm_factor;
            facing_dir_y = facing_r_y / facing_dir_norm_factor;
            facing_dir_z = 0.0;
        }

        // now time save!!!
        this_pmt_info.clear();
        this_pmt_info.push_back(pos_x);
        this_pmt_info.push_back(pos_y);
        this_pmt_info.push_back(pos_z);
        this_pmt_info.push_back(facing_dir_x);
        this_pmt_info.push_back(facing_dir_y);
        this_pmt_info.push_back(facing_dir_z);
        this_pmt_info.push_back(coating_flag);
        this_pmt_info.push_back(pmt_facing_area);
        this_pmt_info.push_back(pmt_side_area);
        pmt_parsed_information_.push_back(this_pmt_info);
        n_pmts_to_simulate += (size_t) 1;

    }
    // so we've parsed our pmt info, but now we need to get our secondary locations
    // let's start with chunking up the sides of the detector
    size_t n_chunks_c = (size_t) cylinder_circumference / desired_chunk_width_;
    size_t n_chunks_z = (size_t) cylinder_height / desired_chunk_height_;

    std::vector<double> possible_circumference_positions(n_chunks_c);
    for (size_t i = 0; i < n_chunks_c; i++){
        double this_circ = cylinder_circumference * i / n_chunks_c;
        possible_circumference_positions[i] = this_circ;
    }

    std::vector<double> possible_z_positions(n_chunks_z - 1);
    for (size_t i = 1; i < n_chunks_z; i++){
        double this_z = cylinder_min_z + ((cylinder_max_z - cylinder_min_z) * i / n_chunks_z);
        possible_z_positions[i-1] = this_z;
    }

    double area_side_chunks = std::abs((possible_circumference_positions[1] - possible_circumference_positions[0]) * (possible_z_positions[1] - possible_z_positions[0]));

    // now let's calculate x and y from these circumference positions
    for (size_t i = 0; i < possible_circumference_positions.size(); i++){
        double loc_theta = possible_circumference_positions[i] / cylinder_radius;
        double loc_x = cylinder_radius * std::cos(loc_theta);
        double loc_y = cylinder_radius * std::sin(loc_theta);

        // now let's iterate over possible z positions
        for (size_t j = 0; j <possible_z_positions.size(); j++){
            double loc_z = possible_z_positions[j];

            // now we need to check if we're on an uncoated pmt...if so we do NOT save
            bool save_this_loc = true;
            double pmt_portion = portion_light_reflected_by_tpb_;
            for (size_t pmt_it = 0; pmt_it < pmt_parsed_information_.size(); pmt_it ++){
                std::vector<double> this_pmt_info = pmt_parsed_information_[pmt_it];
                double pmt_coating_flag = this_pmt_info[6];

                double pmt_x = this_pmt_info[0];
                double pmt_y = this_pmt_info[1];
                double pmt_z = this_pmt_info[2];

                // now let's calculate distance from this location to this pmt
                double loc_dist_to_pmt = std::sqrt(std::pow(pmt_x - loc_x, 2) + std::pow(pmt_y - loc_y, 2) + std::pow(pmt_z - loc_z, 2));

                // now check if we're on top a pmt
                if (loc_dist_to_pmt <= pmt_radius){
                        // oops! we're on a pmt! if it's an uncoated pmt, save_this_loc = false, if it's a coated pmt, TPB portion = 0.5
                        if (pmt_coating_flag == 1.0){
                            pmt_portion *= 0.5;
                        }
                        else{
                            save_this_loc = false;
                        }
                    }
                }

            // so we have a possible loc and we've finished check if it's on top an uncoated pmt...let's save (maybe)!
            if (save_this_loc){
                double facing_radius = std::pow(loc_x, 2) + std::pow(loc_y, 2);
                double facing_r_x = - loc_x / facing_radius;
                double facing_r_y = - loc_y / facing_radius;
                double facing_dir_norm_factor = std::sqrt(std::pow(facing_r_x, 2) + std::pow(facing_r_y, 2));
                double facing_dir_x = facing_r_x / facing_dir_norm_factor;
                double facing_dir_y = facing_r_y / facing_dir_norm_factor;
                double facing_dir_z = 0.0;

                // now a vector to save things to
                this_loc_info.clear();
                this_loc_info.push_back(loc_x);
                this_loc_info.push_back(loc_y);
                this_loc_info.push_back(loc_z);
                this_loc_info.push_back(facing_dir_x);
                this_loc_info.push_back(facing_dir_y);
                this_loc_info.push_back(facing_dir_z);
                this_loc_info.push_back(pmt_portion);
                this_loc_info.push_back(area_side_chunks);
                this_loc_info.push_back(area_side_chunks * chunk_side_area_factor);
                locations_to_check_information_.push_back(this_loc_info);

                // we are also going to pre-compute the light yield from 1 visible photon from this loc to every pmt
                // as well as the travel time for that 1 photon to go from this loc to every pmt
                loc_to_pmt_photon_yields.clear();
                loc_to_pmt_photon_propagation_times.clear();
                secondary_loc_to_pmt_propagation(full_acceptance, c_cm_per_nsec, vis_index_of_refraction, pmt_quantum_efficiency,
                                                 visible_absorption_length_, pmt_parsed_information_, loc_x, loc_y, loc_z,
                                                 loc_to_pmt_photon_yields, loc_to_pmt_photon_propagation_times);
                locations_to_check_to_pmt_yield_.push_back(loc_to_pmt_photon_yields);
                locations_to_check_to_pmt_travel_time_.push_back(loc_to_pmt_photon_propagation_times);
            }
        }

    }


    // time to chunk up the top and bottom of the detector
    std::vector<double> possible_x_positions(n_chunks_top_ - 1);
    for (size_t i = 1; i < n_chunks_top_; i++){
        double x = cylinder_min_x + ((cylinder_max_x - cylinder_min_x) * i / n_chunks_top_);
        possible_x_positions[i-1] = x;
    }

    std::vector<double> possible_y_positions(n_chunks_top_ - 1);
    for (size_t i = 1; i < n_chunks_top_; i++){
        double y = cylinder_min_y + ((cylinder_max_y - cylinder_min_y) * i / n_chunks_top_);
        possible_y_positions[i-1] = y;
    }

    double area_face_chunks = std::abs((possible_x_positions[1] - possible_x_positions[0]) * (possible_y_positions[1] - possible_y_positions[0]));

    // now loop over our possible x and y pos
    for (size_t x_it = 0; x_it < possible_x_positions.size(); x_it ++){
        for (size_t y_it = 0; y_it < possible_y_positions.size(); y_it ++){
            // now let's make the sure this x, y combination is physics
            double this_x = possible_x_positions[x_it];
            double this_y = possible_y_positions[y_it];
            double radius = std::sqrt(std::pow(this_x, 2) + std::pow(this_y, 2));
            if (radius >= cylinder_max_x){
                continue;
            }

            // ok so this is a valid position on the face of the detector!
            // let's check if either the top or bottom pos is on a PMT
            double z_top = cylinder_max_z;
            double z_bottom = cylinder_min_z;

            bool save_top_loc = true;
            double pmt_portion_top = portion_light_reflected_by_tpb_;
            bool save_bottom_loc = true;
            double pmt_portion_bottom = portion_light_reflected_by_tpb_;
            for (size_t pmt_it = 0; pmt_it < pmt_parsed_information_.size(); pmt_it ++){
                std::vector<double> this_pmt_info = pmt_parsed_information_[pmt_it];
                double pmt_coating_flag = this_pmt_info[6];

                double pmt_x = this_pmt_info[0];
                double pmt_y = this_pmt_info[1];
                double pmt_z = this_pmt_info[2];

                // now let's calculate distance from this location to this pmt
                double top_loc_dist_to_pmt = std::sqrt(std::pow(pmt_x - this_x, 2) + std::pow(pmt_y - this_y, 2) + std::pow(pmt_z - z_top, 2));
                double bottom_loc_dist_to_pmt = std::sqrt(std::pow(pmt_x - this_x, 2) + std::pow(pmt_y - this_y, 2) + std::pow(pmt_z - z_bottom, 2));

                // now check if we're on top a pmt
                if (top_loc_dist_to_pmt <= pmt_radius){
                    // oops! we're on a pmt! if it's an uncoated pmt, save_this_loc = false, if it's a coated pmt, TPB portion = 0.5
                    if (pmt_coating_flag == 1.0){
                        pmt_portion_top *= 0.5;
                    }
                    else{
                        save_top_loc = false;
                    }
                }

                // now check if we're on bottom a pmt
                if (bottom_loc_dist_to_pmt <= pmt_radius){
                    // oops! we're on a pmt! if it's an uncoated pmt, save_this_loc = false, if it's a coated pmt, TPB portion = 0.5
                    if (pmt_coating_flag == 1.0){
                        pmt_portion_bottom *= 0.5;
                    }
                    else{
                        save_bottom_loc = false;
                    }
                }

            }
            // ok we've finished checking if we're on pmts...let's save
            if (save_top_loc){
                double facing_dir_x = 0.0;
                double facing_dir_y = 0.0;
                double facing_dir_z = -1.0;

                // now a vector to save things to
                this_loc_info.clear();
                this_loc_info.push_back(this_x);
                this_loc_info.push_back(this_y);
                this_loc_info.push_back(z_top);
                this_loc_info.push_back(facing_dir_x);
                this_loc_info.push_back(facing_dir_y);
                this_loc_info.push_back(facing_dir_z);
                this_loc_info.push_back(pmt_portion_top);
                this_loc_info.push_back(area_face_chunks);
                this_loc_info.push_back(area_face_chunks * chunk_side_area_factor);
                locations_to_check_information_.push_back(this_loc_info);

                // we are also going to pre-compute the light yield from 1 visible photon from this loc to every pmt
                // as well as the travel time for that 1 photon to go from this loc to every pmt
                loc_to_pmt_photon_yields.clear();
                loc_to_pmt_photon_propagation_times.clear();
                secondary_loc_to_pmt_propagation(full_acceptance, c_cm_per_nsec, vis_index_of_refraction, pmt_quantum_efficiency,
                                                 visible_absorption_length_, pmt_parsed_information_, this_x, this_y, z_top,
                                                 loc_to_pmt_photon_yields, loc_to_pmt_photon_propagation_times);
                locations_to_check_to_pmt_yield_.push_back(loc_to_pmt_photon_yields);
                locations_to_check_to_pmt_travel_time_.push_back(loc_to_pmt_photon_propagation_times);
            }

            if (save_bottom_loc){
                double facing_dir_x = 0.0;
                double facing_dir_y = 0.0;
                double facing_dir_z = 1.0;

                // now a vector to save things to
                this_loc_info.clear();
                this_loc_info.push_back(this_x);
                this_loc_info.push_back(this_y);
                this_loc_info.push_back(z_bottom);
                this_loc_info.push_back(facing_dir_x);
                this_loc_info.push_back(facing_dir_y);
                this_loc_info.push_back(facing_dir_z);
                this_loc_info.push_back(pmt_portion_bottom);
                this_loc_info.push_back(area_face_chunks);
                this_loc_info.push_back(area_face_chunks * chunk_side_area_factor);
                locations_to_check_information_.push_back(this_loc_info);

                // we are also going to pre-compute the light yield from 1 visible photon from this loc to every pmt
                // as well as the travel time for that 1 photon to go from this loc to every pmt
                loc_to_pmt_photon_yields.clear();
                loc_to_pmt_photon_propagation_times.clear();
                secondary_loc_to_pmt_propagation(full_acceptance, c_cm_per_nsec, vis_index_of_refraction, pmt_quantum_efficiency,
                                                 visible_absorption_length_, pmt_parsed_information_, this_x, this_y, z_bottom,
                                                 loc_to_pmt_photon_yields, loc_to_pmt_photon_propagation_times);
                locations_to_check_to_pmt_yield_.push_back(loc_to_pmt_photon_yields);
                locations_to_check_to_pmt_travel_time_.push_back(loc_to_pmt_photon_propagation_times);
            }


        }
    }

    // while we're pre-computing things, let's also pre-compute verticies for our ensemble of sodium events
    // let's make some random number generators that we will pass to our functions
    std::random_device rd;
    std::mt19937 gen(rd());
    // Create a uniform distribution between 0 and 1
    std::uniform_real_distribution<double> dis_0_1(0.0, 1.0);
    // Create a uniform distribution between -1 and 1
    std::uniform_real_distribution<double> dis_neg_1_1(-1.0, 1.0);
    // Create a uniform distribution between 0 and 2pi
    std::uniform_real_distribution<double> dis_angle(0.0, 2.0*M_PI);

    // let's throw some sodium events and get a list of verticies to simulate!
    bool escaped;
    double final_x;
    double final_y;
    double final_z;
    std::vector<double> this_vertex (3);

    for (size_t event_it = 0; event_it < 10 * n_events_to_simulate_; event_it ++){
        if (total_events_that_escaped < n_events_to_simulate_){
            get_1275kev_photon(source_diameter,
                               source_rod_lower_end_cap,
                               pos_rad,
                               decay_constant,
                               detector_radius,
                               detector_lower_end_cap,
                               gen,
                               dis_angle,
                               dis_0_1,
                               dis_neg_1_1,
                               escaped,
                               final_x,
                               final_y,
                               final_z);
            if (escaped){
                total_events_that_escaped += 1;
                // let's save this vertex!
                this_vertex.clear();
                this_vertex.push_back(final_x);
                this_vertex.push_back(final_y);
                this_vertex.push_back(final_z);
                verticies_to_simuate_.push_back(this_vertex);
            }
        }
    }

    geo_seen = true;
    PushFrame(frame);
}

void LAr_scintillation_timing(double const & time, double const & R_s, double const & R_t, double const & tau_s, double const & tau_t, double const & tau_rec, double & resulting_light){
    double singlet = (R_s /tau_s) * std::exp(- time / tau_s);
    double triplet = (R_t /tau_t) * std::exp(- time / tau_t);
    double recombination = (1 - R_s - R_t) / (tau_rec * std::pow((1 + time / tau_rec), 2));
    resulting_light = singlet + triplet + recombination;

}

void TPB_emission(double const & time, double const & R_TPB, double const & tau_TPB, double & resulting_light){
    double prompt = ((1 - R_TPB) / tau_TPB) * std::exp(- time / tau_TPB);
    resulting_light = prompt;
}

void generalized_extreme_value_pdf(double const & time, double const & sigma, double const & mu, double const & xi, double & resulting_value){
    resulting_value = 0;
    double s = (time - mu) / sigma;
    if (xi == 0){
        resulting_value = std::exp(-s) * std::exp(-std::exp(-s));
    }
    if (xi != 0 and xi * s > -1){
        resulting_value = std::pow((1 + xi * s), (-(1 + 1/xi))) * std::exp(-std::pow((1 + xi * s), (-1/xi)));
    }
}

void linear_interpolation(double const & desired_time,
                        std::vector<double> const & all_times,
                        std::vector<double> const & all_data,
                        double & desired_value){
    auto it_below = std::lower_bound(all_times.begin(), all_times.end(), desired_time);

    if (it_below == all_times.begin()) {
        desired_value = all_data.front();
    }
    else if (it_below == all_times.end()) {
        desired_value = all_data.back();
    }
    else {
        int idx_below = std::distance(all_times.begin(), it_below) - 1;
        double time_below = all_times[idx_below];
        double data_value_below = all_data[idx_below];

        double time_above = all_times[idx_below + 1];
        double data_value_above = all_data[idx_below + 1];

        desired_value = data_value_below + (desired_time - time_below) * (data_value_above - data_value_below) / (time_above - time_below);
    }
}

void convolve_light_with_gev(double const & this_time_step,
                            std::vector<double> const & light_profile,
                            std::vector<double> const & light_times,
                            double const & delta,
                            size_t const & n_convolution_chunks,
                            double const & sigma,
                            double const & mu,
                            double const & xi,
                            double & result){
    double total_val = 0;
    double delta_convolution_time_step = this_time_step / n_convolution_chunks;

    for (size_t i = 0; i <= n_convolution_chunks; i++){
        double t2 = this_time_step * i / n_convolution_chunks;
        double gev_result;
        generalized_extreme_value_pdf(this_time_step - t2, sigma, mu, xi, gev_result);
        double interpolated_light_value;
        linear_interpolation(t2, light_times, light_profile, interpolated_light_value);
        double this_val = gev_result * interpolated_light_value * delta_convolution_time_step;
        total_val += this_val;
    }

    result = total_val * delta;

}

void LAr_scintillation_light_integration(double const & this_time_step,
                                        double const & delta,
                                        size_t const & n_convolution_chunks,
                                        double const & R_s,
                                        double const & R_t,
                                        double const & tau_s,
                                        double const & tau_t,
                                        double const & tau_rec,
                                        double const & R_TPB,
                                        double const & tau_TPB,
                                        double & result){

    double total_val = 0;
    double delta_convolution_time_step = this_time_step / n_convolution_chunks;

    for (size_t i = 0; i <= n_convolution_chunks; i++){
        double t2 = this_time_step * i / n_convolution_chunks;
        double TPB_result;
        TPB_emission(this_time_step - t2, R_TPB, tau_TPB, TPB_result);
        double LAr_scint_result;
        LAr_scintillation_timing(t2, R_s, R_t, tau_s, tau_t, tau_rec, LAr_scint_result);
        double light_val = TPB_result * LAr_scint_result * delta_convolution_time_step;
        total_val += light_val;
    }

    result = total_val * delta;

}

void get_total_light_profile(double const & R_s,
                            double const & R_t,
                            double const & tau_s,
                            double const & tau_t,
                            double const & tau_rec,
                            double const & R_TPB,
                            double const & tau_TPB,
                            double const & sigma,
                            double const & mu,
                            double const & xi,
                            size_t const & n_convolution_chunks,
                            std::vector<double> const & times,
                            std::vector<double> & final_light_profile){
    // times is a vector of times that we are integrating over
    double delta = times[1] - times[0];
    std::vector<double> LAr_light_profile (times.size());
    for (size_t t = 0; t < times.size(); t++){
        LAr_scintillation_light_integration(times[t], delta, n_convolution_chunks, R_s, R_t, tau_s, tau_t, tau_rec, R_TPB, tau_TPB, LAr_light_profile[t]);
    }
    // let's normalize our light profile
    double total_light_profile_value = 0;
    for (size_t i = 0; i < LAr_light_profile.size(); i++){
        total_light_profile_value += LAr_light_profile[i];
    }

    std::vector<double> normalized_LAr_light_profile (LAr_light_profile.size());
    for (size_t i = 0; i < LAr_light_profile.size(); i++){
        normalized_LAr_light_profile[i] = LAr_light_profile[i] / total_light_profile_value;
    }

    // now let's convolve our light profile with a gev
    std::vector<double> light_convolved_with_gev (times.size());
    for (size_t t = 0; t < times.size(); t++){
        convolve_light_with_gev(times[t], normalized_LAr_light_profile, times, delta, n_convolution_chunks, sigma, mu, xi, light_convolved_with_gev[t]);
    }

    // let's normalize our final light profile
    double total_light_convolved_with_gev_value = 0;
    for (size_t i = 0; i < light_convolved_with_gev.size(); i++){
        total_light_convolved_with_gev_value += light_convolved_with_gev[i];
    }

    for (size_t i = 0; i < light_convolved_with_gev.size(); i++){
        final_light_profile[i] = light_convolved_with_gev[i] / total_light_convolved_with_gev_value;
    }

}

//time for the meaty functions...actually propagating light!
void vertex_to_pmt_propagation(double const & full_acceptance,
                               double const & c_cm_per_nsec,
                               double const & uv_index_of_refraction,
                               double const & vis_index_of_refraction,
                               double const & quantum_efficiency,
                               double const & n_photons_produced,
                               double const & absorption_length,
                               std::vector<std::vector<double>> const & pmt_parsed_information_,
                               bool const & is_visible,
                               double const & vertex_x,
                               double const & vertex_y,
                               double const & vertex_z,
                               std::vector<double> & pmt_photon_yields,
                               std::vector<double> & pmt_photon_propagation_times){

    // we pass this function two vectors to store yields and time offsets for charge in each ptm
    // let's define some things
    double default_val = 0.0;
    double efficiency;
    double omega;
    double distance_travelled;
    double photons_in_this_pmt;
    double travel_time;
    double pmt_x_loc;
    double pmt_y_loc;
    double pmt_z_loc;
    double facing_dir_x;
    double facing_dir_y;
    double facing_dir_z;
    double coating_flag;
    double pmt_facing_area;
    double pmt_side_area;

    // let's start by looping over our pmt parsed information
    for (size_t pmt_it = 0; pmt_it < pmt_parsed_information_.size(); pmt_it ++){
        efficiency = 1.0; // this describes light that goes into pmts...so for UV light on coated it will be changed to 0.5
                          // for vis light on coated, it's 0.5
                          // for UV light on unocated, it's 0
                          // for vis light on unocated, it's 1.0

        pmt_x_loc = pmt_parsed_information_.at(pmt_it).at(0);
        pmt_y_loc = pmt_parsed_information_.at(pmt_it).at(1);
        pmt_z_loc = pmt_parsed_information_.at(pmt_it).at(2);
        facing_dir_x = pmt_parsed_information_.at(pmt_it).at(3);
        facing_dir_y = pmt_parsed_information_.at(pmt_it).at(4);
        facing_dir_z = pmt_parsed_information_.at(pmt_it).at(5);
        coating_flag = pmt_parsed_information_.at(pmt_it).at(6); // coating flag == 1.0 for coated pmts and 0.0 for uncoated!
        pmt_facing_area = pmt_parsed_information_.at(pmt_it).at(7);
        pmt_side_area = pmt_parsed_information_.at(pmt_it).at(8);

        if (is_visible == false and coating_flag == 0.0){
            // this is the case where we are propagating UV light and we have an uncoated pmt...so we won't see anything
            pmt_photon_yields.push_back(default_val);
            pmt_photon_propagation_times.push_back(default_val);
            continue;
        }

        else{
            // let's get solid angle and distance from vertex to this pmt
            get_solid_angle_and_distance_vertex_to_location(vertex_x, vertex_y, vertex_z,
                                                            pmt_x_loc, pmt_y_loc, pmt_z_loc,
                                                            facing_dir_x, facing_dir_y, facing_dir_z,
                                                            pmt_facing_area, pmt_side_area, omega, distance_travelled);

           // now let's do a check for the efficiency
            if ((is_visible == false) or (is_visible and coating_flag == 1.0)){
                // selecting UV light on coated pmts
                // and visible light on coated pmts
                efficiency = 0.5;
            }

            // let's also get pmt eff based on incident angle
            double pmt_ang_eff;
            get_pmt_efficiency(vertex_x, vertex_y, vertex_z,
                               pmt_x_loc, pmt_y_loc, pmt_z_loc,
                               facing_dir_x, facing_dir_y, facing_dir_z, pmt_ang_eff);
 
            efficiency *= pmt_ang_eff;

            // now call function to get light yields
            get_light_yield(omega, full_acceptance, distance_travelled, n_photons_produced, absorption_length, efficiency, photons_in_this_pmt);
            photons_in_this_pmt *= quantum_efficiency;

            // ok so we have the photons seen by this pmt, let's get the propagation times
            if (is_visible){
                travel_time = distance_travelled / (c_cm_per_nsec / vis_index_of_refraction); // units of nsec
            }
            else{
                travel_time = distance_travelled / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec
            }

            // now all that's left to do is save!
            pmt_photon_yields.push_back(photons_in_this_pmt);
            pmt_photon_propagation_times.push_back(travel_time);
        }

    }

}

// I bet we can combine vertex_to_TPB_propagation and TPB_to_PMT_propagation functions into 1...
void vertex_to_TPB_to_PMT_propagation(double const & full_acceptance,
                            double const & c_cm_per_nsec,
                            double const & uv_index_of_refraction,
                            double const & UV_absorption_length,
                            double const & vertex_x,
                            double const & vertex_y,
                            double const & vertex_z,
                            double const & n_photons_produced,
                            std::vector<std::vector<double>> const & pmt_parsed_information_,
                            std::vector<std::vector<double>> const & locations_to_check_information_,
                            std::vector<std::vector<double>> const & locations_to_check_to_pmt_yield_,
                            std::vector<std::vector<double>> const & locations_to_check_to_pmt_travel_time_,
                            std::vector<std::vector<double>> & cumulative_pmt_photon_yields,
                            std::vector<std::vector<double>> & cumulative_pmt_photon_propagation_times){
    // define some things
    double omega;
    double distance_travelled;
    double photons_at_secondary_location;
    double travel_time;
    double loc_x;
    double loc_y;
    double loc_z;
    double facing_dir_x;
    double facing_dir_y;
    double facing_dir_z;
    double pmt_portion;
    double facing_area;
    double side_area;
    std::vector<double> loc_to_pmt_photon_yields;
    std::vector<double> loc_to_pmt_photon_propagation_times;

    // let's start by looping our our secondary locations parsed information
    for (size_t loc_it = 0; loc_it < locations_to_check_information_.size(); loc_it ++){

        loc_x = locations_to_check_information_.at(loc_it).at(0);
        loc_y = locations_to_check_information_.at(loc_it).at(1);
        loc_z = locations_to_check_information_.at(loc_it).at(2);
        facing_dir_x = locations_to_check_information_.at(loc_it).at(3);
        facing_dir_y = locations_to_check_information_.at(loc_it).at(4);
        facing_dir_z = locations_to_check_information_.at(loc_it).at(5);
        pmt_portion = locations_to_check_information_.at(loc_it).at(6); // portion_light_reflected_by_tpb_ if we are NOT on a pmt, otherwise 0.5 * portion_light_reflected_by_tpb_ 
        facing_area = locations_to_check_information_.at(loc_it).at(7);
        side_area = locations_to_check_information_.at(loc_it).at(8);

        // let's get solid angle and distance from vertex to this loc
        get_solid_angle_and_distance_vertex_to_location(vertex_x, vertex_y, vertex_z,
                                                        loc_x, loc_y, loc_z,
                                                        facing_dir_x, facing_dir_y, facing_dir_z,
                                                        facing_area, side_area, omega, distance_travelled);
        // now call function to get light yields
        // we are propagating scint light from vertex to TPB locs, so want to use UV absorption length
        get_light_yield(omega, full_acceptance, distance_travelled, n_photons_produced, UV_absorption_length, pmt_portion, photons_at_secondary_location);

        // now let's get the travel time as well
        travel_time = distance_travelled / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec

        // now we've already pre-computed the light yield for 1 photon going from every secondary loc to each pmt
        // as well as that photon's travel time
        loc_to_pmt_photon_yields.clear();
        loc_to_pmt_photon_propagation_times.clear();
        loc_to_pmt_photon_yields = locations_to_check_to_pmt_yield_.at(loc_it);
        loc_to_pmt_photon_propagation_times = locations_to_check_to_pmt_travel_time_.at(loc_it);

        // now we need to go through our vectors containing 1 element for every pmt and push_back to correct vector in our cumulative vector
        for (size_t j = 0; j < loc_to_pmt_photon_yields.size(); j++){
            // j is index of pmt
            cumulative_pmt_photon_yields.at(j).push_back(loc_to_pmt_photon_yields.at(j) * photons_at_secondary_location);
            cumulative_pmt_photon_propagation_times.at(j).push_back(loc_to_pmt_photon_propagation_times.at(j) + travel_time);
        }
    }

}

void manual_binning(std::vector<double> const & bin_centers,
                    double const & bin_width,
                    std::vector<double> const & times_to_bin,
                    std::vector<double> const & charges_to_bin,
                    std::vector<double> & this_pmt_binned_charges){
    // function to bin our charges
    size_t current_bin_idx;
    // we are going to loop over the times_to_bin and figure out which bin idx they correspond to
    for (size_t data_it = 0; data_it < charges_to_bin.size(); data_it ++){
        current_bin_idx = (size_t) times_to_bin.at(data_it) / bin_width;
        this_pmt_binned_charges.at(current_bin_idx) += charges_to_bin.at(data_it);
    }
}

void get_yields_per_pmt(std::vector<double> const & direct_photon_yields,
                        std::vector<double> const & direct_photon_propagation_times,
                        std::vector<std::vector<double>> const & indirect_photon_yields,
                        std::vector<std::vector<double>> const & indirect_photon_propagation_times,
                        std::vector<double> const & light_times,
                        std::vector<double> const & light_profile,
                        std::vector<double> const & bin_centers,
                        double const & bin_width,
                        double const & noise_rate_per_time_bin,
                        std::vector<std::vector<double>> & binned_charges) {

    // final function!!! just need to put everthing together
    // at the end of the day, we want to take the light seen in every pmt, bin the charges, and save to binned_charges (dims = n_pmts x n_time_bins)

    binned_charges = std::vector<std::vector<double>>(direct_photon_yields.size(), std::vector<double>(bin_centers.size()));

    double direct_yield;
    double direct_time_offset;
    std::vector<double> adjusted_light_yields (light_profile.size());
    std::vector<double> adjusted_light_times (light_times.size());
    std::vector<double> this_pmt_binned_charges(bin_centers.size());

    for (size_t pmt_it = 0; pmt_it < direct_photon_yields.size(); pmt_it ++) {
        std::vector<double> & this_pmt_binned_charges = binned_charges.at(pmt_it);

        direct_yield = direct_photon_yields.at(pmt_it);
        direct_time_offset = direct_photon_propagation_times.at(pmt_it);
        size_t first_direct_time_bin = (size_t) direct_time_offset / bin_width;
        size_t max_direct_time_bin = (size_t) (light_times.back() + direct_time_offset) / bin_width;
        max_direct_time_bin = std::min(max_direct_time_bin, bin_centers.size() - 1);
        max_direct_time_bin += 1;
        size_t num_direct_bins = max_direct_time_bin - first_direct_time_bin;

        for (size_t light_it = 0; light_it < num_direct_bins; ++light_it){
            size_t bin_it = first_direct_time_bin + light_it;
            this_pmt_binned_charges.at(bin_it) += direct_yield * light_profile.at(light_it);
        }

        // while we're on this pmt, let's see how much indirect charge there is
        std::vector<double> const & indirect_yields = indirect_photon_yields.at(pmt_it);
        std::vector<double> const & indirect_time_offsets = indirect_photon_propagation_times.at(pmt_it);

        for (size_t indirect_it = 0; indirect_it < indirect_yields.size(); ++indirect_it) {
            size_t first_indirect_time_bin = (size_t) indirect_time_offsets.at(indirect_it) / bin_width;
            size_t max_indirect_time_bin = (size_t) (light_times.back() + indirect_time_offsets.at(indirect_it)) / bin_width;
            max_indirect_time_bin = std::min(max_indirect_time_bin, bin_centers.size() - 1);
            max_indirect_time_bin += 1;
            size_t num_indirect_bins = max_indirect_time_bin - first_indirect_time_bin;

            double indirect_yield = indirect_yields.at(indirect_it);

            for (size_t light_it = 0; light_it < num_indirect_bins; ++light_it){
                size_t bin_it = first_indirect_time_bin + light_it;
                this_pmt_binned_charges.at(bin_it) += indirect_yield * light_profile.at(light_it);
            }
        }

        // ok so we have all the data we need for this pmt, time to bin!
        this_pmt_binned_charges.clear();
        this_pmt_binned_charges.resize(bin_centers.size());
        manual_binning(bin_centers, bin_width, adjusted_light_times, adjusted_light_yields, this_pmt_binned_charges);
        // now save!!!
        binned_charges.push_back(this_pmt_binned_charges);

    }
}

void put_simulation_steps_together(double const & full_acceptance,
                                   double const & c_cm_per_nsec,
                                   double const & uv_index_of_refraction,
                                   double const & vis_index_of_refraction,
                                   double const & quantum_efficiency,
                                   double const & n_photons_produced,
                                   double const & UV_absorption_length,
                                   double const & vis_absorption_length,
                                   std::vector<std::vector<double>> const & pmt_parsed_information_,
                                   std::vector<std::vector<double>> const & locations_to_check_information_,
                                   std::vector<std::vector<double>> const & locations_to_check_to_pmt_yield_,
                                   std::vector<std::vector<double>> const & locations_to_check_to_pmt_travel_time_,
                                   double const & vertex_x,
                                   double const & vertex_y,
                                   double const & vertex_z,
                                   std::vector<double> const & light_times,
                                   std::vector<double> const & light_profile,
                                   std::vector<double> const & bin_centers,
                                   double const & bin_width,
                                   double const & noise_rate_per_time_bin,
                                   std::vector<std::vector<double>> & binned_charges){

    // this function puts all of our photon simulation code together and returns the binned yields of charge in each pmt
    // some definitions
    bool is_visible;
    std::vector<double> direct_pmt_photon_yields;
    std::vector<double> direct_pmt_photon_propagation_times;
    std::vector<double> secondary_location_photon_yields;
    std::vector<double> secondary_location_photon_propagation_times;
    std::vector<std::vector<double>> cumulative_pmt_photon_yields(200);
    std::vector<std::vector<double>> cumulative_pmt_photon_propagation_times(200);

    // first we propagate light from vertex to pmts directly
    is_visible = false;
    vertex_to_pmt_propagation(full_acceptance, c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, quantum_efficiency,
                              n_photons_produced, UV_absorption_length, pmt_parsed_information_, is_visible, vertex_x, vertex_y, vertex_z,
                              direct_pmt_photon_yields, direct_pmt_photon_propagation_times);

    // now we propagate light from vertex to TPB locs
    vertex_to_TPB_to_PMT_propagation(full_acceptance, c_cm_per_nsec, uv_index_of_refraction, UV_absorption_length,
                                     vertex_x, vertex_y, vertex_z, n_photons_produced, pmt_parsed_information_,
                                     locations_to_check_information_, locations_to_check_to_pmt_yield_, locations_to_check_to_pmt_travel_time_,
                                     cumulative_pmt_photon_yields, cumulative_pmt_photon_propagation_times);

    // now we can put it all together by binning!
    get_yields_per_pmt(direct_pmt_photon_yields, direct_pmt_photon_propagation_times, cumulative_pmt_photon_yields, cumulative_pmt_photon_propagation_times,
                       light_times, light_profile, bin_centers, bin_width, binned_charges);
}

void FrameThread(std::atomic<bool> & running,
                 double const & full_acceptance,
                 double const & c_cm_per_nsec,
                 double const & uv_index_of_refraction,
                 double const & vis_index_of_refraction,
                 double const & quantum_efficiency,
                 double const & n_photons_produced,
                 double const & UV_absorption_length,
                 double const & vis_absorption_length,
                 std::vector<std::vector<double>> const & pmt_parsed_information_,
                 std::vector<std::vector<double>> const & locations_to_check_information_,
                 std::vector<std::vector<double>> const & locations_to_check_to_pmt_yield_,
                 std::vector<std::vector<double>> const & locations_to_check_to_pmt_travel_time_,
                 std::vector<double> * vertex,
                 std::vector<double> const & light_times,
                 std::vector<double> const & light_profile,
                 std::vector<double> const & bin_centers,
                 double const & bin_width,
                 std::vector<std::vector<double>> & binned_charges){

    // call simulation code
    put_simulation_steps_together(full_acceptance, c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, quantum_efficiency,
                                  n_photons_produced, UV_absorption_length, vis_absorption_length, pmt_parsed_information_, locations_to_check_information_,
                                  locations_to_check_to_pmt_yield_, locations_to_check_to_pmt_travel_time_,
                                  vertex->at(0), vertex->at(1), vertex->at(2), light_times, light_profile, bin_centers, bin_width, binned_charges);


    running.store(false);

}

void RunFrameThread(PhotonPropagationJob * job,
                    double const & full_acceptance,
                    double const & c_cm_per_nsec,
                    double const & uv_index_of_refraction,
                    double const & vis_index_of_refraction,
                    double const & quantum_efficiency,
                    double const & n_photons_produced,
                    double const & UV_absorption_length,
                    double const & vis_absorption_length,
                    std::vector<std::vector<double>> const & pmt_parsed_information_,
                    std::vector<std::vector<double>> const & locations_to_check_information_,
                    std::vector<std::vector<double>> const & locations_to_check_to_pmt_yield_,
                    std::vector<std::vector<double>> const & locations_to_check_to_pmt_travel_time_,
                    std::vector<double> const & light_times,
                    std::vector<double> const & light_profile,
                    std::vector<double> const & bin_centers,
                    double const & bin_width){

    job->running.store(true);
    job->thread = std::thread(FrameThread,
                              std::ref(job->running),
                              std::cref(full_acceptance),
                              std::cref(c_cm_per_nsec),
                              std::cref(uv_index_of_refraction),
                              std::cref(vis_index_of_refraction),
                              std::cref(quantum_efficiency),
                              std::cref(n_photons_produced),
                              std::cref(UV_absorption_length),
                              std::cref(vis_absorption_length),
                              std::cref(pmt_parsed_information_),
                              std::cref(locations_to_check_information_),
                              std::cref(locations_to_check_to_pmt_yield_),
                              std::cref(locations_to_check_to_pmt_travel_time_),
                              job->vertex,
                              std::cref(light_times),
                              std::cref(light_profile),
                              std::cref(bin_centers),
                              std::cref(bin_width),
                              std::ref(*job->binned_charges));
}

size_t findNearestIndex(std::vector<double> const & vec, double targetValue) {
 
    size_t nearestIndex = 0;
    double minDifference = std::abs(vec[0] - targetValue);

    for (size_t i = 1; i < vec.size(); ++i) {
        double difference = std::abs(vec[i] - targetValue);
        if (difference < minDifference) {
            minDifference = difference;
            nearestIndex = i;
        }
    }

    return nearestIndex;
}

double PhotonPropagation::GetSimulation(double const & singlet_ratio_,
                                        double const & triplet_ratio_,
                                        double const & singlet_tau_,
                                        double const & triplet_tau_,
                                        double const & recombination_tau_,
                                        double const & TPB_ratio_,
                                        double const & TPB_tau_,
                                        double const & UV_absorption_length_,
                                        double const & n_photons_produced_) {
    // will be used for multi-threading our simulation jobs
    std::deque<PhotonPropagationJob *> free_jobs;
    std::deque<PhotonPropagationJob *> running_jobs;
    std::deque<PhotonPropagationResult> results;
    size_t min_vertex_idx = 0;

    // first let's set up our light profile
    double start_time = 0.0;
    double end_time = 92.0;
    double bin_width = 2.0;
    size_t n_time_bins = (size_t)((end_time - start_time)/bin_width);

    std::vector<double> times(n_time_bins);
    std::vector<double> bin_centers(n_time_bins + 10);

    for (size_t i = 0; i < n_time_bins; i++){
        times.at(i) = (double) i * 2.0;
    }

    for (size_t j = 0; j < bin_centers.size(); j++){
        bin_centers.at(j) = (double) j * 2.0 + 1.0;
    }

    // now let's get out light profile
    std::vector<double> light_profile(n_time_bins);
    get_total_light_profile(singlet_ratio_,
                            triplet_ratio_,
                            singlet_tau_,
                            triplet_tau_,
                            recombination_tau_,
                            TPB_tau_,
                            times,
                            light_profile);


    // ok so we've seen our geometry file, pre-computed the lists of pmt info and locations to check for light propagation
    // we've also pre-compute the verticies of each event we want to simulate
    // and we just computed our light profile....all that's left to do is propagate our light!
    // we have the parameter total_events_that_escaped which says how many events we are simulating

    // so now we can loop over our vector containing vector of verticies to simulate
    // let's make out final vector that we will be saving to
    std::vector<std::vector<double>> summed_over_events_binned_charges(n_pmts_to_simulate, std::vector<double>(bin_centers.size(), 0.0));
                                                                            // dimensions are n_pmts x n_time_bins
    std::vector<std::vector<double>> summed_over_squared_events_binned_charges(n_pmts_to_simulate, std::vector<double>(bin_centers.size(), 0.0));
                                                                            // dimensions are n_pmts x n_time_bins
    // loop over all verticies to simulate
    for (size_t vert_it = 0; vert_it < verticies_to_simuate_.size(); vert_it++){
        // now ready to simulate.. this is where we want to dispatch our threads
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
			    if(not running_jobs[i]->running.load()) {
				    PhotonPropagationJob * job = running_jobs[i];
                    running_jobs.erase(running_jobs.begin() + i);
                    free_jobs.push_back(job);
                    job->thread.join();
                    results[job->vertex_index - min_vertex_idx].done = true;
                } else {
                    PhotonPropagationJob * job = running_jobs[i];
                }
            }

            // Check for any done results and push the corresponding frames
            size_t results_done = 0;
            for(size_t i=0; i<results.size(); ++i) {
                if(results[i].done) {
                    // let's save to all_events_binned_charges
                    for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
                        for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
                            summed_over_events_binned_charges[pmt_it][time_bin_it] += results[i].binned_charges->at(pmt_it)[time_bin_it];
                            summed_over_squared_events_binned_charges[pmt_it][time_bin_it] += std::pow(results[i].binned_charges->at(pmt_it)[time_bin_it], 2.0);
                        }
                    }
                    results[i].binned_charges = nullptr;
                    results[i].vertex_index = 0;
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
            PhotonPropagationJob * job = nullptr;

            if(free_jobs.size() > 0) {
                job = free_jobs.front();
                job->running.store(false);
                free_jobs.pop_front();
            } else if(running_jobs.size() < num_threads) {
                job = new PhotonPropagationJob();
                job->running.store(false);
                job->thread_index = running_jobs.size();
                job->binned_charges = new std::vector<std::vector<double>>();
            }

            if(job != nullptr and results.size() < max_cached_vertices) {
                job->running.store(true);
                running_jobs.push_back(job);
                job->vertex_index = vert_it;
                job->vertex = &verticies_to_simuate_[job->vertex_index];
                results.emplace_back();
                results.back().binned_charges = job->binned_charges;
                results.back().vertex_index = job->vertex_index;
                results.back().done = false;
                RunFrameThread(job, full_acceptance, c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, pmt_quantum_efficiency,
                               n_photons_produced_, UV_absorption_length_, visible_absorption_length_, pmt_parsed_information_, locations_to_check_information_,
                               locations_to_check_to_pmt_yield_, locations_to_check_to_pmt_travel_time_, times, light_profile,bin_centers, bin_width);
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
            if(not running_jobs[i]->running.load()) {
                PhotonPropagationJob * job = running_jobs[i];
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                job->thread.join();
                results[job->vertex_index - min_vertex_idx].done = true;
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results[i].done) {
                // let's save to all_events_binned_charges
                for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
                    for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
                        summed_over_events_binned_charges[pmt_it][time_bin_it] += results[i].binned_charges->at(pmt_it)[time_bin_it];
                        summed_over_squared_events_binned_charges[pmt_it][time_bin_it] += std::pow(results[i].binned_charges->at(pmt_it)[time_bin_it], 2.0);
                    }
                }
                results[i].binned_charges = nullptr;
                results[i].vertex_index = 0;
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

    // instead of just returning the averaged_all_events_binned_charges, we are going to compute the liklihood quickly
    // we are computing the liklihood on a per-pmt basis on a per-time bin basis
    // we've already found the max time of the data and set that equal to 0
    // let's do the same with the simulation
    // then we actually want to compute the liklihood evaluation from 5 nsec before the max time until 70 nsec after the max time

    double simulation_time_of_max;
    double simulation_max_value;
    double total_charge_per_time_bin;
    for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
        total_charge_per_time_bin = 0;
        for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
            total_charge_per_time_bin += summed_over_events_binned_charges[pmt_it][time_bin_it];
        }
        if (total_charge_per_time_bin > simulation_max_value){
            simulation_max_value = total_charge_per_time_bin;
            simulation_time_of_max = bin_centers[time_bin_it];
        }
    }

    // let's subtract off the time of the max val
    std::vector<double> simulation_times;
    for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
        simulation_times.push_back(bin_centers[time_bin_it] - simulation_time_of_max);
    }

    double min_time_to_evaluate = -5.0;  // relative to summed wf peak times
    double max_time_to_evaluate = 70;

    double total_nllh;

    // now we can loop over our simulation, check that we are within the time bands to calculate the nllh, then calculate it!
    double mu;
    double sigma_squared;
    double k;
    size_t this_time_bin_in_data_idx;

    for (size_t time_bin_it = 0; time_bin_it < simulation_times.size(); time_bin_it ++){
        if (simulation_times[time_bin_it] >= min_time_to_evaluate and simulation_times[time_bin_it] <= max_time_to_evaluate){
            // ok we are within the times we want to evaluate!
            // we need to find the corresponding idx in times_of_data_points_
            this_time_bin_in_data_idx = findNearestIndex(times_of_data_points_, simulation_times[time_bin_it]);        
            for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
                // so now we are looping over each pmt, we need to calculate the nllh for each time bin for each pmt
                k = data_series_[pmt_it][this_time_bin_in_data_idx];
                mu = summed_over_events_binned_charges[pmt_it][time_bin_it] * (n_data_samples_ / total_events_that_escaped);
                sigma_squared = summed_over_squared_events_binned_charges[pmt_it][time_bin_it] * std::pow((n_data_samples_ / total_events_that_escaped), 2);
                total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
            }
        }
    }


    return total_nllh;

}
