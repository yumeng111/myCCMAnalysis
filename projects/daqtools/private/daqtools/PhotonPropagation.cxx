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
#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <math.h>

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
#include <daqtools/PhotonPropagation.h>

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

void combined_area_penalty(double const & chunk_center_x,
                           double const & chunk_center_y,
                           double const & x0,
                           double const & x1,
                           double const & y0,
                           double const & y1,
                           bool const & top_bottom_flag,
                           I3Vector<I3Vector<double>> const & relevant_coated_pmt_locs,
                           I3Vector<I3Vector<double>> const & relevant_uncoated_pmt_locs,
                           double const & detector_radius,
                           double const & pmt_radius,
                           double const & half_chunk_hypotenuse,
                           I3Vector<I3Vector<double>> & this_chunk_valid_points, 
                           I3Vector<I3Vector<double>> & this_chunk_invalid_points,
                           I3Vector<I3Vector<double>> & this_chunk_coated_pmt_points){

    // x0, x1, y0, and y1 are the corners of the chunk we are intersting in calcualting overlap with
    // top_bottom_flag is true if point is on top/bottom and false if on side --> dictates if we check if point is within the detector
    // before calling this function, we figure out the relevant coated and uncoated pmt locs and save to vectors

    // granularity for chunking up our chunk
    size_t n_chunks_x_side = 100;
    size_t n_chunks_y_side = 100;

    // we are splitting up this chunk such that we do not have any points on the edges of the chunk bc that will cause over counting
    // first chunking up in x
    std::vector<double> possible_chunk_x_positions(n_chunks_x_side);
    for (size_t i = 0; i < n_chunks_x_side; i++){
        double x = x0 + ((x1 - x0) * (double)i / (double)n_chunks_x_side);
        possible_chunk_x_positions[i] = x;
    }
    double width = std::abs(possible_chunk_x_positions[1] - possible_chunk_x_positions[0]);
    for (size_t i = 0; i < n_chunks_x_side; i++){
        possible_chunk_x_positions[i] += (width / 2.0);
    }

    // now chunking up in y
    std::vector<double> possible_chunk_y_positions(n_chunks_y_side);
    for (size_t i = 0; i < n_chunks_y_side; i++){
        double y = y1 + ((y0 - y1) * (double)i / (double)n_chunks_y_side);
        possible_chunk_y_positions[i] = y;
    }
    double height = std::abs(possible_chunk_y_positions[1] - possible_chunk_y_positions[0]);
    for (size_t i = 0; i < n_chunks_y_side; i++){
        possible_chunk_y_positions[i] += (height / 2.0);
    }

    // now that we are done with that, we can loop over our points in the chunk and check:
    // 1) is this point actually inside the chunk we are considering? should always be yes...but good to check
    // 2) is this point inside the detector dimenions? if no -- this is an invalid point!
    // 3) is this point on top of an uncoated pmt? if yes -- this is an invalid point!
    // 4) is this point on top of a coated pmt? if yes -- this is a valid point, but keep track for tpb portions

    I3Vector<double> this_point;
    for (size_t x_chunk_it = 0; x_chunk_it < possible_chunk_x_positions.size(); x_chunk_it ++){
        for (size_t y_chunk_it = 0; y_chunk_it < possible_chunk_y_positions.size(); y_chunk_it ++){
            // booleans to keep track of this point
            bool valid_point = true;
            bool on_coated_pmt = false;

            double this_point_in_chunk_x = possible_chunk_x_positions.at(x_chunk_it);
            double this_point_in_chunk_y = possible_chunk_y_positions.at(y_chunk_it);

            // first things first, let's make sure this point is actually on our chunk by checking radius
            double this_point_dist_from_center_chunk = std::sqrt(std::pow(this_point_in_chunk_x - chunk_center_x, 2) + std::pow(this_point_in_chunk_y - chunk_center_y, 2)); 
            if (this_point_dist_from_center_chunk > half_chunk_hypotenuse){
                // oops! this point is not actually on our chunk
                continue;
            }

            // now we check if this point is inside the detector
            // we start off assuming all points are valid, so want to check if invalid and flip bool
            if (top_bottom_flag){
                double this_point_dist_to_edge_of_detector = std::sqrt(std::pow(this_point_in_chunk_x, 2) + std::pow(this_point_in_chunk_y, 2));
                if (this_point_dist_to_edge_of_detector > detector_radius){
                    // oops! this point is outside the radius of the detector! definitely invalid
                    valid_point = false;
                }
            }

            // now if this point is within the detector, we need to check if it's on an uncoated pmt --> that will also flip valid bool
            if (valid_point){
                // we we have a list of the un coated pmt x and y locs, let's loop through and calculate distances from this point to all of those
                for(size_t uncoated_pmt_it = 0; uncoated_pmt_it < relevant_uncoated_pmt_locs.size(); uncoated_pmt_it++){
                    double this_uncoated_pmt_x = relevant_uncoated_pmt_locs[uncoated_pmt_it][0];
                    double this_uncoated_pmt_y = relevant_uncoated_pmt_locs[uncoated_pmt_it][1];
                    if (!top_bottom_flag and this_uncoated_pmt_x < -319.0){
                        // these are pmts on the edges of where we start counting circumference
                        // let's translate these c, z positions to the edge edge, check, and then check like normal
                        double adjusted_pmt_x = -this_uncoated_pmt_x;
                        double adjusted_dist = std::sqrt(std::pow(this_point_in_chunk_x - adjusted_pmt_x, 2) + std::pow(this_point_in_chunk_y - this_uncoated_pmt_y, 2));
                        if (adjusted_dist <= pmt_radius){
                            valid_point = false;
                        }
                    }
                    double distance_from_this_point_to_uncoated_pmt = std::sqrt(std::pow(this_point_in_chunk_x - this_uncoated_pmt_x, 2) + std::pow(this_point_in_chunk_y - this_uncoated_pmt_y, 2));
                    if (distance_from_this_point_to_uncoated_pmt <= pmt_radius){
                        // we are on an uncoated pmt!!!
                        valid_point = false;
                    }
                }
            }

            // ok so finished check if on uncoated pmt -- last thing to check is if we are on a coated pmt (only if this is still a valid point of course)!
            if (valid_point){
                // we we have a list of the coated pmt x and y locs, let's loop through and calculate distances from this point to all of those
                for(size_t coated_pmt_it = 0; coated_pmt_it < relevant_coated_pmt_locs.size(); coated_pmt_it++){
                    double this_coated_pmt_x = relevant_coated_pmt_locs[coated_pmt_it][0];
                    double this_coated_pmt_y = relevant_coated_pmt_locs[coated_pmt_it][1];
                    if (!top_bottom_flag and this_coated_pmt_x < -319.0){
                        // these are pmts on the edges of where we start counting circumference
                        // let's translate these c, z positions to the edge edge, check, and then check like normal
                        double adjusted_pmt_x = -this_coated_pmt_x;
                        double adjusted_dist = std::sqrt(std::pow(this_point_in_chunk_x - adjusted_pmt_x, 2) + std::pow(this_point_in_chunk_y - this_coated_pmt_y, 2));
                        if (adjusted_dist <= pmt_radius){
                            on_coated_pmt = true;
                        }
                    }
                    double distance_from_this_point_to_coated_pmt = std::sqrt(std::pow(this_point_in_chunk_x - this_coated_pmt_x, 2) + std::pow(this_point_in_chunk_y - this_coated_pmt_y, 2));
                    if (distance_from_this_point_to_coated_pmt <= pmt_radius){
                        // we are on an coated pmt!!!
                        on_coated_pmt = true;
                    }
                }
            }

            // now we can save this point as appropriate!
            this_point.clear();
            this_point.push_back(this_point_in_chunk_x);
            this_point.push_back(this_point_in_chunk_y);

            if (valid_point){
                this_chunk_valid_points.push_back(this_point);
                if (on_coated_pmt){
                    this_chunk_coated_pmt_points.push_back(this_point);
                }
            }
            if (!valid_point){
                this_chunk_invalid_points.push_back(this_point);
            }

        }
    }


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

void get_511kev_photon(double const & source_diameter,
                       double const & source_rod_lower_end_cap,
                       double const & pos_rad,
                       double const & decay_constant,
                       double const & detector_radius,
                       double const & detector_lower_end_cap,
                       std::mt19937 & gen,
                       std::uniform_real_distribution<double> & dis_angle,
                       std::uniform_real_distribution<double> & dis_0_1,
                       std::uniform_real_distribution<double> & dis_neg_1_1,
                       bool & escaped_photon_1,
                       double & final_x_photon_1,
                       double & final_y_photon_1,
                       double & final_z_photon_1,
                       bool & escaped_photon_2,
                       double & final_x_photon_2,
                       double & final_y_photon_2,
                       double & final_z_photon_2){
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
    check_escape(x, y, z, nx, ny, nz, source_rod_lower_end_cap, detector_radius, decay_constant, pos_rad,
                detector_lower_end_cap, gen, dis_0_1, escaped_photon_1, final_x_photon_1, final_y_photon_1, final_z_photon_1);
    check_escape(x, y, z, -nx, -ny, -nz, source_rod_lower_end_cap, detector_radius, decay_constant, pos_rad,
                detector_lower_end_cap, gen, dis_0_1, escaped_photon_2, final_x_photon_2, final_y_photon_2, final_z_photon_2);
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

        pmt_x_loc = pmt_parsed_information_.at(pmt_it).at(0);
        pmt_y_loc = pmt_parsed_information_.at(pmt_it).at(1);
        pmt_z_loc = pmt_parsed_information_.at(pmt_it).at(2);
        facing_dir_x = pmt_parsed_information_.at(pmt_it).at(3);
        facing_dir_y = pmt_parsed_information_.at(pmt_it).at(4);
        facing_dir_z = pmt_parsed_information_.at(pmt_it).at(5);
        coating_flag = pmt_parsed_information_.at(pmt_it).at(6); // coating flag == 1.0 for coated pmts and 0.0 for uncoated!
        pmt_facing_area = pmt_parsed_information_.at(pmt_it).at(7);
        pmt_side_area = pmt_parsed_information_.at(pmt_it).at(8);

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
PhotonPropagation::PhotonPropagation() {}

void PhotonPropagation::SetData(I3Vector<I3Vector<double>> data_series){
    data_series_ = data_series;
    // data_series_ is index by n_pmts x n_time_bins
    // while we are setting the data, let's also set the time
    // we want the peak of the data time to be at 0 and 2 nsec binning

    times_of_data_points_.clear();
    for (size_t second_dim_it = 0; second_dim_it < data_series_.at(0).size(); second_dim_it ++ ){
        times_of_data_points_.push_back((double)second_dim_it * 2.0);
    }

    max_data_value_ = 0;
    time_of_max_data_value_ = 0;

    // let's first find the time of max, then set that equal to zero
    double total_charge_per_time_bin;
    for (size_t time_bin_it = 0; time_bin_it < times_of_data_points_.size(); time_bin_it ++){
        total_charge_per_time_bin = 0;
        for (size_t pmt_it = 0; pmt_it < data_series_.size(); pmt_it ++){
            total_charge_per_time_bin += data_series_.at(pmt_it).at(time_bin_it);
        }
        if (total_charge_per_time_bin > max_data_value_){
            max_data_value_ = total_charge_per_time_bin;
            time_of_max_data_value_ = times_of_data_points_.at(time_bin_it);
        }
    }

    // ok now we can subtract time_of_max_data_value_ from times_of_data_points_s
    for (size_t time_bin_it = 0; time_bin_it < times_of_data_points_.size(); time_bin_it ++){
        times_of_data_points_.at(time_bin_it) -= time_of_max_data_value_;
    }

    /*// print out data to check
    std::cout << "printing data to check : " << std::endl;
    total_charge_per_time_bin = 0;
    for (size_t time_bin_it = 0; time_bin_it < times_of_data_points_.size(); time_bin_it ++){
        total_charge_per_time_bin = 0;
        for (size_t pmt_it = 0; pmt_it < data_series_.size(); pmt_it ++){
            total_charge_per_time_bin += data_series_.at(pmt_it).at(time_bin_it);
        }
        std::cout << "at time = " << times_of_data_points_.at(time_bin_it) << " summed charge = " << total_charge_per_time_bin << std::endl;
    }*/

}

void PhotonPropagation::SetDataSampleSize(size_t n_data_samples){
    n_data_samples_ = n_data_samples;
}

void PhotonPropagation::SetNThreads(size_t const & n_threads){
    if (n_threads == 0){
        num_threads = std::thread::hardware_concurrency();
    }
    else{
        num_threads = n_threads;
    }
}

size_t PhotonPropagation::GetNFaceChunks(){
    return face_chunks_counter;
}

size_t PhotonPropagation::GetNSideChunks(){
    return side_chunks_counter;
}

size_t PhotonPropagation::GetNSimulatedEvents(){
    return equivalent_events_in_data;
}

void PhotonPropagation::SetZOffset(double source_z_offset){
    source_z_offset_ = source_z_offset;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopUncoatedPMTLocs(){
    return uncoated_pmt_locs_top_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopCoatedPMTLocs(){
    return coated_pmt_locs_top_;
}
I3Vector<I3Vector<double>> PhotonPropagation::GetBottomUncoatedPMTLocs(){
    return uncoated_pmt_locs_bottom_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetBottomCoatedPMTLocs(){
    return coated_pmt_locs_bottom_;
}
I3Vector<I3Vector<double>> PhotonPropagation::GetSideUncoatedPMTLocs(){
    return uncoated_pmt_locs_side_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetSideCoatedPMTLocs(){
    return coated_pmt_locs_side_;
}

double PhotonPropagation::GetTopLocWidth(){
    return final_top_loc_width_;
}

double PhotonPropagation::GetTopLocHeight(){
    return final_top_loc_height_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopLocXY(){
    return top_loc_xy_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopLocInsideCoatedPMTXY(){
    return top_loc_inside_uncoated_pmt_xy_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopLocOutsideCoatedPMTXY(){
    return top_loc_outside_uncoated_pmt_xy_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopLocInsideDetectorRadiusXY(){
    return top_loc_inside_detector_radius_xy_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetTopLocOutsideDetectorRadiusXY(){
    return top_loc_outside_detector_radius_xy_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllTopLocValidDots(){
    return all_top_loc_valid_dots_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllTopLocInvalidDots(){
    return all_top_loc_invalid_dots_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllTopLocCoatedPMTDots(){
    return all_top_loc_coated_pmt_dots_;
}
I3Vector<I3Vector<double>> PhotonPropagation::GetAllBottomLocValidDots(){
    return all_bottom_loc_valid_dots_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllBottomLocInvalidDots(){
    return all_bottom_loc_invalid_dots_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllBottomLocCoatedPMTDots(){
    return all_bottom_loc_coated_pmt_dots_;
}
I3Vector<I3Vector<double>> PhotonPropagation::GetAllSideLocValidDots(){
    return all_side_loc_valid_dots_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllSideLocInvalidDots(){
    return all_side_loc_invalid_dots_;
}

I3Vector<I3Vector<double>> PhotonPropagation::GetAllSideLocCoatedPMTDots(){
    return all_side_loc_coated_pmt_dots_;
}

void PhotonPropagation::GetEventVertices(size_t const & n_events_to_simulate){
    // set our events parameter
    n_events_to_simulate_ = n_events_to_simulate;
    total_events_that_escaped = 0;
    equivalent_events_in_data = 0;

    // now let's empty out some vectors
    if (verticies_to_simuate_1275_.size() > 0){
        for (size_t i = 0; i < verticies_to_simuate_1275_.size(); i++){
            verticies_to_simuate_1275_[i].clear();
        }
        verticies_to_simuate_1275_.clear();
    }
    if (verticies_to_simuate_511_.size() > 0){
        for (size_t i = 0; i < verticies_to_simuate_511_.size(); i++){
            verticies_to_simuate_511_[i].clear();
        }
        verticies_to_simuate_511_.clear();
    }
    if (thread_1275_verticies_.size() > 0){
        for(size_t i = 0; i < thread_1275_verticies_.size(); i++){
            for (size_t j = 0; j < thread_1275_verticies_[i].size(); j++){
                thread_1275_verticies_[i][j].clear();
            }
            thread_1275_verticies_[i].clear();
        }
        thread_1275_verticies_.clear();
    }
    if (thread_511_verticies_.size() > 0){
        for(size_t i = 0; i < thread_511_verticies_.size(); i++){
            for (size_t j = 0; j < thread_511_verticies_[i].size(); j++){
                thread_511_verticies_[i][j].clear();
            }
            thread_511_verticies_[i].clear();
        }
        thread_511_verticies_.clear();
    }
    if (thread_verticies_.size() > 0){
        for(size_t i = 0; i < thread_verticies_.size(); i++){
            for (size_t j = 0; j < thread_verticies_[i].size(); j++){
                thread_verticies_[i][j].clear();
            }
            thread_verticies_[i].clear();
        }
        thread_verticies_.clear();
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
    bool escaped_photon_1;
    double final_x_photon_1;
    double final_y_photon_1;
    double final_z_photon_1;
    bool escaped_photon_2;
    double final_x_photon_2;
    double final_y_photon_2;
    double final_z_photon_2;
    std::vector<double> this_vertex (3);

    while (equivalent_events_in_data < n_events_to_simulate_){
        // so for the high energy bump, we want 3 photons -- isotropic 1275 kev, and 2 back to back 511 kev photons
        // let's call our functions to check if 511 and 1275 kev photons escape
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
        get_511kev_photon(source_diameter,
                          source_rod_lower_end_cap,
                          pos_rad,
                          decay_constant,
                          detector_radius,
                          detector_lower_end_cap,
                          gen,
                          dis_angle,
                          dis_0_1,
                          dis_neg_1_1,
                          escaped_photon_1,
                          final_x_photon_1,
                          final_y_photon_1,
                          final_z_photon_1,
                          escaped_photon_2,
                          final_x_photon_2,
                          final_y_photon_2,
                          final_z_photon_2);
        // now we only want events where all 3 photons escape
        if (escaped and escaped_photon_1 and escaped_photon_2){
            total_events_that_escaped += 3;
            equivalent_events_in_data += 1;
            // let's save this vertex!
            //std::cout << "[" << final_x << ", " << final_y << ", " << final_z  << "]," << std::endl;
            //std::cout << "[" << final_x_photon_1 << ", " << final_y_photon_1 << ", " << final_z_photon_1  << "]," << std::endl;
            //std::cout << "[" << final_x_photon_2 << ", " << final_y_photon_2 << ", " << final_z_photon_2  << "]," << std::endl;
            // first we need to apply our z offset! if not set, the default value is zero
            final_z += source_z_offset_;
            final_z_photon_1 += source_z_offset_;
            final_z_photon_2 += source_z_offset_;
            // now we can save
            this_vertex.clear();
            this_vertex.push_back(final_x);
            this_vertex.push_back(final_y);
            this_vertex.push_back(final_z);
            verticies_to_simuate_1275_.push_back(this_vertex);
            this_vertex.clear();
            this_vertex.push_back(final_x_photon_1);
            this_vertex.push_back(final_y_photon_1);
            this_vertex.push_back(final_z_photon_1);
            verticies_to_simuate_511_.push_back(this_vertex);
            this_vertex.clear();
            this_vertex.push_back(final_x_photon_2);
            this_vertex.push_back(final_y_photon_2);
            this_vertex.push_back(final_z_photon_2);
            verticies_to_simuate_511_.push_back(this_vertex);
        }
    }
    // we should also pre-emptively chunk out our vertices for the threads
    // so we have num_threads to work with and total_events_that_escaped to split over them (being careful to keep track of 1275 and 511 separtely)
    // we also know that we have twice as many 511 events as 1275, so can split threads into 511 and 1275
    if (num_threads == 1){
        std::vector<std::vector<double>> filler_vec;
        for (size_t vertices_1275_it = 0; vertices_1275_it < verticies_to_simuate_1275_.size(); vertices_1275_it ++){
            filler_vec.push_back(verticies_to_simuate_1275_.at(vertices_1275_it));
        }
        for (size_t vertices_511_it = 0; vertices_511_it < verticies_to_simuate_511_.size(); vertices_511_it ++){
            filler_vec.push_back(verticies_to_simuate_511_.at(vertices_511_it));
        }
        thread_verticies_.push_back(filler_vec);
    }
    else {
        size_t total_number_1275_verticies = verticies_to_simuate_1275_.size();
        size_t n_threads_for_1275 = (size_t) num_threads/3;
        size_t n_1275_verticies_per_thread_rounding_up = (size_t) ((total_number_1275_verticies + n_threads_for_1275 - 1) / n_threads_for_1275); // this is rounding up!
        size_t n_1275_verticies_per_thread_rounding_down = n_1275_verticies_per_thread_rounding_up - 1;
        size_t n_1275_threads_rounding_up = (size_t) (total_number_1275_verticies - (n_1275_verticies_per_thread_rounding_down * n_threads_for_1275)) /
                                                     (n_1275_verticies_per_thread_rounding_up - n_1275_verticies_per_thread_rounding_down);
        size_t n_1275_threads_rounding_down = n_threads_for_1275 - n_1275_threads_rounding_up;

        size_t total_number_511_verticies = verticies_to_simuate_511_.size();
        size_t n_threads_for_511 = num_threads - n_threads_for_1275;
        size_t n_511_verticies_per_thread_rounding_up = (size_t) ((total_number_511_verticies  + n_threads_for_511 - 1) / n_threads_for_511); // this is rounding up!
        size_t n_511_verticies_per_thread_rounding_down = n_511_verticies_per_thread_rounding_up - 1;
        size_t n_511_threads_rounding_up = (size_t) (total_number_511_verticies - (n_511_verticies_per_thread_rounding_down * n_threads_for_511)) /
                                                     (n_511_verticies_per_thread_rounding_up - n_511_verticies_per_thread_rounding_down);
        size_t n_511_threads_rounding_down = n_threads_for_511 - n_511_threads_rounding_up;

        // ok so now we can fill in vector of vector of vertex vectors for each thread type
        std::vector<std::vector<double>> verticies_for_one_thread;
        size_t vertices_in_this_clump = 0;
        bool finished_accumulating_threads_rounding_up = false;
        bool finished_accumulating_threads_rounding_down = false;
        size_t total_threads_rounding_up = 0;
        size_t total_threads_rounding_down = 0;

        for (size_t vertices_1275_it = 0; vertices_1275_it < verticies_to_simuate_1275_.size(); vertices_1275_it ++){
            if (total_threads_rounding_up == n_1275_threads_rounding_up){
                finished_accumulating_threads_rounding_up = true;
            }
            if (total_threads_rounding_down == n_1275_threads_rounding_down){
                finished_accumulating_threads_rounding_down = true;
            }

            if (finished_accumulating_threads_rounding_up == false){
                if (vertices_in_this_clump == (n_1275_verticies_per_thread_rounding_up - 1)){
                    verticies_for_one_thread.push_back(verticies_to_simuate_1275_.at(vertices_1275_it));
                    thread_1275_verticies_.push_back(verticies_for_one_thread);
                    vertices_in_this_clump = 0;
                    verticies_for_one_thread.clear();
                    total_threads_rounding_up += 1;
                }
                else {
                    verticies_for_one_thread.push_back(verticies_to_simuate_1275_.at(vertices_1275_it));
                    vertices_in_this_clump += 1;
                }
            }

            if (finished_accumulating_threads_rounding_up and finished_accumulating_threads_rounding_down == false){
                if (vertices_in_this_clump == (n_1275_verticies_per_thread_rounding_down - 1)){
                    verticies_for_one_thread.push_back(verticies_to_simuate_1275_.at(vertices_1275_it));
                    thread_1275_verticies_.push_back(verticies_for_one_thread);
                    vertices_in_this_clump = 0;
                    verticies_for_one_thread.clear();
                    total_threads_rounding_down += 1;
                }
                else {
                    verticies_for_one_thread.push_back(verticies_to_simuate_1275_.at(vertices_1275_it));
                    vertices_in_this_clump += 1;
                }
            }
        }

        verticies_for_one_thread.clear();
        vertices_in_this_clump = 0;
        finished_accumulating_threads_rounding_up = false;
        finished_accumulating_threads_rounding_down = false;
        total_threads_rounding_up = 0;
        total_threads_rounding_down = 0;

        for (size_t vertices_511_it = 0; vertices_511_it < verticies_to_simuate_511_.size(); vertices_511_it ++){
            if (total_threads_rounding_up == n_511_threads_rounding_up){
                finished_accumulating_threads_rounding_up = true;
            }
            if (total_threads_rounding_down == n_511_threads_rounding_down){
                finished_accumulating_threads_rounding_down = true;
            }

            if (finished_accumulating_threads_rounding_up == false){
                if (vertices_in_this_clump == (n_511_verticies_per_thread_rounding_up - 1)){
                    verticies_for_one_thread.push_back(verticies_to_simuate_511_.at(vertices_511_it));
                    thread_511_verticies_.push_back(verticies_for_one_thread);
                    vertices_in_this_clump = 0;
                    verticies_for_one_thread.clear();
                    total_threads_rounding_up += 1;
                }
                else {
                    verticies_for_one_thread.push_back(verticies_to_simuate_511_.at(vertices_511_it));
                    vertices_in_this_clump += 1;
                }
            }

            if (finished_accumulating_threads_rounding_up and finished_accumulating_threads_rounding_down == false){
                if (vertices_in_this_clump == (n_511_verticies_per_thread_rounding_down - 1)){
                    verticies_for_one_thread.push_back(verticies_to_simuate_511_.at(vertices_511_it));
                    thread_511_verticies_.push_back(verticies_for_one_thread);
                    vertices_in_this_clump = 0;
                    verticies_for_one_thread.clear();
                    total_threads_rounding_down += 1;
                }
                else {
                    verticies_for_one_thread.push_back(verticies_to_simuate_511_.at(vertices_511_it));
                    vertices_in_this_clump += 1;
                }
            }
        }

        // yayy! thread_511_verticies_ and thread_1275_verticies_ are all filled in correctly!
        // let's combine into one list to give threads

        for (size_t vert_1275_it = 0; vert_1275_it < thread_1275_verticies_.size(); vert_1275_it ++){
           thread_verticies_.push_back(thread_1275_verticies_.at(vert_1275_it));
        }

        for (size_t vert_511_it = 0; vert_511_it < thread_511_verticies_.size(); vert_511_it ++){
           thread_verticies_.push_back(thread_511_verticies_.at(vert_511_it));
        }

        // last check, go through thread_verticies_ and count how many events there
        //size_t event_counter_double_check = 0;
        //size_t thread_counter_double_check = 0;
        //for (size_t thread_it = 0; thread_it < thread_verticies_.size(); thread_it++){
        //    thread_counter_double_check += 1;
        //    for (size_t vert_it = 0; vert_it < thread_verticies_[thread_it].size(); vert_it++){
        //        event_counter_double_check += 1;
        //    }
        //}
        //std::cout << "thread_counter_double_check = " << thread_counter_double_check << std::endl;
        //std::cout << "event_counter_double_check = " << event_counter_double_check << std::endl;
    }
}

void PhotonPropagation::GetPMTInformation(I3FramePtr frame){
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    I3Map<CCMPMTKey, CCMOMGeo> const & pmt_geo = geo.pmt_geo;
    std::vector<double> this_pmt_info (9);

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

    // let's save some pmt locs for debugging
    I3Vector<double> this_pmt_loc;
    for (size_t pmt_it = 0; pmt_it < pmt_parsed_information_.size(); pmt_it ++){
        std::vector<double> this_pmt_info = pmt_parsed_information_.at(pmt_it);
        double pmt_coating_flag = this_pmt_info.at(6);

        double pmt_x = this_pmt_info.at(0);
        double pmt_y = this_pmt_info.at(1);
        double pmt_z = this_pmt_info.at(2);

        this_pmt_loc.clear();
        this_pmt_loc.push_back(pmt_x);
        this_pmt_loc.push_back(pmt_y);

        if (pmt_z == 58.0){
            if (pmt_coating_flag == 0.0){
                uncoated_pmt_locs_top_.push_back(this_pmt_loc);
            }
            if (pmt_coating_flag == 1.0){
                coated_pmt_locs_top_.push_back(this_pmt_loc);
            }
        }
        if (pmt_z == -58.0){
            if (pmt_coating_flag == 0.0){
                uncoated_pmt_locs_bottom_.push_back(this_pmt_loc);
            }
            if (pmt_coating_flag == 1.0){
                coated_pmt_locs_bottom_.push_back(this_pmt_loc);
            }
        }
        if (pmt_z > -58.0 and pmt_z < 58.0){
            // for pmts on the sides of the detector, we actually want to save pmt_c and pmt_z
            double pmt_theta = -atan2(pmt_y, pmt_x);
            pmt_theta = fmod(M_PI + pmt_theta, 2.0 * M_PI) - M_PI;
            double pmt_c = cylinder_radius * pmt_theta;

            this_pmt_loc.clear();
            this_pmt_loc.push_back(pmt_c);
            this_pmt_loc.push_back(pmt_z);

            if (pmt_coating_flag == 0.0){
                uncoated_pmt_locs_side_.push_back(this_pmt_loc);
            }
            if (pmt_coating_flag == 1.0){
                coated_pmt_locs_side_.push_back(this_pmt_loc);
            }
        }
    }


}


void PhotonPropagation::GetSecondaryLocs(double const & desired_chunk_width, double const & desired_chunk_height) {
    // set our chunk vars
    desired_chunk_width_ = desired_chunk_width;
    desired_chunk_height_ = desired_chunk_height;

    // reset the counters
    face_chunks_counter = 0;
    side_chunks_counter = 0;

    double area_top = 0;
    double area_sides = 0;
    double area_bottom = 0;
    double area_top_valid = 0;
    double area_sides_valid = 0;
    double area_bottom_valid = 0;
    // finally, let's reset some vectors that we're saving info to
    if (locations_to_check_information_.size() > 0){
        // our secondary loc vectors are already filled! let's empty!
        for (size_t i = 0; i < locations_to_check_information_.size(); i++){
            locations_to_check_information_[i].clear();
        }
        locations_to_check_information_.clear();
    }

    if (locations_to_check_to_pmt_yield_.size() > 0){
       for (size_t i = 0; i < locations_to_check_to_pmt_yield_.size(); i++){
            locations_to_check_to_pmt_yield_[i].clear();
        }
        locations_to_check_to_pmt_yield_.clear();
    }

    if (locations_to_check_to_pmt_travel_time_.size() > 0){
       for (size_t i = 0; i < locations_to_check_to_pmt_travel_time_.size(); i++){
            locations_to_check_to_pmt_travel_time_[i].clear();
        }
        locations_to_check_to_pmt_travel_time_.clear();
    }
    if (all_top_loc_valid_dots_.size() > 0){
        for (size_t i = 0; i < all_top_loc_valid_dots_.size(); i++){
            all_top_loc_valid_dots_[i].clear();
        }
        all_top_loc_valid_dots_.clear();
    }
    if (all_top_loc_invalid_dots_.size() > 0){
        for (size_t i = 0; i < all_top_loc_invalid_dots_.size(); i++){
            all_top_loc_invalid_dots_[i].clear();
        }
        all_top_loc_invalid_dots_.clear();
    }
    if (all_top_loc_coated_pmt_dots_.size() > 0){
        for (size_t i = 0; i < all_top_loc_coated_pmt_dots_.size(); i++){
            all_top_loc_coated_pmt_dots_[i].clear();
        }
        all_top_loc_coated_pmt_dots_.clear();
    }
    if (all_bottom_loc_valid_dots_.size() > 0){
        for (size_t i = 0; i < all_bottom_loc_valid_dots_.size(); i++){
            all_bottom_loc_valid_dots_[i].clear();
        }
        all_bottom_loc_valid_dots_.clear();
    }
    if (all_bottom_loc_invalid_dots_.size() > 0){
        for (size_t i = 0; i < all_bottom_loc_invalid_dots_.size(); i++){
            all_bottom_loc_invalid_dots_[i].clear();
        }
        all_bottom_loc_invalid_dots_.clear();
    }
    if (all_bottom_loc_coated_pmt_dots_.size() > 0){
        for (size_t i = 0; i < all_bottom_loc_coated_pmt_dots_.size(); i++){
            all_bottom_loc_coated_pmt_dots_[i].clear();
        }
        all_bottom_loc_coated_pmt_dots_.clear();
    }
    if (all_side_loc_valid_dots_.size() > 0){
        for (size_t i = 0; i < all_side_loc_valid_dots_.size(); i++){
            all_side_loc_valid_dots_[i].clear();
        }
        all_side_loc_valid_dots_.clear();
    }
    if (all_side_loc_invalid_dots_.size() > 0){
        for (size_t i = 0; i < all_side_loc_invalid_dots_.size(); i++){
            all_side_loc_invalid_dots_[i].clear();
        }
        all_side_loc_invalid_dots_.clear();
    }
    if (all_side_loc_coated_pmt_dots_.size() > 0){
        for (size_t i = 0; i < all_side_loc_coated_pmt_dots_.size(); i++){
            all_side_loc_coated_pmt_dots_[i].clear();
        }
        all_side_loc_coated_pmt_dots_.clear();
    }

    std::vector<double> this_loc_info (9);
    std::vector<double> loc_to_pmt_photon_yields;
    std::vector<double> loc_to_pmt_photon_propagation_times;
    // so we've parsed our pmt info, but now we need to get our secondary locations
    // let's start with chunking up the sides of the detector
    size_t n_chunks_c = (size_t) cylinder_circumference / desired_chunk_width_;
    size_t n_chunks_z = (size_t) cylinder_height / desired_chunk_height_;

    std::vector<double> possible_circumference_positions(n_chunks_c);
    for (size_t i = 0; i < n_chunks_c; i++){
        double this_circ = -(cylinder_circumference / 2.0) + ((cylinder_circumference) * (double)i / (double)n_chunks_c);
        possible_circumference_positions.at(i) = this_circ;
    }

    std::vector<double> possible_z_positions(n_chunks_z);
    for (size_t i = 0; i < n_chunks_z; i++){
        double this_z = (-cylinder_height / 2.0) + ((cylinder_height) * (double)i / (double)n_chunks_z);
        possible_z_positions[i] = this_z;
    }

    double actual_chunk_width = std::abs(possible_circumference_positions.at(1) - possible_circumference_positions.at(0));
    double actual_chunk_height = std::abs(possible_z_positions.at(1) - possible_z_positions.at(0));
    double area_side_chunks = actual_chunk_width * actual_chunk_height;

    double chunk_theta = (M_PI / 4.0); // half angle of corner of chunk in radians (should always be pi/4 or ~0.78)
    double chunk_diagonal_distance = actual_chunk_height / std::sin(chunk_theta);
    double chunk_half_diagonal_distance = chunk_diagonal_distance / 2.0;

    for (size_t i = 0; i < possible_circumference_positions.size(); i++){
        possible_circumference_positions[i] += (actual_chunk_width / 2.0);
    }
    for (size_t i = 0; i < possible_z_positions.size(); i++){
        possible_z_positions[i] += (actual_chunk_height / 2.0);
    }


    // setting up some variables we will need
    double small_percent_pmt_occluded = 0.0;
    I3Vector<I3Vector<double>> this_chunk_valid_points;
    I3Vector<I3Vector<double>> this_chunk_invalid_points;
    I3Vector<I3Vector<double>> this_chunk_coated_pmt_points;

    bool top_bottom_flag = false;
    size_t total_valid_points;
    size_t total_invalid_points;
    size_t total_coated_pmt_points;
    double ratio_valid;
    double ratio_coated;

    // now let's calculate x and y from these circumference positions
    for (size_t i = 0; i < possible_circumference_positions.size(); i++){
        double loc_c = possible_circumference_positions.at(i);
        double loc_theta = loc_c / cylinder_radius;
        double loc_x = cylinder_radius * std::cos(loc_theta);
        double loc_y = -cylinder_radius * std::sin(loc_theta);

        // now let's iterate over possible z positions
        for (size_t j = 0; j < possible_z_positions.size(); j++){
            area_sides += area_side_chunks;
            double loc_z = possible_z_positions.at(j);

            // let's also get the edges of our chunks for calculating overlap
            double x0 = loc_c - (actual_chunk_width / 2);
            double x1 = loc_c + (actual_chunk_width / 2);
            double y0 = loc_z + (actual_chunk_height / 2);
            double y1 = loc_z - (actual_chunk_height / 2);

            // let's call our combined_area_penalty function
            // first clear out some vectors
            if (this_chunk_valid_points.size() > 0){
                for(size_t i = 0; i < this_chunk_valid_points.size(); i++){
                    this_chunk_valid_points[i].clear();
                }
                this_chunk_valid_points.clear();
            }
            if (this_chunk_invalid_points.size() > 0){
                for(size_t i = 0; i < this_chunk_invalid_points.size(); i++){
                    this_chunk_invalid_points[i].clear();
                }
                this_chunk_invalid_points.clear();
            }
            if (this_chunk_coated_pmt_points.size() > 0){
                for(size_t i = 0; i < this_chunk_coated_pmt_points.size(); i++){
                    this_chunk_coated_pmt_points[i].clear();
                }
                this_chunk_coated_pmt_points.clear();
            }

            // now call combined_area_penalty, passing the appropriate pmts
            combined_area_penalty(loc_c, loc_z, x0, x1, y0, y1, top_bottom_flag, coated_pmt_locs_side_, uncoated_pmt_locs_side_, cylinder_radius,
                                  pmt_radius, chunk_half_diagonal_distance, this_chunk_valid_points, this_chunk_invalid_points, this_chunk_coated_pmt_points);

            // now save for debugging
            if (this_chunk_valid_points.size() > 0){
                for(size_t i = 0; i < this_chunk_valid_points.size(); i++){
                    all_side_loc_valid_dots_.push_back(this_chunk_valid_points[i]);
                }
            }
            if (this_chunk_invalid_points.size() > 0){
                for(size_t i = 0; i < this_chunk_invalid_points.size(); i++){
                    all_side_loc_invalid_dots_.push_back(this_chunk_invalid_points[i]);
                }
            }
            if (this_chunk_coated_pmt_points.size() > 0){
                for(size_t i = 0; i < this_chunk_coated_pmt_points.size(); i++){
                    all_side_loc_coated_pmt_dots_.push_back(this_chunk_coated_pmt_points[i]);
                }
            }
            // now go through our lists of valid and invalid points and tally up
            total_valid_points = this_chunk_valid_points.size();
            total_invalid_points = this_chunk_invalid_points.size();
            total_coated_pmt_points = this_chunk_coated_pmt_points.size();

            ratio_valid = (double)total_valid_points / ((double)total_valid_points + (double)total_invalid_points);
            ratio_coated = (double)total_coated_pmt_points / (double)total_valid_points;

            // now let's save if enough of this chunk is valid
            if (ratio_valid > small_percent_pmt_occluded){
                double valid_area = area_side_chunks * ratio_valid;
                double pmt_portion = portion_light_reflected_by_tpb_;
                if (ratio_coated > 0.0){
                    pmt_portion *= (0.5 * ratio_coated);
                }
                area_sides_valid += valid_area;

                // calculate facing direction
                double facing_radius = std::pow(loc_x, 2) + std::pow(loc_y, 2);
                double facing_r_x = - loc_x / facing_radius;
                double facing_r_y = - loc_y / facing_radius;
                double facing_dir_norm_factor = std::sqrt(std::pow(facing_r_x, 2) + std::pow(facing_r_y, 2));
                double facing_dir_x = facing_r_x / facing_dir_norm_factor;
                double facing_dir_y = facing_r_y / facing_dir_norm_factor;
                double facing_dir_z = 0.0;

                // now a vector to save things to
                side_chunks_counter += 1;
                this_loc_info.clear();
                this_loc_info.push_back(loc_x);
                this_loc_info.push_back(loc_y);
                this_loc_info.push_back(loc_z);
                this_loc_info.push_back(facing_dir_x);
                this_loc_info.push_back(facing_dir_y);
                this_loc_info.push_back(facing_dir_z);
                this_loc_info.push_back(pmt_portion);
                this_loc_info.push_back(valid_area);
                this_loc_info.push_back(valid_area * chunk_side_area_factor);
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
    size_t n_chunks_x = (size_t) (cylinder_max_x - cylinder_min_x) / desired_chunk_width_;
    size_t n_chunks_y = (size_t) (cylinder_max_y - cylinder_min_y) / desired_chunk_height_;

    std::vector<double> possible_x_positions(n_chunks_x);
    for (size_t i = 0; i < n_chunks_x; i++){
        double x = cylinder_min_x + ((cylinder_max_x - cylinder_min_x) * (double)i / (double)n_chunks_x);
        possible_x_positions[i] = x;
    }

    std::vector<double> possible_y_positions(n_chunks_y);
    for (size_t i = 0; i < n_chunks_y; i++){
        double y = cylinder_min_y + ((cylinder_max_y - cylinder_min_y) * (double)i / (double)n_chunks_y);
        possible_y_positions[i] = y;
    }

    actual_chunk_width = std::abs(possible_x_positions.at(1) - possible_x_positions.at(0));
    actual_chunk_height = std::abs(possible_y_positions.at(1) - possible_y_positions.at(0));
    double area_face_chunks = actual_chunk_width * actual_chunk_height;

    chunk_theta = (M_PI / 4.0); // half angle of corner of chunk in radians (should always be pi/4 or ~0.78)
    chunk_diagonal_distance = actual_chunk_height / std::sin(chunk_theta);
    chunk_half_diagonal_distance = chunk_diagonal_distance / 2.0;
    double max_radius_partially_enclosed = cylinder_radius + chunk_half_diagonal_distance;

    for(size_t i = 0; i < possible_x_positions.size(); i++){
        possible_x_positions[i] += (actual_chunk_width / 2.0);
    }
    for(size_t i = 0; i < possible_y_positions.size(); i++){
        possible_y_positions[i] += (actual_chunk_height / 2.0);
    }

    // define some things we will need
    // first clear out some vectors
    if (this_chunk_valid_points.size() > 0){
        for(size_t i = 0; i < this_chunk_valid_points.size(); i++){
            this_chunk_valid_points[i].clear();
        }
        this_chunk_valid_points.clear();
    }
    if (this_chunk_invalid_points.size() > 0){
        for(size_t i = 0; i < this_chunk_invalid_points.size(); i++){
            this_chunk_invalid_points[i].clear();
        }
        this_chunk_invalid_points.clear();
    }
    if (this_chunk_coated_pmt_points.size() > 0){
        for(size_t i = 0; i < this_chunk_coated_pmt_points.size(); i++){
            this_chunk_coated_pmt_points[i].clear();
        }
        this_chunk_coated_pmt_points.clear();
    }

    double z_top = cylinder_max_z;
    double z_bottom = cylinder_min_z;
    std::vector<double> possible_z = {z_top, z_bottom};
    top_bottom_flag = true;

    for (size_t x_it = 0; x_it < possible_x_positions.size(); x_it ++){
        for (size_t y_it = 0; y_it < possible_y_positions.size(); y_it ++){
            // now let's make the sure this x, y combination is on the circle
            double this_x = possible_x_positions.at(x_it);
            double this_y = possible_y_positions.at(y_it);
            double radius = std::sqrt(std::pow(this_x, 2) + std::pow(this_y, 2));
            if (radius > max_radius_partially_enclosed){
                // these chunks are definitely outside our circle
                continue;
            }

            // let's also get the edges of our chunks for calculating overlap
            double x0 = this_x - (actual_chunk_width / 2);
            double x1 = this_x + (actual_chunk_width / 2);
            double y0 = this_y + (actual_chunk_height / 2);
            double y1 = this_y - (actual_chunk_height / 2);

            // now let's loop over the two possibilties of being on top or bottom of the detector
            for(size_t z_it = 0; z_it < possible_z.size(); z_it ++){
                double this_z = possible_z[z_it];
                // let's call our combined_area_penalty function
                // first clear out some vectors
                if (this_chunk_valid_points.size() > 0){
                    for(size_t i = 0; i < this_chunk_valid_points.size(); i++){
                        this_chunk_valid_points[i].clear();
                    }
                    this_chunk_valid_points.clear();
                }
                if (this_chunk_invalid_points.size() > 0){
                    for(size_t i = 0; i < this_chunk_invalid_points.size(); i++){
                        this_chunk_invalid_points[i].clear();
                    }
                    this_chunk_invalid_points.clear();
                }
                if (this_chunk_coated_pmt_points.size() > 0){
                    for(size_t i = 0; i < this_chunk_coated_pmt_points.size(); i++){
                        this_chunk_coated_pmt_points[i].clear();
                    }
                    this_chunk_coated_pmt_points.clear();
                }

                // now call combined_area_penalty, passing the appropriate pmts
                if (this_z == z_top){
                    area_top += area_face_chunks;
                    combined_area_penalty(this_x, this_y, x0, x1, y0, y1, top_bottom_flag, coated_pmt_locs_top_, uncoated_pmt_locs_top_, cylinder_radius,
                                          pmt_radius, chunk_half_diagonal_distance, this_chunk_valid_points, this_chunk_invalid_points, this_chunk_coated_pmt_points);
                    // now save for debugging
                    if (this_chunk_valid_points.size() > 0){
                        for(size_t i = 0; i < this_chunk_valid_points.size(); i++){
                            all_top_loc_valid_dots_.push_back(this_chunk_valid_points[i]);
                        }
                    }
                    if (this_chunk_invalid_points.size() > 0){
                        for(size_t i = 0; i < this_chunk_invalid_points.size(); i++){
                            all_top_loc_invalid_dots_.push_back(this_chunk_invalid_points[i]);
                        }
                    }
                    if (this_chunk_coated_pmt_points.size() > 0){
                        for(size_t i = 0; i < this_chunk_coated_pmt_points.size(); i++){
                            all_top_loc_coated_pmt_dots_.push_back(this_chunk_coated_pmt_points[i]);
                        }
                    }
                }
                if (this_z == z_bottom){
                    area_bottom += area_face_chunks;
                    combined_area_penalty(this_x, this_y, x0, x1, y0, y1, top_bottom_flag, coated_pmt_locs_bottom_, uncoated_pmt_locs_bottom_, cylinder_radius,
                                          pmt_radius, chunk_half_diagonal_distance, this_chunk_valid_points, this_chunk_invalid_points, this_chunk_coated_pmt_points);
                    // now save for debugging
                    if (this_chunk_valid_points.size() > 0){
                        for(size_t i = 0; i < this_chunk_valid_points.size(); i++){
                            all_bottom_loc_valid_dots_.push_back(this_chunk_valid_points[i]);
                        }
                    }
                    if (this_chunk_invalid_points.size() > 0){
                        for(size_t i = 0; i < this_chunk_invalid_points.size(); i++){
                            all_bottom_loc_invalid_dots_.push_back(this_chunk_invalid_points[i]);
                        }
                    }
                    if (this_chunk_coated_pmt_points.size() > 0){
                        for(size_t i = 0; i < this_chunk_coated_pmt_points.size(); i++){
                            all_bottom_loc_coated_pmt_dots_.push_back(this_chunk_coated_pmt_points[i]);
                        }
                    }
                }

                // now go through our lists of valid and invalid points and tally up
                total_valid_points = this_chunk_valid_points.size();
                total_invalid_points = this_chunk_invalid_points.size();
                total_coated_pmt_points = this_chunk_coated_pmt_points.size();

                ratio_valid = (double)total_valid_points / ((double)total_valid_points + (double)total_invalid_points);
                ratio_coated = (double)total_coated_pmt_points / (double)total_valid_points;

                // now let's save if enough of this chunk is valid
                if (ratio_valid > small_percent_pmt_occluded){
                    double valid_area = area_face_chunks * ratio_valid;
                    double pmt_portion = portion_light_reflected_by_tpb_;
                    if (ratio_coated > 0.0){
                        pmt_portion *= (0.5 * ratio_coated);
                    }

                    double facing_dir_x = 0.0;
                    double facing_dir_y = 0.0;
                    double facing_dir_z;
                    if (this_z == z_top){
                        facing_dir_z = -1.0;
                        area_top_valid += valid_area;
                    }
                    if (this_z == z_bottom){
                        facing_dir_z = 1.0;
                        area_bottom_valid += valid_area;
                    }

                    // now a vector to save things to
                    face_chunks_counter += 1;
                    this_loc_info.clear();
                    this_loc_info.push_back(this_x);
                    this_loc_info.push_back(this_y);
                    this_loc_info.push_back(this_z);
                    this_loc_info.push_back(facing_dir_x);
                    this_loc_info.push_back(facing_dir_y);
                    this_loc_info.push_back(facing_dir_z);
                    this_loc_info.push_back(pmt_portion);
                    this_loc_info.push_back(valid_area);
                    this_loc_info.push_back(valid_area * chunk_side_area_factor);
                    locations_to_check_information_.push_back(this_loc_info);

                    // we are also going to pre-compute the light yield from 1 visible photon from this loc to every pmt
                    // as well as the travel time for that 1 photon to go from this loc to every pmt
                    loc_to_pmt_photon_yields.clear();
                    loc_to_pmt_photon_propagation_times.clear();
                    secondary_loc_to_pmt_propagation(full_acceptance, c_cm_per_nsec, vis_index_of_refraction, pmt_quantum_efficiency,
                                                     visible_absorption_length_, pmt_parsed_information_, this_x, this_y, this_z,
                                                     loc_to_pmt_photon_yields, loc_to_pmt_photon_propagation_times);
                    locations_to_check_to_pmt_yield_.push_back(loc_to_pmt_photon_yields);
                    locations_to_check_to_pmt_travel_time_.push_back(loc_to_pmt_photon_propagation_times);
                }
            }
        }
    }
    std::cout << "total area on top face = " << area_top << " and valid area = " << area_top_valid << std::endl;
    std::cout << "total area on bottom face = " << area_bottom << " and valid area = " << area_bottom_valid << std::endl;
    std::cout << "total area on sides face = " << area_sides << " and valid area = " << area_sides_valid << std::endl;
    std::cout << "total area = " << area_top + area_bottom + area_sides << " and valid area = " << area_top_valid + area_bottom_valid + area_sides_valid << std::endl;
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
        double time_below = all_times.at(idx_below);
        double data_value_below = all_data.at(idx_below);

        double time_above = all_times.at(idx_below + 1);
        double data_value_above = all_data.at(idx_below + 1);

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
    double delta = times.at(1) - times.at(0);
    std::vector<double> LAr_light_profile (times.size());
    for (size_t t = 0; t < times.size(); t++){
        LAr_scintillation_light_integration(times.at(t), delta, n_convolution_chunks, R_s, R_t, tau_s, tau_t, tau_rec, R_TPB, tau_TPB, LAr_light_profile.at(t));
    }
    // let's normalize our light profile
    double total_light_profile_value = 0;
    for (size_t i = 0; i < LAr_light_profile.size(); i++){
        total_light_profile_value += LAr_light_profile.at(i);
    }

    std::vector<double> normalized_LAr_light_profile (LAr_light_profile.size());
    for (size_t i = 0; i < LAr_light_profile.size(); i++){
        normalized_LAr_light_profile.at(i) = LAr_light_profile.at(i) / total_light_profile_value;
    }

    // now let's convolve our light profile with a gev
    std::vector<double> light_convolved_with_gev (times.size());
    for (size_t t = 0; t < times.size(); t++){
        convolve_light_with_gev(times.at(t), normalized_LAr_light_profile, times, delta, n_convolution_chunks, sigma, mu, xi, light_convolved_with_gev.at(t));
    }

    // let's normalize our final light profile
    double total_light_convolved_with_gev_value = 0;
    for (size_t i = 0; i < light_convolved_with_gev.size(); i++){
        total_light_convolved_with_gev_value += light_convolved_with_gev.at(i);
    }

    for (size_t i = 0; i < light_convolved_with_gev.size(); i++){
        final_light_profile.at(i) = light_convolved_with_gev.at(i) / total_light_convolved_with_gev_value;
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
        // first add a lil noise
        for(size_t bin_it = 0; bin_it < this_pmt_binned_charges.size(); bin_it++){
            this_pmt_binned_charges.at(bin_it) += noise_rate_per_time_bin;
        }
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
                       light_times, light_profile, bin_centers, bin_width, noise_rate_per_time_bin, binned_charges);
}

void FrameThread(bool const & vertex_1275_flag,
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
                 std::vector<std::vector<double>> * vector_of_vertices,
                 std::vector<double> const & light_times,
                 std::vector<double> const & light_profile,
                 std::vector<double> const & bin_centers,
                 double const & bin_width,
                 double const & noise_rate_per_time_bin,
                 std::vector<std::vector<std::vector<double>>> & vector_of_vertices_binned_charges,
                 std::vector<std::vector<std::vector<double>>> & vector_of_vertices_binned_charges_squared) {

    // call simulation code
    std::vector<std::vector<double>> binned_charges;
    double n_photons_produced_this_event;
    if (vertex_1275_flag){
        n_photons_produced_this_event = n_photons_produced * 1.275;
    } else {
        n_photons_produced_this_event = n_photons_produced * 0.511;
    }

    for (size_t vertex_it = 0; vertex_it < vector_of_vertices->size(); vertex_it ++){
        if (binned_charges.size() > 0){
            for (size_t i = 0; i < binned_charges.size(); i++){
                binned_charges.at(i).clear();
            }
            binned_charges.clear();
        }
        put_simulation_steps_together(full_acceptance, c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, quantum_efficiency,
                                      n_photons_produced_this_event, UV_absorption_length, vis_absorption_length, pmt_parsed_information_, locations_to_check_information_,
                                      locations_to_check_to_pmt_yield_, locations_to_check_to_pmt_travel_time_,
                                      vector_of_vertices->at(vertex_it).at(0), vector_of_vertices->at(vertex_it).at(1), vector_of_vertices->at(vertex_it).at(2),
                                      light_times, light_profile, bin_centers, bin_width, noise_rate_per_time_bin, binned_charges);

        // now let's add binned_charges to vector_of_vertices_summed_binned_charges!
        for (size_t pmt_it = 0; pmt_it < binned_charges.size(); pmt_it ++){
            for (size_t time_bin_it = 0; time_bin_it < binned_charges[0].size(); time_bin_it ++){
                vector_of_vertices_binned_charges[pmt_it][time_bin_it][vertex_it] = binned_charges[pmt_it][time_bin_it];
                vector_of_vertices_binned_charges_squared[pmt_it][time_bin_it][vertex_it] = std::pow(binned_charges[pmt_it][time_bin_it], 2);
            }
        }

    }

}


size_t findNearestIndex(std::vector<double> const & vec, double targetValue) {

    size_t nearestIndex = 0;
    double minDifference = std::abs(vec.at(0) - targetValue);

    for (size_t i = 1; i < vec.size(); ++i) {
        double difference = std::abs(vec.at(i) - targetValue);
        if (difference < minDifference) {
            minDifference = difference;
            nearestIndex = i;
        }
    }

    return nearestIndex;
}

I3Vector<I3Vector<I3Vector<double>>> PhotonPropagation::GetSimulation(double const & singlet_ratio_,
//I3Vector<I3Vector<double>> PhotonPropagation::GetSimulation(double const & singlet_ratio_,
//double PhotonPropagation::GetSimulation(double const & singlet_ratio_,
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
    size_t ending_event_idx = 0;

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
                            TPB_ratio_,
                            TPB_tau_,
                            smearing_sigma_,
                            smearing_mu_,
                            smearing_xi_,
                            n_convolution_chunks_,
                            times,
                            light_profile);


    // ok so we've seen our geometry file, pre-computed the lists of pmt info and locations to check for light propagation
    // we've also pre-compute the verticies of each event we want to simulate
    // and we just computed our light profile....all that's left to do is propagate our light!
    // we have the parameter total_events_that_escaped which says how many events we are simulating

    // so now we can loop over our vector containing vector of verticies to simulate
    // let's make out final vector that we will be saving to
    I3Vector<I3Vector<I3Vector<double>>> events_binned_charges(n_pmts_to_simulate, I3Vector<I3Vector<double>>(bin_centers.size(), I3Vector<double>(total_events_that_escaped)));
                                                                            // dimensions are n_pmts x n_time_bins x n_events
    std::vector<std::vector<std::vector<double>>> events_squared_binned_charges(n_pmts_to_simulate, std::vector<std::vector<double>>(bin_centers.size(), std::vector<double>(total_events_that_escaped)));
                                                                            // dimensions are n_pmts x n_time_bins x n_events

    //std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    // loop over thread_verticies_ and dispatch one thread per vector of vertices
    // note -- we added 1275kev vertices to thread_verticies_ first, so can use that to keep track of if we're dealing with a 1275 or 511 event
    for (size_t vert_it = 0; vert_it < thread_verticies_.size(); vert_it ++){
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
				    PhotonPropagationJob * job = running_jobs.at(i);
                    running_jobs.erase(running_jobs.begin() + i);
                    free_jobs.push_back(job);
                    //job->thread.join();
                    results.at(job->vector_of_vertices_index - min_vertex_idx).done = true;
                } else {
                    PhotonPropagationJob * job = running_jobs.at(i);
                }
            }
            // Check for any done results and push the corresponding frames
            size_t results_done = 0;
            for(size_t i=0; i<results.size(); ++i) {
                if(results.at(i).done) {
                    // let's save to all_events_binned_charges
                    for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
                        for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
                            for (size_t event_it = results.at(i).event_start_idx; event_it < results.at(i).event_end_idx; event_it ++){
                                events_binned_charges.at(pmt_it).at(time_bin_it).at(event_it) =
                                    results.at(i).vector_of_vertices_binned_charges->at(pmt_it).at(time_bin_it).at(event_it - results.at(i).event_start_idx);
                                events_squared_binned_charges.at(pmt_it).at(time_bin_it).at(event_it) =
                                    results.at(i).vector_of_vertices_binned_charges_squared->at(pmt_it).at(time_bin_it).at(event_it - results.at(i).event_start_idx);
                            }
                        }
                    }
                    results.at(i).vector_of_vertices_binned_charges = nullptr;
                    results.at(i).vector_of_vertices_binned_charges_squared = nullptr;
                    results.at(i).event_start_idx = 0;
                    results.at(i).event_end_idx = 0;
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
            }

            if(job != nullptr and results.size() < max_cached_vertices) {
                job->running.store(true);
                running_jobs.push_back(job);
                job->vector_of_vertices_index = vert_it;
                if (job->vector_of_vertices_index < thread_1275_verticies_.size()){
                    job->vertex_1275_flag = true;
                }
                else{
                    job->vertex_1275_flag = false;
                }
                job->vector_of_vertices = &thread_verticies_.at(job->vector_of_vertices_index);
                job->event_start_idx = ending_event_idx;
                job->event_end_idx = job->event_start_idx + job->vector_of_vertices->size();
                ending_event_idx = job->event_end_idx;
                job->vector_of_vertices_binned_charges =
                    std::make_shared<std::vector<std::vector<std::vector<double>>>>(n_pmts_to_simulate,
                            std::vector<std::vector<double>>(bin_centers.size(), std::vector<double>(job->event_end_idx - job->event_start_idx)));
                job->vector_of_vertices_binned_charges_squared =
                    std::make_shared<std::vector<std::vector<std::vector<double>>>>(n_pmts_to_simulate,
                            std::vector<std::vector<double>>(bin_centers.size(), std::vector<double>(job->event_end_idx - job->event_start_idx)));
                results.emplace_back();
                results.back().vector_of_vertices_binned_charges = job->vector_of_vertices_binned_charges;
                results.back().vector_of_vertices_binned_charges_squared = job->vector_of_vertices_binned_charges_squared;
                results.back().event_start_idx = job->event_start_idx;
                results.back().event_end_idx = job->event_end_idx;
                results.back().done = false;
                FrameThread(job->vertex_1275_flag, full_acceptance, c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction,
                            pmt_quantum_efficiency, n_photons_produced_, UV_absorption_length_, visible_absorption_length_,
                            pmt_parsed_information_, locations_to_check_information_, locations_to_check_to_pmt_yield_,
                            locations_to_check_to_pmt_travel_time_, job->vector_of_vertices, times, light_profile, bin_centers, bin_width, noise_rate_per_time_bin,
                            *job->vector_of_vertices_binned_charges, *job->vector_of_vertices_binned_charges_squared);
                job->running.store(false);
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
                PhotonPropagationJob * job = running_jobs.at(i);
                running_jobs.erase(running_jobs.begin() + i);
                free_jobs.push_back(job);
                //job->thread.join();
                results.at(job->vector_of_vertices_index - min_vertex_idx).done = true;
            }
        }

        // Check for any done results and push the corresponding frames
        size_t results_done = 0;
        for(size_t i=0; i<results.size(); ++i) {
            if(results.at(i).done) {
                // let's save to all_events_binned_charges
                for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
                    for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
                        for (size_t event_it = results.at(i).event_start_idx; event_it < results.at(i).event_end_idx; event_it ++){
                            events_binned_charges.at(pmt_it).at(time_bin_it).at(event_it) =
                                results.at(i).vector_of_vertices_binned_charges->at(pmt_it).at(time_bin_it).at(event_it - results.at(i).event_start_idx);
                            events_squared_binned_charges.at(pmt_it).at(time_bin_it).at(event_it) =
                                results.at(i).vector_of_vertices_binned_charges_squared->at(pmt_it).at(time_bin_it).at(event_it - results.at(i).event_start_idx);
                        }
                    }
                }
                results.at(i).vector_of_vertices_binned_charges = nullptr;
                results.at(i).vector_of_vertices_binned_charges_squared = nullptr;
                results.at(i).event_start_idx = 0;
                results.at(i).event_end_idx = 0;
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
        for(PhotonPropagationJob * obj : free_jobs){
            //delete obj->vector_of_vertices;
            delete obj;
        }
    }

    // let's sum yields across events
    I3Vector<I3Vector<double>> summed_over_events_binned_charges(n_pmts_to_simulate, I3Vector<double>(bin_centers.size(), 0.0)); // dimensions are n_pmts x n_time_bins
    std::vector<std::vector<double>> summed_over_squared_events_binned_charges(n_pmts_to_simulate, std::vector<double>(bin_centers.size(), 0.0)); // dimensions are n_pmts x n_time_bins
    double this_pmt_this_time_summed_across_events;
    double this_pmt_this_time_summed_across_events_squared;
    for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
        for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
            this_pmt_this_time_summed_across_events = 0;
            this_pmt_this_time_summed_across_events_squared = 0;
            for (size_t event_it = 0; event_it < total_events_that_escaped; event_it ++){
                this_pmt_this_time_summed_across_events += events_binned_charges[pmt_it][time_bin_it][event_it];
                this_pmt_this_time_summed_across_events_squared += events_squared_binned_charges[pmt_it][time_bin_it][event_it];
            }
            summed_over_events_binned_charges[pmt_it][time_bin_it] = this_pmt_this_time_summed_across_events;
            summed_over_squared_events_binned_charges[pmt_it][time_bin_it] = this_pmt_this_time_summed_across_events_squared;
        }
    }
    return events_binned_charges;

    //return summed_over_events_binned_charges;
    //return events_binned_charges;
    //std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    //std::cout << "finished simulating " << total_events_that_escaped << " events in " << std::chrono::duration_cast<std::chrono::milliseconds> (end - begin).count() << ".at(ms)" << std::endl;

    // let's print out the averaged summed wf just to make sure i didnt fuck anything up
    double total_charge_in_this_time_bin;
    std::cout << "for " << total_events_that_escaped << " total events, or " << equivalent_events_in_data << " equivalent events in data, avg summed wf = " << std::endl;
    for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
        total_charge_in_this_time_bin = 0;
        for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
            total_charge_in_this_time_bin += summed_over_events_binned_charges.at(pmt_it).at(time_bin_it);
        }
        total_charge_in_this_time_bin /= total_events_that_escaped;
        std::cout << total_charge_in_this_time_bin << ", " << std::endl;
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
            total_charge_per_time_bin += summed_over_events_binned_charges.at(pmt_it).at(time_bin_it);
        }
        if (total_charge_per_time_bin > simulation_max_value){
            simulation_max_value = total_charge_per_time_bin;
            simulation_time_of_max = bin_centers.at(time_bin_it);
        }
    }

    // let's subtract off the time of the max val
    std::vector<double> simulation_times;
    for (size_t time_bin_it = 0; time_bin_it < bin_centers.size(); time_bin_it ++){
        simulation_times.push_back(bin_centers.at(time_bin_it) - simulation_time_of_max);
    }
    // print out simulation to check
    //std::cout << "printing simulation to check : " << std::endl;
    //total_charge_per_time_bin = 0;
    //for (size_t time_bin_it = 0; time_bin_it < simulation_times.size(); time_bin_it ++){
    //    total_charge_per_time_bin = 0;
    //    for (size_t pmt_it = 0; pmt_it < summed_over_events_binned_charges.size(); pmt_it ++){
    //        total_charge_per_time_bin += summed_over_events_binned_charges[pmt_it][time_bin_it];
    //    }
    //    std::cout << "at time = " << simulation_times[time_bin_it] << " summed charge = " << total_charge_per_time_bin << std::endl;
    //}

    double min_time_to_evaluate = -6.0;  // relative to summed wf peak times
    double max_time_to_evaluate = 60;

    double total_nllh = 0.0;

    // now we can loop over our simulation, check that we are within the time bands to calculate the nllh, then calculate it!
    double mu;
    double sigma_squared;
    double k;
    size_t this_time_bin_in_data_idx;
    bool simulation_data_present;
    size_t this_time_simulation_idx;
    size_t this_time_data_idx;
    size_t simulation_bins_without_data = 0;

    // now instead of looping over simulation times and seeing if they are within our roi, we loop over roi and find corresponding simulation time (or noise if not filled)
    for (double time_roi_it = min_time_to_evaluate; time_roi_it <= max_time_to_evaluate; time_roi_it += 2.0){
        simulation_data_present = true;
        // let's check if our time is out of bounds of our simulation
        if (time_roi_it < simulation_times.front() or time_roi_it > simulation_times.back()){
            // ok we have no simulation data for this time, let's use our noise rate * n events
            simulation_data_present = false;
            simulation_bins_without_data += 1;
        } else {
            this_time_simulation_idx = findNearestIndex(simulation_times, time_roi_it);
        }
        this_time_data_idx = findNearestIndex(times_of_data_points_, time_roi_it);
        for (size_t pmt_it = 0; pmt_it < n_pmts_to_simulate; pmt_it ++){
            // so now we are looping over each pmt, we need to calculate the nllh for each time bin for each pmt
            k = data_series_[pmt_it][this_time_data_idx];
            if (simulation_data_present){
                mu = summed_over_events_binned_charges[pmt_it][this_time_simulation_idx] * (n_data_samples_ / equivalent_events_in_data);
                sigma_squared = summed_over_squared_events_binned_charges[pmt_it][this_time_simulation_idx] * std::pow((n_data_samples_ / equivalent_events_in_data), 2);
            } else {
                mu = (noise_rate_per_time_bin * equivalent_events_in_data) * (n_data_samples_ / equivalent_events_in_data);
                sigma_squared = std::pow(noise_rate_per_time_bin * equivalent_events_in_data, 2) * std::pow((n_data_samples_ / equivalent_events_in_data), 2);
            }
            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
        }

    }

    size_t max_simulation_bins_without_data = 3;

    if (simulation_bins_without_data > max_simulation_bins_without_data){
        // oops! too many bins without data! set the nllh to - inf
        total_nllh = (double)(-std::numeric_limits<double>::infinity());
    }


    //return total_nllh;

}
