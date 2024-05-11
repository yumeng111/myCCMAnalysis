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
#include <analytic-light-yields/YieldsPerPMT.h>

void combined_area_penalty(double const & chunk_center_x,
                           double const & chunk_center_y,
                           double const & x0,
                           double const & x1,
                           double const & y0,
                           double const & y1,
                           bool const & top_bottom_flag,
                           std::vector<std::vector<double>> const & relevant_coated_pmt_locs,
                           std::vector<std::vector<double>> const & relevant_uncoated_pmt_locs,
                           double const & detector_radius,
                           double const & pmt_radius,
                           double const & half_chunk_hypotenuse,
                           std::vector<std::vector<double>> & this_chunk_valid_points,
                           std::vector<std::vector<double>> & this_chunk_invalid_points,
                           std::vector<std::vector<double>> & this_chunk_coated_pmt_points) {

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

    std::vector<double> this_point;
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

void get_light_yield(double const & omega,
                     double const & distance_travelled,
                     double const & absorption_length,
                     double & light_yield) {

    light_yield = omega * std::exp(- distance_travelled / absorption_length);
}

YieldsPerPMT::YieldsPerPMT(I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height) :
    portion_light_reflected_by_tpb_(portion_light_reflected_by_tpb)
{
    SetGeoFrame(geo_frame);
    SetChunks(desired_chunk_width, desired_chunk_height);
}

void YieldsPerPMT::SetGeoFrame(I3FramePtr geo_frame) {
    this->geo_frame = geo_frame;
    SetPMTInformation(geo_frame);
}

void YieldsPerPMT::SetPMTInformation(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }

    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    std::map<CCMPMTKey, CCMOMGeo> const & pmt_geo = geo.pmt_geo;
    std::vector<double> this_pmt_loc;

    for(std::pair<CCMPMTKey const, CCMOMGeo> const & it : pmt_geo) {
        bool this_pmt_status = true;
        bool coating_status = true;
        bool side_status = false;

        CCMPMTType omtype = it.second.omtype;
        double coating_flag;
        if (omtype == coated_omtype) {
            coating_flag = 1.0; // coating flag = 1 ==> pmt is coated!
        }
        else if (omtype == uncoated_omtype) {
            coating_flag = 0.0; // uncoated pmt!
            coating_status = false;
        }
        else {
            continue; // this pmt is not an 8in! we dont need it!
        }
        // now let's compute our facing direction
        I3Position position = it.second.position;
        double pos_x = position.GetX() / I3Units::cm;
        double pos_y = position.GetY() / I3Units::cm;
        double pos_z = position.GetZ() / I3Units::cm;

        double facing_dir_x;
        double facing_dir_y;
        double facing_dir_z;

        if (pos_z > 50.0) {
            // region 0 pmts!
            facing_dir_x = 0.0;
            facing_dir_y = 0.0;
            facing_dir_z = -1.0;
        }
        if (pos_z < -50.0) {
            // region 6 pmts!
            facing_dir_x = 0.0;
            facing_dir_y = 0.0;
            facing_dir_z = 1.0;
        }
        if (pos_z < 50.0 and pos_z > -50.0) {
            double facing_radius = std::pow(pos_x, 2) + std::pow(pos_y, 2);
            double facing_r_x = pos_x / facing_radius;
            double facing_r_y = pos_y / facing_radius;
            double facing_dir_norm_factor = std::sqrt(std::pow(facing_r_x, 2) + std::pow(facing_r_y, 2));
            facing_dir_x = facing_r_x / facing_dir_norm_factor;
            facing_dir_y = facing_r_y / facing_dir_norm_factor;
            facing_dir_z = 0.0;
            side_status = true;
        }

        // now time save!!!
        yields_pmt_info this_pmt_info;
        this_pmt_info.key = it.first;
        this_pmt_info.pmt_x = pos_x;
        this_pmt_info.pmt_y = pos_y;
        this_pmt_info.pmt_z = pos_z;
        this_pmt_info.facing_dir_x = facing_dir_x;
        this_pmt_info.facing_dir_y = facing_dir_y;
        this_pmt_info.facing_dir_z = facing_dir_z;
        this_pmt_info.coating_flag = coating_flag;
        this_pmt_info.pmt_facing_area = pmt_facing_area;
        this_pmt_info.pmt_side_area = pmt_side_area;
        pmt_parsed_information_.push_back(this_pmt_info);

        // let's save a few extra things useful for chunking up detector
        this_pmt_loc.clear();
        this_pmt_loc.push_back(pos_x);
        this_pmt_loc.push_back(pos_y);
        if (pos_z == 58.0){
            if (coating_flag == 0.0){
                uncoated_pmt_locs_top_.push_back(this_pmt_loc);
            }
            if (coating_flag == 1.0){
                coated_pmt_locs_top_.push_back(this_pmt_loc);
            }
        }
        if (pos_z == -58.0){
            if (coating_flag == 0.0){
                uncoated_pmt_locs_bottom_.push_back(this_pmt_loc);
            }
            if (coating_flag == 1.0){
                coated_pmt_locs_bottom_.push_back(this_pmt_loc);
            }
        }
        if (pos_z > -58.0 and pos_z < 58.0){
            // for pmts on the sides of the detector, we actually want to save pmt_c and pos_z
            double pmt_theta = -atan2(pos_y, pos_x);
            pmt_theta = fmod(M_PI + pmt_theta, 2.0 * M_PI) - M_PI;
            double pmt_c = cylinder_radius * pmt_theta;

            this_pmt_loc.clear();
            this_pmt_loc.push_back(pmt_c);
            this_pmt_loc.push_back(pos_z);

            if (coating_flag == 0.0){
                uncoated_pmt_locs_side_.push_back(this_pmt_loc);
            }
            if (coating_flag == 1.0){
                coated_pmt_locs_side_.push_back(this_pmt_loc);
            }
        }
    }
}

void YieldsPerPMT::SetChunks(double desired_chunk_width, double desired_chunk_height) {

    this->desired_chunk_width_ = desired_chunk_width;
    this->desired_chunk_height_ = desired_chunk_height;

    // so we've parsed our pmt info, but now we need to get our secondary locations
    // let's start with chunking up the sides of the detector
    size_t n_chunks_c = (size_t) cylinder_circumference / desired_chunk_width;
    size_t n_chunks_z = (size_t) cylinder_height / desired_chunk_height;

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
    std::vector<std::vector<double>> this_chunk_valid_points;
    std::vector<std::vector<double>> this_chunk_invalid_points;
    std::vector<std::vector<double>> this_chunk_coated_pmt_points;

    bool top_bottom_flag = false;
    size_t total_valid_points;
    size_t total_invalid_points;
    size_t total_coated_pmt_points;
    double ratio_valid;
    double ratio_coated;

    // now let's calculate x and y from these circumference positions
    for (size_t i = 0; i < possible_circumference_positions.size(); i++) {
        double loc_c = possible_circumference_positions.at(i);
        double loc_theta = loc_c / cylinder_radius;
        double loc_x = cylinder_radius * std::cos(loc_theta);
        double loc_y = -cylinder_radius * std::sin(loc_theta);

        // now let's iterate over possible z positions
        for (size_t j = 0; j < possible_z_positions.size(); j++){
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

            // now go through our lists of valid and invalid points and tally up
            total_valid_points = this_chunk_valid_points.size();
            total_invalid_points = this_chunk_invalid_points.size();
            total_coated_pmt_points = this_chunk_coated_pmt_points.size();

            ratio_valid = (double)total_valid_points / ((double)total_valid_points + (double)total_invalid_points);
            ratio_coated = (double)total_coated_pmt_points / (double)total_valid_points;

            // now let's save if enough of this chunk is valid
            if (ratio_valid > small_percent_pmt_occluded){
                std::vector<double> side_chunk_xy;
                side_chunk_xy.push_back(loc_c);
                side_chunk_xy.push_back(loc_z);

                double valid_area = area_side_chunks * ratio_valid;
                double pmt_portion = portion_light_reflected_by_tpb_;
                if (ratio_coated > 0.0){
                    pmt_portion *= (0.5 * ratio_coated);
                }
                // calculate facing direction
                double facing_radius = std::pow(loc_x, 2) + std::pow(loc_y, 2);
                double facing_r_x = loc_x / facing_radius;
                double facing_r_y = loc_y / facing_radius;
                double facing_dir_norm_factor = std::sqrt(std::pow(facing_r_x, 2) + std::pow(facing_r_y, 2));
                double facing_dir_x = facing_r_x / facing_dir_norm_factor;
                double facing_dir_y = facing_r_y / facing_dir_norm_factor;
                double facing_dir_z = 0.0;

                // now save!!!
                secondary_loc_info this_loc_info;
                this_loc_info.loc_x = loc_x;
                this_loc_info.loc_y = loc_y;
                this_loc_info.loc_z = loc_z;
                this_loc_info.facing_dir_x = facing_dir_x;
                this_loc_info.facing_dir_y = facing_dir_y;
                this_loc_info.facing_dir_z = facing_dir_z;
                this_loc_info.pmt_portion = pmt_portion;
                this_loc_info.loc_facing_area = valid_area;
                this_loc_info.loc_side_area = valid_area * chunk_side_area_factor;
                locations_to_check_information_.push_back(this_loc_info);
            }
        }
    }

    // time to chunk up the top and bottom of the detector
    size_t n_chunks_x = (size_t) (cylinder_max_x - cylinder_min_x) / desired_chunk_width;
    size_t n_chunks_y = (size_t) (cylinder_max_y - cylinder_min_y) / desired_chunk_height;

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
                    combined_area_penalty(this_x, this_y, x0, x1, y0, y1, top_bottom_flag, coated_pmt_locs_top_, uncoated_pmt_locs_top_, cylinder_radius,
                                          pmt_radius, chunk_half_diagonal_distance, this_chunk_valid_points, this_chunk_invalid_points, this_chunk_coated_pmt_points);
                }
                if (this_z == z_bottom){
                    combined_area_penalty(this_x, this_y, x0, x1, y0, y1, top_bottom_flag, coated_pmt_locs_bottom_, uncoated_pmt_locs_bottom_, cylinder_radius,
                                          pmt_radius, chunk_half_diagonal_distance, this_chunk_valid_points, this_chunk_invalid_points, this_chunk_coated_pmt_points);
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
                        std::vector<double> top_chunk_xy;
                        top_chunk_xy.push_back(this_x);
                        top_chunk_xy.push_back(this_y);
                        facing_dir_z = -1.0;
                    }
                    if (this_z == z_bottom){
                        std::vector<double> bottom_chunk_xy;
                        bottom_chunk_xy.push_back(this_x);
                        bottom_chunk_xy.push_back(this_y);
                        facing_dir_z = 1.0;
                    }

                    secondary_loc_info this_loc_info;
                    this_loc_info.loc_x = this_x;
                    this_loc_info.loc_y = this_y;
                    this_loc_info.loc_z = this_z;
                    this_loc_info.facing_dir_x = facing_dir_x;
                    this_loc_info.facing_dir_y = facing_dir_y;
                    this_loc_info.facing_dir_z = facing_dir_z;
                    this_loc_info.pmt_portion = pmt_portion;
                    this_loc_info.loc_facing_area = valid_area;
                    this_loc_info.loc_side_area = valid_area * chunk_side_area_factor;
                    locations_to_check_information_.push_back(this_loc_info);
                }
            }
        }
    }
}

// let's make a function to calculate how much light from a vertex goes into each PMT

void vertex_to_pmt_propagation(double const & c_cm_per_nsec,
                               double const & uv_index_of_refraction,
                               double const & vis_index_of_refraction,
                               double const & n_photons_produced,
                               double const & time_offset,
                               double const & absorption_length,
                               std::vector<yields_pmt_info> const & pmt_parsed_information_,
                               bool const & is_visible,
                               I3Position const & vertex,
                               boost::shared_ptr<PhotonYieldSummarySeriesMap> & all_pmt_yields_map,
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

        // check if this is a pmt key we should fite
        bool fit = false;
        for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++) {
            if (keys_to_fit.at(pmt_it) == pmt_key){
                fit = true;
            }
        }
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
            (*all_pmt_yields_map)[pmt_key] = PhotonYieldSummarySeries ();
        }

        else {
            double omega;
            double distance_travelled;
            // let's get solid angle and distance from vertex to this pmt
            get_solid_angle_and_distance_vertex_to_location(vertex.GetX() / I3Units::cm, vertex.GetY() / I3Units::cm, vertex.GetZ() / I3Units::cm,
                                                            pmt_x_loc, pmt_y_loc, pmt_z_loc,
                                                            facing_dir_x, facing_dir_y, facing_dir_z,
                                                            pmt_facing_area, pmt_side_area, omega, distance_travelled);

            double photons_in_this_pmt;
            // now call function to get light yields
            get_light_yield(omega, distance_travelled, absorption_length, photons_in_this_pmt);
            photons_in_this_pmt *= n_photons_produced;

            double travel_time;
            // ok so we have the photons seen by this pmt, let's get the propagation times
            if (is_visible) {
                travel_time = distance_travelled / (c_cm_per_nsec / vis_index_of_refraction); // units of nsec
            }
            else{
                travel_time = distance_travelled / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec
            }

            // now all that's left to do is save!
            PhotonYieldSummary this_pmt_yield_summary;
            this_pmt_yield_summary.time = travel_time + time_offset;
            this_pmt_yield_summary.yield = photons_in_this_pmt;
            if (is_visible){
                this_pmt_yield_summary.photon_source = PhotonYieldSummary::PhotonSource::TPBFoil;
            }
            else {
                this_pmt_yield_summary.photon_source = PhotonYieldSummary::PhotonSource::Vertex;
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

void vertex_to_TPB_to_PMT_propagation(double const & c_cm_per_nsec,
                                      double const & uv_index_of_refraction,
                                      double const & vis_index_of_refraction,
                                      double const & n_photons_produced,
                                      double const & UV_absorption_length,
                                      double const & vis_absorption_length,
                                      I3Position const & vertex,
                                      std::vector<yields_pmt_info> const & pmt_parsed_information_,
                                      std::vector<secondary_loc_info> const & locations_to_check_information_,
                                      boost::shared_ptr<PhotonYieldSummarySeriesMap> & all_pmt_yields_map,
                                      std::vector<CCMPMTKey> const & keys_to_fit) {
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
        get_light_yield(omega, distance_travelled, UV_absorption_length, photons_at_secondary_location);
        photons_at_secondary_location *= n_photons_produced;
        photons_at_secondary_location *= pmt_portion;

        // now let's get the travel time as well
        travel_time = distance_travelled / (c_cm_per_nsec / uv_index_of_refraction); // units of nsec

        // ok so now we have time offset and yields at this loc, let's calculate what gets into our PMTs
        vertex_to_pmt_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction,
                                  photons_at_secondary_location, travel_time, vis_absorption_length,
                                  pmt_parsed_information_, true, I3Position(loc_x * I3Units::cm, loc_y * I3Units::cm, loc_z * I3Units::cm),
                                  all_pmt_yields_map, keys_to_fit);

        // since vertex_to_pmt_propagation saves info -- we're done!
    }
}

void PutSimulationStepsTogether(HESodiumEvent const & soidum_event,
                                double const & c_cm_per_nsec,
                                double const & uv_index_of_refraction,
                                double const & vis_index_of_refraction,
                                double const & UV_absorption_length,
                                double const & vis_absorption_length,
                                std::vector<yields_pmt_info> const & pmt_parsed_information_,
                                std::vector<secondary_loc_info> const & locations_to_check_information_,
                                boost::shared_ptr<PhotonYieldSummarySeriesMap> & all_pmt_yields_map,
                                std::vector<CCMPMTKey> const & keys_to_fit) {

    // ok so let's grab our vertex locations from our soidum_event
    I3Position photon_vertex = soidum_event.photon_vertex;
    I3Position electron_vertex = soidum_event.electron_vertex;
    I3Position positron_vertex = soidum_event.positron_vertex;

    double n_photons_1275 = 1.275;
    double n_photons_511 = 0.511;
    double default_time_offset = 0.0;

    // now let's call our functions to get direct vertex --> PMT yields (UV scint light)
    vertex_to_pmt_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_1275,
                              default_time_offset, UV_absorption_length, pmt_parsed_information_,
                              false, photon_vertex, all_pmt_yields_map, keys_to_fit);
    vertex_to_pmt_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_511,
                              default_time_offset, UV_absorption_length, pmt_parsed_information_,
                              false, electron_vertex, all_pmt_yields_map, keys_to_fit);
    vertex_to_pmt_propagation(c_cm_per_nsec, uv_index_of_refraction, vis_index_of_refraction, n_photons_511,
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

boost::shared_ptr<PhotonYieldSummarySeriesMap> YieldsPerPMT::GetAllYields(HESodiumEvent const & event_vertex, double UV_absorption_length, std::vector<CCMPMTKey> const & keys_to_fit) {

    // so we have an event we want to simulate
    // probably want to to multi-thread at some point, but for now we will just loop over all high energy sodium events

    // let's also make our PhotonYieldSummary map to save to
    boost::shared_ptr<PhotonYieldSummarySeriesMap> all_pmt_yields_map_ = boost::make_shared<PhotonYieldSummarySeriesMap> ();

    // ok we've set up the geometry stuff, now we can take the vertex and calculate light yields in each pmt
    PutSimulationStepsTogether(event_vertex, c_cm_per_nsec_, uv_index_of_refraction_,
                               vis_index_of_refraction_, UV_absorption_length, vis_absorption_length_,
                               pmt_parsed_information_, locations_to_check_information_, all_pmt_yields_map_, keys_to_fit);

    return all_pmt_yields_map_;
}


