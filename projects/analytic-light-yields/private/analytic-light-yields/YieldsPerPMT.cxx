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
#include "dataclasses/physics/HESodiumEvent.h"
#include <analytic-light-yields/SodiumVertexDistribution.h>

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


void YieldsPerPMT::GetPlottingInformation(CCMPMTKey key, size_t n_events_to_simulate, double z_position){
    // let's grab some event vertices
    std::shared_ptr<SodiumVertexDistribution> sodium_events_constructor = std::make_shared<SodiumVertexDistribution> ();
    boost::shared_ptr<HESodiumEventSeries> event_vertices = sodium_events_constructor->GetEventVertices(n_events_to_simulate, z_position);

    // now get all yields
    size_t n_threads = 0;
    double uv_absorption = 2800.0;
    double max_time = 100.0;
    double photons_per_mev = 10000.0;
    std::map<CCMPMTKey, std::vector<double>> binned_yields;
    std::map<CCMPMTKey, std::vector<double>> binned_square_yields;

    GetAllYields<double>(n_threads, event_vertices, uv_absorption, {key}, max_time, photons_per_mev, binned_yields, binned_square_yields);

    // now grab yields for this pmt
    yields_plotting = binned_yields.at(key);
}


