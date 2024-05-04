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
#include "dataclasses/physics/PhotonYieldSummary.h"

struct pmt_info {

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
    double portion_light_reflected_by_tpb_ = 1.0;
    double vis_absorption_length_ = 2000.0;
    double desired_chunk_width_ = 20.0; // use 5 for finer binning
    double desired_chunk_height_ = 20.0; // use 5 for finer binning
    std::vector<pmt_info> pmt_parsed_information_;
    std::vector<secondary_loc_info> locations_to_check_information_;

    unsigned int coated_omtype = (unsigned int)10;
    unsigned int uncoated_omtype = (unsigned int)20;
    I3Vector<I3Vector<double>> uncoated_pmt_locs_top_;
    I3Vector<I3Vector<double>> coated_pmt_locs_top_;
    I3Vector<I3Vector<double>> uncoated_pmt_locs_bottom_;
    I3Vector<I3Vector<double>> coated_pmt_locs_bottom_;
    I3Vector<I3Vector<double>> uncoated_pmt_locs_side_;
    I3Vector<I3Vector<double>> coated_pmt_locs_side_;

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

    bool set_pmt_info_ = false;
    bool set_secondary_locs_ = false;

public:
    YieldsPerPMT();
    void GetPMTInformation(I3FramePtr frame);
    void GetSecondaryLocs(double const & desired_chunk_width, double const & desired_chunk_height) ;
    boost::shared_ptr<PhotonYieldSummarySeriesMap> GetAllYields(boost::shared_ptr<HESodiumEventSeries> const & event_vertices,  I3FramePtr geo_frame, double const & UV_absorption_length);
};
#endif
