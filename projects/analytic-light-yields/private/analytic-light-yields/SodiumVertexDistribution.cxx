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

#include <icetray/I3Units.h>
#include <analytic-light-yields/SodiumVertexDistribution.h>

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

void get_1275kev_photon(double const & z_offset,
                        double const & source_diameter,
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
    double z = z_offset;

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

void get_511kev_photon(double const & z_offset,
                       double const & source_diameter,
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
    double z = z_offset;

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

SodiumVertexDistribution::SodiumVertexDistribution() {}

boost::shared_ptr<HESodiumEventSeries> SodiumVertexDistribution::GetEventVertices(size_t const & n_events_to_simulate, double const & z_position){
    // set our events parameter
    equivalent_events_in_data = 0;

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

    while (equivalent_events_in_data < n_events_to_simulate){
        // so for the high energy bump, we want 3 photons -- isotropic 1275 kev, and 2 back to back 511 kev photons
        // let's call our functions to check if 511 and 1275 kev photons escape
        get_1275kev_photon(z_position,
                           source_diameter,
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
        get_511kev_photon(z_position,
                          source_diameter,
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
            equivalent_events_in_data += 1;
            // now we can save!!!
            HESodiumEvent single_sodium_event;
            single_sodium_event.photon_vertex = I3Position(final_x * I3Units::cm, final_y * I3Units::cm, final_z * I3Units::cm);
            single_sodium_event.electron_vertex = I3Position(final_x_photon_1 * I3Units::cm, final_y_photon_1 * I3Units::cm, final_z_photon_1 * I3Units::cm);
            single_sodium_event.positron_vertex = I3Position(final_x_photon_2 * I3Units::cm, final_y_photon_2 * I3Units::cm, final_z_photon_2 * I3Units::cm);
            event_vertices->push_back(single_sodium_event);
        }
    }

    return event_vertices;
}


