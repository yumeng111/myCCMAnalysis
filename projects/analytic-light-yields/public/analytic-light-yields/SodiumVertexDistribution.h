#ifndef SodiumVertexDistribution_H
#define SodiumVertexDistribution_H

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

#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/physics/HESodiumEvent.h"

class SodiumVertexDistribution {

    // now some geometry things for throwing source events
    double cylinder_height = 48.4 * 2.54; // units of cm!
    double cylinder_radius = 40.0 * 2.54; // units of cm!
    double cylinder_min_z = - cylinder_height / 2.0;
    double rod_diameter = 1.0;
    double source_diameter = 0.8;
    double rod_width = (rod_diameter - source_diameter ) / 2;
    double source_inset = -0.25;
    double decay_constant = 10.0;
    double source_rod_lower_end_cap = - source_inset;
    double detector_lower_end_cap = cylinder_min_z;
    double detector_radius = cylinder_radius;
    double pos_rad = source_diameter / 2.0;

    // let's keep track of events we're simulation
    size_t equivalent_events_in_data;

public:
    SodiumVertexDistribution();
    boost::shared_ptr<HESodiumEventSeries> GetEventVertices(size_t const & n_events_to_simulate, double const & z_position);
};
#endif
