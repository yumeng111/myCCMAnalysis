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
#include <analytic-light-yields/GenerateExpectation.h>

GenerateExpectation::GenerateExpectation() :
    keys_to_fit(I3VectorCCMPMTKey()), geo_frame(geo_frame), n_sodium_events(n_sodium_events), desired_chunk_width(desired_chunk_width), desired_chunk_height(desired_chunk_height) {}

GenerateExpectation::GenerateExpectation(I3VectorCCMPMTKey keys_to_fit, size_t n_sodium_events, I3FramePtr geo_frame, double portion_light_reflected_by_tpb, double desired_chunk_width, double desired_chunk_height) :
    keys_to_fit(keys_to_fit), geo_frame(geo_frame), n_sodium_events(n_sodium_events), desired_chunk_width(desired_chunk_width), desired_chunk_height(desired_chunk_height) {}

void GenerateExpectation::GetSodiumVertices(size_t n_events_to_simulate, double z_position) {
    sodium_events_constructor = std::make_shared<SodiumVertexDistribution> ();
    event_vertices = sodium_events_constructor->GetEventVertices(n_events_to_simulate, z_position);
}

void GenerateExpectation::GetYieldsAndOffsets(double uv_absorption) {
    yields_per_pmt_per_event.clear();
    binned_yields.clear();
    binned_square_yields.clear();
    yields_and_offset_constructor = std::make_shared<YieldsPerPMT>(geo_frame, portion_light_reflected_by_tpb, desired_chunk_width, desired_chunk_height);
    // now loop over events and get map between CCMPMTKey and PhotonYieldSummarySeries
    for (size_t sodium_it = 0; sodium_it < event_vertices->size(); ++sodium_it) {
        boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_per_event = yields_and_offset_constructor->GetAllYields(event_vertices->at(sodium_it), uv_absorption, keys_to_fit);
        yields_per_pmt_per_event.push_back(yields_per_event);
        for (PhotonYieldSummarySeriesMap::const_iterator i = yields_per_event->begin(); i != yields_per_event->end(); i++) {
            I3Vector<PhotonYieldSummary> const & yields = i->second;
            if(yields.size() == 0) {
                continue;
            }
        }
    }
}

void GenerateExpectation::ComputeBinnedYield(CCMPMTKey key, double max_time) {
    size_t n_bins = max_time / 2.0;
    binned_yields[key] = std::vector<double>(n_bins, 0.0);
    binned_square_yields[key] = std::vector<double>(n_bins, 0.0);
    std::vector<double> & binned_yields_per_pmt = binned_yields[key];
    std::vector<double> & binned_square_yields_per_pmt = binned_square_yields[key];

    for (size_t sodium_it = 0; sodium_it < event_vertices->size(); ++sodium_it) {
        boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_per_event = yields_per_pmt_per_event.at(sodium_it);
        PhotonYieldSummarySeriesMap::const_iterator i = yields_per_event->find(key);
        if (i == yields_per_event->end()) {
            continue;
        }
        I3Vector<PhotonYieldSummary> const & yields = i->second;
        if(yields.size() == 0) {
            continue;
        }
        for(size_t yield_it = 0; yield_it < yields.size(); ++yield_it) {
            PhotonYieldSummary const & yield = yields.at(yield_it);
            size_t bin_idx = yield.time / 2.0;
            if(bin_idx >= n_bins) {
                continue;
            }
            binned_yields_per_pmt.at(bin_idx) += yield.yield;
            binned_square_yields_per_pmt.at(bin_idx) += yield.yield * yield.yield;
        }
    }
    binned_yields[key] = binned_yields_per_pmt;
    binned_square_yields[key] = binned_square_yields_per_pmt;
}

