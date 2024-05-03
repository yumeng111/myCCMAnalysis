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

GenerateExpectation::GenerateExpectation() {}

boost::shared_ptr<I3MapPMTKeyVectorDouble> GenerateExpectation::GetExpectation(boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_map, I3Vector<double> LAr_light_profile){

    // so let's recap where we are
    // 1) we generated sodium event vertices
    // 2) we took those vertices and calculated yields and time offsets in PMTs
    // 3) we calculated our light profile
    // so all that's left to do is take our yields and apply to LAr light profile then bin charge
    // after binning, we are going to be ready to compare to data!

    double light_profile_start_time = 0.0;
    size_t n_light_profile_time_bins = LAr_light_profile.size();
    I3Vector<double> light_profile_times;
    for (size_t i = 0; i < n_light_profile_time_bins; i++){
        light_profile_times.push_back(light_profile_start_time + 2.0 * (double)i);
    }

    // let's also make the time binning for our expectation
    double expectation_start_time = 0.0;
    size_t n_expectation_time_bins = n_light_profile_time_bins + 10; // making it a bit longer than light profile to account for time shifts
    I3Vector<double> expectation_times;
    for (size_t i = 0; i < n_expectation_time_bins; i++){
        expectation_times.push_back(expectation_start_time + 2.0 * (double)i);
    }

    CCMPMTKey pmt_key;
    double time_offset;
    double relative_yields;
    double light_time;
    double light_val;
    size_t bin_idx;
    for (PhotonYieldSummarySeriesMap::const_iterator i = yields_map->begin(); i != yields_map->end(); i++) {
        pmt_key = i->first;

        // check if this pmt is in final_binned_yields_ -- if not, add it
        if (final_binned_yields_->find(pmt_key) == final_binned_yields_->end()) {
            (*final_binned_yields_)[pmt_key] = std::vector<double> (n_expectation_time_bins, 0.0);
        }

        // now let's loop over all the PhotonYieldSummary for this PMT
        for(PhotonYieldSummary const & this_yield: i->second) {
            time_offset = this_yield.time;
            relative_yields = this_yield.yield;

            // now apply time offsets and yields to our light profile
            for (size_t light_time_it = 0; light_time_it < n_light_profile_time_bins; light_time_it++){
                light_time = light_profile_times.at(light_time_it) + time_offset;
                light_val = LAr_light_profile.at(light_time_it) * relative_yields;

                // now binning!
                bin_idx = (size_t) (light_time / 2.0);
                (*final_binned_yields_)[pmt_key].at(bin_idx) += light_val;
            }
        }
    }

    return final_binned_yields_;

}







