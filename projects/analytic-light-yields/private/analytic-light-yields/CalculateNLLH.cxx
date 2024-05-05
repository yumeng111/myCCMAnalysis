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
#include <analytic-light-yields/CalculateNLLH.h>

CalculateNLLH::CalculateNLLH() {}

double CalculateNLLH::GetNLLH(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
                              boost::shared_ptr<I3MapPMTKeyVectorDouble> nuisance_paramters, boost::shared_ptr<I3MapPMTKeyVectorDouble> data,
                              double const & n_data_events, size_t time_bin_offset){

    // let's check to see if we've initialized our GenerateExpectation constructor
    if (gen_expectation == nullptr){
        gen_expectation = new GenerateExpectation();
    }
    // let's grab our expectation
    boost::shared_ptr<I3MapPMTKeyVectorDouble> pred = gen_expectation->GetExpectation(analytic_light_yield_setup, geo_frame);

    // ok so we have data and we have prediction
    // before we can calculate our LLH, we need to figure our the timing
    // let's try aligning on event start time but allowing a time offset to float
    // then we will calculate our llh from 6nsec after event start until 70nsec after event start

    size_t llh_bins_from_start = 5;
    size_t llh_time_bins = 35;

    std::vector<double> this_pmt_pred;
    std::vector<double> this_pmt_data;
    size_t corresponding_data_bin;
    for (I3MapPMTKeyVectorDouble::const_iterator i = pred->begin(); i != pred->end(); i++) {
        // looping over each pmt key in our pred map between CCMPMTKey and std::vector<double>
        this_pmt_pred = i->second; 
        this_pmt_data = data->at(i->first);

        // now loop over the time bins in our pred
        for (size_t time_bin_it = llh_bins_from_start; time_bin_it < this_pmt_pred.size(); time_bin_it++){
            corresponding_data_bin = time_bin_it + time_bin_offset;    
        }
    }


    return 0.0;

}







