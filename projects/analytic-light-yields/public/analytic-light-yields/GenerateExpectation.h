#ifndef GenerateExpectation_H
#define GenerateExpectation_H

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
#include "dataclasses/I3Map.h"
#include "dataclasses/physics/PhotonYieldSummary.h"

class GenerateExpectation {
    // make map of final binned yields to return
    boost::shared_ptr<I3MapPMTKeyVectorDouble> final_binned_yields_ = boost::make_shared<I3MapPMTKeyVectorDouble> ();

public:
    GenerateExpectation();
    boost::shared_ptr<I3MapPMTKeyVectorDouble> GetExpectation(boost::shared_ptr<PhotonYieldSummarySeriesMap> yields_map, I3Vector<double> LAr_light_profile);
};
#endif
