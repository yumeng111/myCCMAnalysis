#ifndef StandAloneUVAbsMinimizer_H
#define StandAloneUVAbsMinimizer_H

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
#include <memory>

#include "icetray/CCMPMTKey.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "analytic-light-yields/G4YieldsPerPMT.h"

class StandAloneUVAbsMinimizer {

public:
    StandAloneUVAbsMinimizer();
    
    std::vector<double> MinimizeUVAbsorption(std::vector<CCMPMTKey> keys_to_fit, I3MapPMTKeyDouble data);
    double GrabScalingSeed(CCMPMTKey key, double data,  double uv_abs, std::shared_ptr<G4YieldsPerPMT> yields_constructor);
    
    I3MapPMTKeyDouble best_fit_summed_yields;
    I3MapPMTKeyDouble GetBestFitPred() {return best_fit_summed_yields;};
    
    void MakeSpatialDistros(std::vector<CCMPMTKey> keys_to_fit, std::vector<double> uv_abs, double scaling);
    std::vector<I3MapPMTKeyDouble> best_fit_summed_yields_multiple_uv_abs;
    I3MapPMTKeyDouble GetBestFitPredMultipleUVAbs(size_t idx) {return best_fit_summed_yields_multiple_uv_abs.at(idx); };

};
#endif
