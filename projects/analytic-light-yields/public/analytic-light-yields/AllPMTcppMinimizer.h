#ifndef AllPMTcppMinimizer_H
#define AllPMTcppMinimizer_H

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

class AllPMTcppMinimizer {

public:
    AllPMTcppMinimizer();
    
    // some vectors of maps between CCMPMTKey and vector of doubles to save best fits to
    std::vector<I3MapPMTKeyVectorDouble> data;
    std::vector<I3MapPMTKeyVectorDouble> pred;
    std::vector<I3MapPMTKeyVectorDouble> times;

    std::vector<double> MultiplePMTMinimization(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble PMT_efficiencies,
                                                I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                                I3MapPMTKeyDouble time_offset1, I3MapPMTKeyDouble time_offset2, I3MapPMTKeyDouble time_offset3, I3MapPMTKeyDouble time_offset4,
                                                std::vector<std::string> data_file_names,
                                                std::vector<double> z_offsets, size_t n_sodium_events);        

    I3MapPMTKeyVectorDouble GetBestFitData(size_t dataset_idx) {return data.at(dataset_idx);};
    I3MapPMTKeyVectorDouble GetBestFitPred(size_t dataset_idx) {return pred.at(dataset_idx);};
    I3MapPMTKeyVectorDouble GetBestFitTimes(size_t dataset_idx) {return times.at(dataset_idx);};

};
#endif
