#ifndef ZOffsetcppMinimizer_H
#define ZOffsetcppMinimizer_H

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
#include "analytic-light-yields/CalculateNLLH.h"
#include "simclasses/CCMMCPE.h"
#include "simclasses/PhotonSummary.h"

class ZOffsetcppMinimizer {

public:
    ZOffsetcppMinimizer();
    
    // some vectors of maps between CCMPMTKey and vector of doubles to save best fits to
    std::vector<I3MapPMTKeyVectorDouble> data;
    std::vector<I3MapPMTKeyVectorDouble> pred;
    std::vector<I3MapPMTKeyVectorDouble> times;

    std::vector<double> GrabNormSeed(CCMPMTKey key, double baseline_efficiency, double LPmu, double LPsigma, double LPscale, 
                                     std::vector<double> uv_abs_lengths, std::vector<double> z_offsets, std::vector<size_t> n_sodium_events,
                                     std::vector<I3MapPMTKeyDouble> time_offsets, std::vector<std::string> data_file_names); 

    std::vector<double> FitParameters(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                      I3MapPMTKeyDouble time_offset1, I3MapPMTKeyDouble time_offset2, I3MapPMTKeyDouble time_offset3, I3MapPMTKeyDouble time_offset4,
                                      std::vector<std::string> data_file_names, std::vector<double> z_offsets, std::vector<size_t> n_sodium_events,
                                      I3MapPMTKeyDouble pmt_efficiency, std::vector<bool> fit_flags);

    I3MapPMTKeyVectorDouble GetBestFitData(size_t dataset_idx) {return data.at(dataset_idx);};
    I3MapPMTKeyVectorDouble GetBestFitPred(size_t dataset_idx) {return pred.at(dataset_idx);};
    I3MapPMTKeyVectorDouble GetBestFitTimes(size_t dataset_idx) {return times.at(dataset_idx);};

};
#endif
