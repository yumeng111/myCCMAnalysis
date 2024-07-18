#ifndef cppMinimizer_H
#define cppMinimizer_H

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

#include "icetray/CCMPMTKey.h"
#include "dataclasses/I3Vector.h"
#include "simclasses/CCMMCPE.h"
#include "simclasses/PhotonSummary.h"

class cppMinimizer {

public:
    cppMinimizer();
    void OnePMTOneDataSetMinimization(std::vector<CCMPMTKey> all_keys_to_fit, std::vector<std::string> data_file_names, std::vector<double> z_offsets,
                                                     double norm_seed, size_t n_sodium_events, bool use_g4_yields);
    std::vector<double> OnePMTMultipleDataSetMinimization(CCMPMTKey this_key, std::vector<std::string> data_file_names, std::vector<double> z_offsets,
                        std::vector<double> time_offsets, std::vector<double> norm_seeds, size_t n_sodium_events, bool use_g4_yields);
    
    // some vectors of maps between CCMPMTKey and vector of doubles to save best fits to
    std::vector<I3MapPMTKeyVectorDouble> data;
    std::vector<I3MapPMTKeyVectorDouble> pred;
    std::vector<I3MapPMTKeyVectorDouble> times;
    
    I3MapPMTKeyVectorDouble GetBestFitData(size_t dataset_idx) {return data.at(dataset_idx);};
    I3MapPMTKeyVectorDouble GetBestFitPred(size_t dataset_idx) {return pred.at(dataset_idx);};
    I3MapPMTKeyVectorDouble GetBestFitTimes(size_t dataset_idx) {return times.at(dataset_idx);};

    std::vector<std::vector<double>> all_data_to_return;
    std::vector<double> GetBestFitAtPMTIdx(size_t pmt_idx) { return all_data_to_return.at(pmt_idx); };

};
#endif
