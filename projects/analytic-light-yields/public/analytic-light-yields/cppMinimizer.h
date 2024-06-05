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

class cppMinimizer {

public:
    cppMinimizer();
    std::vector<double> OnePMTOneDataSetMinimization(CCMPMTKey this_key, std::vector<std::string> data_file_names, std::vector<double> z_offsets);
    std::vector<double> OnePMTMultipleDataSetMinimization(CCMPMTKey this_key, std::vector<std::string> data_file_names, std::vector<double> z_offsets);
};
#endif
