#ifndef LArScintillationLightProfile_H
#define LArScintillationLightProfile_H

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

class LArScintillationLightProfile {
public:
    LArScintillationLightProfile();
    I3Vector<double> GetFullLightProfile(double const & singlet_ratio, double const & triplet_ratio, double const & singlet_tau,
                                         double const & triplet_tau, double const & recombination_tau, double const & TPB_tau);
    I3Vector<double> GetSimplifiedLightProfile(double const & singlet_ratio, double const & triplet_ratio, double const & singlet_tau,
                                         double const & triplet_tau, double const & TPB_tau);
    I3Vector<double> GetSimplifiedLightProfileDeriv(double const & singlet_ratio, double const & triplet_ratio, double const & singlet_tau,
                                                    double const & triplet_tau, double const & TPB_tau, std::string deriv_variable);
};
#endif
