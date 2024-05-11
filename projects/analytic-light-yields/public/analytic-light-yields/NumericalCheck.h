#ifndef NumericalCheck_H
#define NumericalCheck_H

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
#include <memory>

#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/physics/PhotonYieldSummary.h"
#include "dataclasses/physics/HESodiumEvent.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"
#include <analytic-light-yields/SodiumVertexDistribution.h>
#include <analytic-light-yields/YieldsPerPMT.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>
#include <analytic-light-yields/CalculateNLLH.h>

class NumericalCheck {
    std::shared_ptr<CalculateNLLH> CalculateNLLH_constructor = nullptr;
    std::vector<I3MapPMTKeyVectorDouble> dubug_info;
public:
    NumericalCheck();
    void CheckSimplifiedLightProfile(AnalyticLightYieldGenerator analytic_light_yield_setup, double light_time_offset);
    void CheckNLLHDerivs(double const & k, double const & mu, double const & sigma_squared, double const & DmuDtheta, double const & Dsigma_squaredDtheta);
    // I3MapPMTKeyVectorDouble CheckData(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, std::vector<CCMPMTKey> keys_to_fit,
    //                                                                  I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset);
    // I3MapPMTKeyVectorDouble CheckPred(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
    //                                                                   std::vector<CCMPMTKey> keys_to_fit,
    //                                                                  I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset);
    // I3MapPMTKeyVectorDouble CheckSigma2(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
    //                                                                   std::vector<CCMPMTKey> keys_to_fit,
    //                                                                  I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset);
};
#endif
