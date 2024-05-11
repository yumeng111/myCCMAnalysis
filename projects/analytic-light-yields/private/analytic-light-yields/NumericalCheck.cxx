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
#include "dataclasses/geometry/CCMGeometry.h"
#include "analytic-light-yields/CalculateNLLH.h"
#include "analytic-light-yields/NumericalCheck.h"
#include "analytic-light-yields/autodiff.h"

NumericalCheck::NumericalCheck() {}

void NumericalCheck::CheckSimplifiedLightProfile(AnalyticLightYieldGenerator analytic_light_yield_setup, double light_time_offset) {
    // given light parameters, print out times and values of our simplified light profile


    typedef phys_tools::autodiff::FD<6, double> AD;

    AD Rs(analytic_light_yield_setup.Rs, 0);
    AD Rt(analytic_light_yield_setup.Rt, 1);
    AD tau_s(analytic_light_yield_setup.tau_s, 2);
    AD tau_t(analytic_light_yield_setup.tau_t, 3);
    AD tau_TPB(analytic_light_yield_setup.tau_TPB, 4);
    AD time_offset(light_time_offset, 5);

    std::cout << "getting light profile for Rs = " << Rs.value() << ", Rt = " << Rt.value() << ", tau_s = " << tau_s.value() << ", tau_t = " << tau_t.value() << ", and tau_TPB  = " << tau_TPB.value() << std::endl;

    size_t n_bins = 250;
    std::vector<AD> times(n_bins);
    for (size_t i = 0; i < n_bins; i++) {
        times[i] = i * 2.0 + light_time_offset;
    }

    I3Vector<AD> light_profile(n_bins);
    get_light_profile_no_recombination(Rs, Rt, tau_s, tau_t, tau_TPB, times, light_profile);

    // copying over time binning from LArScintillationLightProfile
    for (size_t i = 0; i < n_bins; i++){
        std::cout << "at time = " <<  (double) i * 2.0 << " light profile = " << light_profile.at(i).value() << std::endl;
    }

    std::cout << "and now for derivatives : " << std::endl;
    for (size_t i = 0; i < n_bins; i++){
        std::cout << "at time = " <<  (double) i * 2.0 << " dRs = " << light_profile.at(i).derivative(0)
                                                       << ", dRt = " << light_profile.at(i).derivative(1)
                                                       << ", dtau_s = " << light_profile.at(i).derivative(2)
                                                       << ", dtau_t = " << light_profile.at(i).derivative(3)
                                                       << ", and dtau_TPB = " << light_profile.at(i).derivative(4) << std::endl;
    }
}

void NumericalCheck::CheckNLLHDerivs(double const & k, double const & mu, double const & sigma_squared, double const & DmuDtheta, double const & Dsigma_squaredDtheta){

    double nllh = MCLLH::LEff()(k, mu, sigma_squared);
    double deriv = MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
    std::cout << "for k = " << k << ", mu = " << mu << ", sigma_squared = " << sigma_squared
        << ", DmuDtheta = " << DmuDtheta << ", and Dsigma_squaredDtheta = " << Dsigma_squaredDtheta << " : " << std::endl;

    std::cout << "nllh = " << nllh << " and gradient = " << deriv << std::endl;
}

//I3MapPMTKeyVectorDouble NumericalCheck::CheckData(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
//                                                                std::vector<CCMPMTKey> keys_to_fit,
//                                                               I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset){
//
//    if (CalculateNLLH_constructor == nullptr){
//        CalculateNLLH_constructor = std::make_shared<CalculateNLLH> ();
//    }
//    // let's grab our nllh
//    dubug_info = CalculateNLLH_constructor->DebugDatavsPred(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
//    return dubug_info[0];
//}
//
//I3MapPMTKeyVectorDouble NumericalCheck::CheckPred(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
//                                                                std::vector<CCMPMTKey> keys_to_fit,
//                                                               I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset){
//
//    if (CalculateNLLH_constructor == nullptr){
//        CalculateNLLH_constructor = std::make_shared<CalculateNLLH> ();
//        dubug_info = CalculateNLLH_constructor->DebugDatavsPred(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
//    }
//
//    return dubug_info[1];
//}
//
//
//I3MapPMTKeyVectorDouble NumericalCheck::CheckSigma2(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
//                                                                std::vector<CCMPMTKey> keys_to_fit,
//                                                               I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset){
//    if (CalculateNLLH_constructor == nullptr){
//        CalculateNLLH_constructor = std::make_shared<CalculateNLLH> ();
//        dubug_info = CalculateNLLH_constructor->DebugDatavsPred(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
//    }
//
//    return dubug_info[2];
//}

