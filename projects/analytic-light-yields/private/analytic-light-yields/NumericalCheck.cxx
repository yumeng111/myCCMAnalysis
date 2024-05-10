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
#include <analytic-light-yields/NumericalCheck.h>

NumericalCheck::NumericalCheck() {}

void NumericalCheck::CheckSimplifiedLightProfile(AnalyticLightYieldGenerator analytic_light_yield_setup, double light_time_offset){
    // given light parameters, print out times and values of our simplified light profile
    if (LAr_scintillation_light_constructor == nullptr){
        LAr_scintillation_light_constructor = std::make_shared<LArScintillationLightProfile> ();
    }
    
    double Rs = analytic_light_yield_setup.Rs;
    double Rt = analytic_light_yield_setup.Rt;
    double tau_s = analytic_light_yield_setup.tau_s;
    double tau_t = analytic_light_yield_setup.tau_t;
    double tau_TPB = analytic_light_yield_setup.tau_TPB;

    std::cout << "getting light profile for Rs = " << Rs << ", Rt = " << Rt << ", tau_s = " << tau_s << ", tau_t = " << tau_t << ", and tau_TPB  = " << tau_TPB << std::endl;

    I3Vector<double> light_profile;
    light_profile = LAr_scintillation_light_constructor->GetSimplifiedLightProfile(Rs, Rt, tau_s, tau_t, tau_TPB, light_time_offset);
   
    // copying over time binning from LArScintillationLightProfile 
    double start_time = 0.0;
    double end_time = 100.0;
    double bin_width = 2.0;
    size_t n_time_bins = (size_t)((end_time - start_time)/bin_width);

    for (size_t i = 0; i < n_time_bins; i++){
        std::cout << "at time = " <<  (double) i * 2.0 << " light profile = " << light_profile.at(i) << std::endl;
    }

    // now let's also check our derivatives 
    I3Vector<double> light_profile_derivRs;
    light_profile_derivRs = LAr_scintillation_light_constructor->GetSimplifiedLightProfileDeriv(Rs, Rt, tau_s, tau_t, tau_TPB, "Rs", light_time_offset);
    I3Vector<double> light_profile_derivRt;
    light_profile_derivRt = LAr_scintillation_light_constructor->GetSimplifiedLightProfileDeriv(Rs, Rt, tau_s, tau_t, tau_TPB, "Rt", light_time_offset);
    I3Vector<double> light_profile_derivtau_s;
    light_profile_derivtau_s = LAr_scintillation_light_constructor->GetSimplifiedLightProfileDeriv(Rs, Rt, tau_s, tau_t, tau_TPB, "tau_s", light_time_offset);
    I3Vector<double> light_profile_derivtau_t;
    light_profile_derivtau_t = LAr_scintillation_light_constructor->GetSimplifiedLightProfileDeriv(Rs, Rt, tau_s, tau_t, tau_TPB, "tau_t", light_time_offset);
    I3Vector<double> light_profile_derivtau_TPB;
    light_profile_derivtau_TPB = LAr_scintillation_light_constructor->GetSimplifiedLightProfileDeriv(Rs, Rt, tau_s, tau_t, tau_TPB, "tau_TPB", light_time_offset);

    std::cout << "and now for derivatives : " << std::endl;
    for (size_t i = 0; i < n_time_bins; i++){
        std::cout << "at time = " <<  (double) i * 2.0 << " dRs = " << light_profile_derivRs.at(i) 
                                                       << ", dRt = " << light_profile_derivRt.at(i) 
                                                       << ", dtau_s = " << light_profile_derivtau_s.at(i)
                                                       << ", dtau_t = " << light_profile_derivtau_t.at(i)
                                                       << ", and dtau_TPB = " << light_profile_derivtau_TPB.at(i) << std::endl;
    }

}

void NumericalCheck::CheckNLLHDerivs(double const & k, double const & mu, double const & sigma_squared, double const & DmuDtheta, double const & Dsigma_squaredDtheta){

    double nllh = MCLLH::LEff()(k, mu, sigma_squared);
    double deriv = MCLLH::LEffDeriv()(k, mu, sigma_squared, DmuDtheta, Dsigma_squaredDtheta);
    std::cout << "for k = " << k << ", mu = " << mu << ", sigma_squared = " << sigma_squared
        << ", DmuDtheta = " << DmuDtheta << ", and Dsigma_squaredDtheta = " << Dsigma_squaredDtheta << " : " << std::endl;

    std::cout << "nllh = " << nllh << " and gradient = " << deriv << std::endl;
}

I3MapPMTKeyVectorDouble NumericalCheck::CheckData(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, 
                                                                std::vector<CCMPMTKey> keys_to_fit,
                                                               I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset){

    if (CalculateNLLH_constructor == nullptr){
        CalculateNLLH_constructor = std::make_shared<CalculateNLLH> ();
    }
    // let's grab our nllh
    dubug_info = CalculateNLLH_constructor->DebugDatavsPred(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
    return dubug_info[0]; 
}

I3MapPMTKeyVectorDouble NumericalCheck::CheckPred(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, 
                                                                std::vector<CCMPMTKey> keys_to_fit,
                                                               I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset){

    if (CalculateNLLH_constructor == nullptr){
        CalculateNLLH_constructor = std::make_shared<CalculateNLLH> ();
        dubug_info = CalculateNLLH_constructor->DebugDatavsPred(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
    }
    
    return dubug_info[1]; 
}


I3MapPMTKeyVectorDouble NumericalCheck::CheckSigma2(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame, 
                                                                std::vector<CCMPMTKey> keys_to_fit,
                                                               I3FramePtr data_frame, double single_pmt_norm, double single_pmt_offset, double light_time_offset){
    if (CalculateNLLH_constructor == nullptr){
        CalculateNLLH_constructor = std::make_shared<CalculateNLLH> ();
        dubug_info = CalculateNLLH_constructor->DebugDatavsPred(analytic_light_yield_setup, geo_frame, data_frame, keys_to_fit, single_pmt_norm, single_pmt_offset, light_time_offset);
    }
    
    return dubug_info[2]; 
}

