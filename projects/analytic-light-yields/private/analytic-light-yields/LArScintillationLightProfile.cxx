#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/python.hpp>

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
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>

double mpmath_gammainc(double const & a, double const & x) {
    namespace bp = boost::python;
    bp::object builtins = bp::import("builtins");
    bp::object mpmath_embed;
    try {
        mpmath_embed = bp::import("mpmath");
    } catch( bp::error_already_set& e ) {
        log_error("Unable to import mpmath");
        PyErr_Clear();
        return -1.0;
    }
    bp::object gammainc_function = mpmath_embed.attr("gammainc");
    bp::object py_result = gammainc_function(a, x);
    bp::object py_real_result = py_result.attr("real");
    bp::object py_float_result = builtins.attr("float")(py_real_result);
    double result = bp::extract<double>(py_float_result);
    return result;
}

void get_total_light_profile(double const & R_s,
                            double const & R_t,
                            double const & tau_s,
                            double const & tau_t,
                            double const & tau_rec,
                            double const & tau_TPB,
                            std::vector<double> const & times,
                            std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double fraction_recombination = (1.0 - R_s - R_t);
    double gamma_const = mpmath_gammainc(-1, -tau_rec / tau_TPB);
    double coeff_one = R_s / (tau_s - tau_TPB);
    double coeff_two = R_t / (tau_t - tau_TPB);
    double coeff_three = fraction_recombination * tau_rec / (tau_TPB * tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_singlet = std::exp(-t / tau_s);
        double exp_triplet = std::exp(-t / tau_t);
        double exp_recombination = std::exp(-(t+tau_rec) / tau_TPB);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);
        double gamma_t = mpmath_gammainc(-1, -(t + tau_rec) / tau_TPB);

        double one = coeff_one * (exp_singlet - exp_prompt_TPB);
        double two = coeff_two * (exp_triplet - exp_prompt_TPB);
        double three = coeff_three * exp_recombination * (gamma_t - gamma_const);

        double y = one + two + three;

        final_light_profile.at(time_it) = y;
    }
}

void get_light_profile_no_recombination(double const & R_s,
                                        double const & R_t,
                                        double const & tau_s,
                                        double const & tau_t,
                                        double const & tau_TPB,
                                        std::vector<double> const & times,
                                        std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = R_s / (tau_s - tau_TPB);
    double coeff_two = R_t / (tau_t - tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_singlet = std::exp(-t / tau_s);
        double exp_triplet = std::exp(-t / tau_t);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * (exp_singlet - exp_prompt_TPB);
        double two = coeff_two * (exp_triplet - exp_prompt_TPB);

        double y = one + two;

        final_light_profile.at(time_it) = y;
    }
}

void D_light_profile_no_recombination_dtau_s(double const & R_s,
                                            double const & R_t,
                                            double const & tau_s,
                                            double const & tau_t,
                                            double const & tau_TPB,
                                            std::vector<double> const & times,
                                            std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = -R_s / std::pow((tau_s - tau_TPB), 2.0);
    double coeff_two = R_s / ((tau_s - tau_TPB) * std::pow(tau_s, 2.0));

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_singlet = std::exp(-t / tau_s);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * (exp_singlet - exp_prompt_TPB);
        double two = coeff_two * t * (exp_singlet);

        double y = one + two;

        final_light_profile.at(time_it) = y;
    }
}

void D_light_profile_no_recombination_dtau_t(double const & R_s,
                                            double const & R_t,
                                            double const & tau_s,
                                            double const & tau_t,
                                            double const & tau_TPB,
                                            std::vector<double> const & times,
                                            std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = -R_t / std::pow((tau_t - tau_TPB), 2.0);
    double coeff_two = R_t / ((tau_t - tau_TPB) * std::pow(tau_t, 2.0));

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_triplet = std::exp(-t / tau_t);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * (exp_triplet - exp_prompt_TPB);
        double two = coeff_two * t * (exp_triplet);

        double y = one + two;

        final_light_profile.at(time_it) = y;
    }
}

void D_light_profile_no_recombination_dtau_TPB(double const & R_s,
                                            double const & R_t,
                                            double const & tau_s,
                                            double const & tau_t,
                                            double const & tau_TPB,
                                            std::vector<double> const & times,
                                            std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = R_s / std::pow((tau_s - tau_TPB), 2.0);
    double coeff_two = R_t / std::pow((tau_t- tau_TPB), 2.0);
    double coeff_three = -R_s / ((tau_s - tau_TPB) * std::pow(tau_TPB, 2.0));
    double coeff_four = -R_t / ((tau_t - tau_TPB) * std::pow(tau_TPB, 2.0));

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_singlet = std::exp(-t / tau_s);
        double exp_triplet = std::exp(-t / tau_t);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * (exp_singlet - exp_prompt_TPB);
        double two = coeff_two * (exp_triplet - exp_prompt_TPB);
        double three = coeff_three * t * (exp_prompt_TPB);
        double four = coeff_four * t * (exp_prompt_TPB);

        double y = one + two + three + four;

        final_light_profile.at(time_it) = y;
    }
}

void D_light_profile_no_recombination_dR_s(double const & R_s,
                                            double const & R_t,
                                            double const & tau_s,
                                            double const & tau_t,
                                            double const & tau_TPB,
                                            std::vector<double> const & times,
                                            std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = 1.0 / (tau_s - tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_singlet = std::exp(-t / tau_s);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * (exp_singlet - exp_prompt_TPB);

        double y = one;

        final_light_profile.at(time_it) = y;
    }
}

void D_light_profile_no_recombination_dR_t(double const & R_s,
                                            double const & R_t,
                                            double const & tau_s,
                                            double const & tau_t,
                                            double const & tau_TPB,
                                            std::vector<double> const & times,
                                            std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = 1.0 / (tau_t - tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_triplet = std::exp(-t / tau_t);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * (exp_triplet - exp_prompt_TPB);

        double y = one;

        final_light_profile.at(time_it) = y;
    }
}

void D_light_profile_no_recombination_dtoffset(double const & R_s,
                                               double const & R_t,
                                               double const & tau_s,
                                               double const & tau_t,
                                               double const & tau_TPB,
                                               std::vector<double> const & times,
                                               std::vector<double> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    double coeff_one = R_s / (tau_s - tau_TPB);
    double coeff_two = R_t / (tau_t - tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        double const & t = times.at(time_it);
        double exp_singlet = std::exp(-t / tau_s);
        double exp_triplet = std::exp(-t / tau_t);
        double exp_prompt_TPB = std::exp(-t / tau_TPB);

        double one = coeff_one * ((exp_singlet / tau_s) - (exp_prompt_TPB / tau_TPB));
        double two = coeff_two * ((exp_triplet / tau_t) - (exp_prompt_TPB / tau_TPB));

        double y = one + two;

        final_light_profile.at(time_it) = y;
    }
}

LArScintillationLightProfile::LArScintillationLightProfile() {}

I3Vector<double> LArScintillationLightProfile::GetFullLightProfile(double const & singlet_ratio,
                                                                   double const & triplet_ratio,
                                                                   double const & singlet_tau,
                                                                   double const & triplet_tau,
                                                                   double const & recombination_tau,
                                                                   double const & TPB_tau) {

    // let's set up our light profile
    double start_time = 0.0;
    double end_time = 100.0;
    double bin_width = 2.0;
    size_t n_time_bins = (size_t)((end_time - start_time)/bin_width);

    std::vector<double> times(n_time_bins);

    for (size_t i = 0; i < n_time_bins; i++){
        times.at(i) = (double) i * 2.0;
    }

    // now let's get out light profile
    I3Vector<double> light_profile(n_time_bins);
    get_total_light_profile(singlet_ratio,
                            triplet_ratio,
                            singlet_tau,
                            triplet_tau,
                            recombination_tau,
                            TPB_tau,
                            times,
                            light_profile);

    return light_profile;
}

I3Vector<double> LArScintillationLightProfile::GetSimplifiedLightProfile(double const & singlet_ratio,
                                                                         double const & triplet_ratio,
                                                                         double const & singlet_tau,
                                                                         double const & triplet_tau,
                                                                         double const & TPB_tau) {

    // let's set up our light profile
    double start_time = 0.0;
    double end_time = 100.0;
    double bin_width = 2.0;
    size_t n_time_bins = (size_t)((end_time - start_time)/bin_width);

    std::vector<double> times(n_time_bins);

    for (size_t i = 0; i < n_time_bins; i++){
        times.at(i) = (double) i * 2.0;
    }

    // now let's get out light profile
    I3Vector<double> light_profile(n_time_bins);
    get_light_profile_no_recombination(singlet_ratio,
                                       triplet_ratio,
                                       singlet_tau,
                                       triplet_tau,
                                       TPB_tau,
                                       times,
                                       light_profile);

    return light_profile;
}

I3Vector<double> LArScintillationLightProfile::GetSimplifiedLightProfileDeriv(double const & singlet_ratio,
                                                                              double const & triplet_ratio,
                                                                              double const & singlet_tau,
                                                                              double const & triplet_tau,
                                                                              double const & TPB_tau,
                                                                              std::string deriv_variable) {

    // let's set up our light profile
    double start_time = 0.0;
    double end_time = 100.0;
    double bin_width = 2.0;
    size_t n_time_bins = (size_t)((end_time - start_time)/bin_width);

    std::vector<double> times(n_time_bins);

    for (size_t i = 0; i < n_time_bins; i++){
        times.at(i) = (double) i * 2.0;
    }

    // now let's get out light profile
    I3Vector<double> light_profile(n_time_bins);
    if (deriv_variable == "Rs"){
        D_light_profile_no_recombination_dR_s(singlet_ratio, triplet_ratio, singlet_tau, triplet_tau, TPB_tau, times, light_profile);
    }
    else if (deriv_variable == "Rt"){
        D_light_profile_no_recombination_dR_t(singlet_ratio, triplet_ratio, singlet_tau, triplet_tau, TPB_tau, times, light_profile);
    }
    else if (deriv_variable == "tau_s"){
        D_light_profile_no_recombination_dtau_s(singlet_ratio, triplet_ratio, singlet_tau, triplet_tau, TPB_tau, times, light_profile);
    }
    else if (deriv_variable == "tau_t"){
        D_light_profile_no_recombination_dtau_t(singlet_ratio, triplet_ratio, singlet_tau, triplet_tau, TPB_tau, times, light_profile);
    }
    else if (deriv_variable == "tau_TPB"){
        D_light_profile_no_recombination_dtau_TPB(singlet_ratio, triplet_ratio, singlet_tau, triplet_tau, TPB_tau, times, light_profile);
    }
    else if (deriv_variable == "t_offset"){
        D_light_profile_no_recombination_dtoffset(singlet_ratio, triplet_ratio, singlet_tau, triplet_tau, TPB_tau, times, light_profile);
    }

    return light_profile;
}


