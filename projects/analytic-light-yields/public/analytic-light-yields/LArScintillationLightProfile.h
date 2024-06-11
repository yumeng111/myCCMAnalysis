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
#include "analytic-light-yields/autodiff.h"

template<typename T>
struct mpmath_gammainc {
    T operator() (double const & a, T const & x) const {
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
        T result = bp::extract<T>(py_float_result);
        return result;
    }
};

template <unsigned int nVars, typename T>
struct mpmath_gammainc<phys_tools::autodiff::FD<nVars, T>> {
    mpmath_gammainc<T> base;
    using AD = phys_tools::autodiff::FD<nVars, T>;

    AD operator() (double const & a, AD const & x) const {
        AD result(base(a, x.value()));
        T gradient[nVars];
        T dGammadX = - exp(-x.value()) * std::pow(x.value(), a-1);
        for(unsigned int i = 0; i < nVars; i++) {
            gradient[i] = dGammadX * x.derivative(i);
        }
        result.copyGradient(gradient);
        return result;
    }
};

template<typename T> T normal_distribution(T const & time, T const & mu, T const & sigma) {
    T scale = 1.0 / (sigma * std::sqrt(2.0 * M_PI));
    T z = (time - mu) / sigma;
    T e = exp(-0.5 * z * z);

    return scale * e;
}

template<typename T> T add_late_pulse_bump(T const & mu,
                                           T const & sigma,
                                           T const & scale,
                                           T const & max_light_prof,
                                           T const & time,
                                           T const & light_profile) {
    T light_prof_w_gauss = light_profile;
    T gaussian_peak_value = normal_distribution(mu, mu, sigma);
    if (time > 10.0){
        light_prof_w_gauss += (scale * max_light_prof * (normal_distribution(time, mu, sigma) / gaussian_peak_value));
    }
    return light_prof_w_gauss;
}

template<typename T>
void get_total_light_profile (T const & R_s,
                              T const & R_t,
                              T const & tau_s,
                              T const & tau_t,
                              T const & tau_rec,
                              T const & tau_TPB,
                              T const & late_pulse_mu, 
                              T const & late_pulse_sigma,
                              T const & late_pulse_scale, 
                              std::vector<T> const & times,
                              std::vector<T> & final_light_profile) {

    mpmath_gammainc<T> gammainc;

    // times is a vector of times to calculate the light profile for
    T fraction_recombination = (1.0 - R_s - R_t);
    T gamma_const = gammainc(-1, -tau_rec / tau_TPB);
    T coeff_one = R_s / (tau_s - tau_TPB);
    T coeff_two = R_t / (tau_t - tau_TPB);
    T coeff_three = fraction_recombination * tau_rec / (tau_TPB * tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        T const & t = times.at(time_it);
        if(t < 0) {
            final_light_profile.at(time_it) = 1e-18 * exp(t / 10.0);
            continue;
        }
        T exp_singlet = exp(-t / tau_s);
        T exp_triplet = exp(-t / tau_t);
        T exp_recombination = exp(-(t+tau_rec) / tau_TPB);
        T exp_prompt_TPB = exp(-t / tau_TPB);
        T gamma_t = gammainc(-1, -(t + tau_rec) / tau_TPB);

        T one = coeff_one * (exp_singlet - exp_prompt_TPB);
        T two = coeff_two * (exp_triplet - exp_prompt_TPB);
        T three = coeff_three * exp_recombination * (gamma_t - gamma_const);

        T y = one + two + three;

        final_light_profile.at(time_it) = y;
    }
}

template<typename T>
void get_light_profile_no_recombination(T const & R_s,
                                        T const & tau_s,
                                        T const & tau_t,
                                        T const & tau_TPB,
                                        T const & late_pulse_mu, 
                                        T const & late_pulse_sigma,
                                        T const & late_pulse_scale, 
                                        std::vector<T> const & times,
                                        std::vector<T> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    T R_t = 1.0 - R_s;
    T coeff_one = R_s / (tau_s - tau_TPB);
    T coeff_two = R_t / (tau_t - tau_TPB);
    // let's loop over times and calculate the light profile at each time
    // note -- we are going to calculate the light profile then normalize so that maximum is equal to 1.0
    // that way late pulse scale can be a percent of the max light profile
    std::vector<T> just_light_prof(final_light_profile.size());
    T max_light_prof;
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        T const & t = times.at(time_it);
        if(t <= 0) {
            just_light_prof.at(time_it) = 1e-18 * exp(t / 10.0);
            max_light_prof = just_light_prof.at(time_it);
            continue;
        }
        T exp_singlet = exp(-t / tau_s);
        T exp_triplet = exp(-t / tau_t);
        T exp_prompt_TPB = exp(-t / tau_TPB);

        T one = coeff_one * (exp_singlet - exp_prompt_TPB);
        T two = coeff_two * (exp_triplet - exp_prompt_TPB);

        T y = one + two;
        just_light_prof.at(time_it) = y; 
        if (y > max_light_prof){
            max_light_prof = y;
        }
    }

    // ok so we have just the light profile
    // now let's loop over times again and divide by the max value then add in our gaussian
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        T const & t = times.at(time_it);
        final_light_profile.at(time_it) = add_late_pulse_bump(late_pulse_mu, late_pulse_sigma, late_pulse_scale, max_light_prof, t, just_light_prof.at(time_it));
    }

}

#endif
