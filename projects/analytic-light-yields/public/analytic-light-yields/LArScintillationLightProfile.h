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
struct mpmath_gammainc<FD<nVars, T>> {
    mpmath_gammainc<T> base;
    using AD = FD<nVars, T>;

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

template<typename T>
void get_total_light_profile (T const & R_s,
                            T const & R_t,
                            T const & tau_s,
                            T const & tau_t,
                            T const & tau_rec,
                            T const & tau_TPB,
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
        T exp_singlet = std::exp(-t / tau_s);
        T exp_triplet = std::exp(-t / tau_t);
        T exp_recombination = std::exp(-(t+tau_rec) / tau_TPB);
        T exp_prompt_TPB = std::exp(-t / tau_TPB);
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
                                        T const & R_t,
                                        T const & tau_s,
                                        T const & tau_t,
                                        T const & tau_TPB,
                                        std::vector<T> const & times,
                                        std::vector<T> & final_light_profile) {

    // times is a vector of times to calculate the light profile for
    T coeff_one = R_s / (tau_s - tau_TPB);
    T coeff_two = R_t / (tau_t - tau_TPB);

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++) {
        T const & t = times.at(time_it);
        if(t < 0) {
            final_light_profile.at(time_it) = 1e-18 * exp(t / 10.0);
            continue;
        }
        T exp_singlet = std::exp(-t / tau_s);
        T exp_triplet = std::exp(-t / tau_t);
        T exp_prompt_TPB = std::exp(-t / tau_TPB);

        T one = coeff_one * (exp_singlet - exp_prompt_TPB);
        T two = coeff_two * (exp_triplet - exp_prompt_TPB);

        T y = one + two;

        final_light_profile.at(time_it) = y;
    }
}

#endif
