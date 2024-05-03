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
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <analytic-light-yields/LArScintillationLightProfile.h>

double gamma_m1(double const & x){

    double E1 = -boost::math::expint(-x);
    return std::pow(x, -1.0) * std::exp(-x) - E1;

}

void get_total_light_profile(double const & R_s,
                            double const & R_t,
                            double const & tau_s,
                            double const & tau_t,
                            double const & tau_rec,
                            double const & tau_TPB,
                            std::vector<double> const & times,
                            std::vector<double> & final_light_profile){

    // times is a vector of times to calculate the light profile for
    double fraction_recombination = (1.0 - R_s - R_t);
    double exp_singlet;
    double exp_triplet;
    double exp_recombination;
    double exp_prompt_TPB;
    double gamma_const;
    double gamma_t;
    double log_t;
    double one;
    double two;
    double three;
    double t;

    // let's loop over times and calculate the light profile at each time
    for (size_t time_it = 0; time_it < times.size(); time_it++){

        t = times.at(time_it);
        exp_singlet = std::exp(-t / tau_s);
        exp_triplet = std::exp(-t / tau_t);
        exp_recombination = std::exp(-(t+tau_rec) / tau_TPB);
        exp_prompt_TPB = std::exp(-t / tau_TPB);

        gamma_const = gamma_m1(tau_rec / tau_TPB);
        gamma_t = gamma_m1((t + tau_rec)/tau_TPB);
        log_t = 2.0 * std::log(t / tau_rec + 1.0);

        one = R_s * (exp_singlet - exp_prompt_TPB) / (tau_s - tau_TPB);
        two = R_t * (exp_triplet - exp_prompt_TPB) / (tau_t - tau_TPB);
        three = fraction_recombination * exp_recombination * tau_rec * (gamma_t - gamma_const + log_t) / (tau_TPB * tau_TPB);

        final_light_profile.at(time_it) = one + two + three;
    }

}


I3Vector<double> LArScintillationLightProfile::GetLightProfile(double const & singlet_ratio,
                                                               double const & triplet_ratio,
                                                               double const & singlet_tau,
                                                               double const & triplet_tau,
                                                               double const & recombination_tau,
                                                               double const & TPB_tau) {

    // let's set up our light profile
    double start_time = 0.0;
    double end_time = 1200.0;
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

