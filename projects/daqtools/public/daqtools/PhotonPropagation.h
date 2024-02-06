#ifndef PhotonPropagation_H
#define PhotonPropagation_H

#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

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

#include <icetray/ctpl.h>
#include <icetray/open.h>
#include <icetray/I3Frame.h>
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
#include "icetray/robust_statistics.h"
#include <dataclasses/geometry/CCMGeometry.h>

struct PhotonPropagationJob {
    std::atomic<bool> running = false;
    std::vector<std::vector<double>>* vector_of_vertices = nullptr;
    bool vertex_1275_flag = false;
    size_t vector_of_vertices_index = 0;
    size_t event_start_idx = 0;
    size_t event_end_idx = 0;
    std::shared_ptr<std::vector<std::vector<std::vector<double>>>> vector_of_vertices_binned_charges = nullptr;
    std::shared_ptr<std::vector<std::vector<std::vector<double>>>> vector_of_vertices_binned_charges_squared = nullptr;
};

struct PhotonPropagationResult {
    size_t event_start_idx = 0;
    size_t event_end_idx = 0;
    std::shared_ptr<std::vector<std::vector<std::vector<double>>>> vector_of_vertices_binned_charges = nullptr;
    std::shared_ptr<std::vector<std::vector<std::vector<double>>>> vector_of_vertices_binned_charges_squared = nullptr;
    bool done = false;
};

namespace MCLLH {
namespace detail {

// compute the sum using the Kahan summation algorithm
template <class InIt>
typename std::iterator_traits<InIt>::value_type accumulate(InIt begin, InIt end) {
	typedef typename std::iterator_traits<InIt>::value_type real;
	real sum = real(0);
	real running_error = real(0);
	real temp;
	real difference;

	for (; begin != end; ++begin) {
		difference = *begin;
		difference -= running_error;
		temp = sum;
		temp += difference;
		running_error = temp;
		running_error -= sum;
		running_error -= difference;
		sum = std::move(temp);
	}
	return sum;
};

// compute log(1+x) without losing precision for small values of x
template<typename T>
T LogOnePlusX(T x)
{
	if (x <= -1.0)
	{
		std::stringstream os;
		os << "Invalid input argument (" << x
		   << "); must be greater than -1.0";
		throw std::invalid_argument( os.str() );
	}

	if (fabs(x) > 1e-4)
	{
		// x is large enough that the obvious evaluation is OK
		return log(1.0 + x);
	}

	// Use Taylor approx. log(1 + x) = x - x^2/2 with error roughly x^3/3
	// Since |x| < 10^-4, |x|^3 < 10^-12, relative error less than 10^-8
	T x2 = x*x;
	T x3 = x2*x;
	T x4 = x3*x;
	return x-x2/2.0+x3/3.0-x4/4.0;
};

template<typename DataType>
DataType get_mu(std::vector<DataType> const & wi) {
    return accumulate(wi.begin(), wi.end());
};

template<typename DataType>
DataType get_sigma2(std::vector<DataType> const & wi) {
    std::vector<double> w2i(wi.size());
    std::transform(wi.begin(), wi.end(), w2i.begin(), [](double w)->double{return w*w;});
    return accumulate(w2i.begin(), w2i.end());
};

} // namespace detail

struct poissonLikelihood{
	template <typename T>
	T operator()(double dataCount, T const & lambda, T const & w2_sum) const{
		if(lambda==0)
			return(dataCount==0?0:-std::numeric_limits<T>::infinity());
		T sum(lambda);
		sum+=lgamma(dataCount+1);
		return(dataCount*log(lambda)-sum);
	}
};

struct gammaPriorPoissonLikelihood {
	template<typename T>
	T operator()(double k, T const & alpha, T const & beta) {
		std::vector<T> items(5);
		items[0] = alpha*log(beta);
		items[1] = lgamma(k+alpha);
		items[2] = -lgamma(k+1);
		items[3] = -(k+alpha)*detail::LogOnePlusX(beta);
		items[4] = -lgamma(alpha);
		return detail::accumulate(items.begin(), items.end());
	}
};

struct LMean {
    template<typename T>
    T operator()(double k, T const & w_sum, T const & w2_sum) const {
        if(w_sum <= 0 || w2_sum < 0) {
            return(k==0?0:-std::numeric_limits<T>::infinity());
        }

        if(w2_sum == 0) {
            return poissonLikelihood()(k, w_sum, w2_sum);
        }

        T zero(0);
        if(w_sum == zero) {
            if(k == 0) {
                return zero;
            }
            else {
                return T(-std::numeric_limits<double>::infinity());
            }
        }

        T alpha = w_sum*w_sum/w2_sum;
        T beta = w_sum/w2_sum;
        T L = gammaPriorPoissonLikelihood()(k, alpha, beta);

        return L;
    }
};

struct LMode {
	template<typename T>
	T operator()(double k, T const & w_sum, T const & w2_sum) const {
		if(w_sum <= 0 || w2_sum < 0) {
			return(k==0?0:-std::numeric_limits<T>::infinity());
		}

		if(w2_sum == 0) {
			return poissonLikelihood()(k, w_sum, w2_sum);
		}

		T zero(0);
		if(w_sum == zero) {
			if(k == 0) {
				return zero;
			}
			else {
				return T(-std::numeric_limits<double>::infinity());
			}
		}

		const T & mu = w_sum;
		T mu2 = mu*mu;
		const T & sigma2 = w2_sum;

		T beta = (mu + sqrt(mu2+sigma2*4.0))/(sigma2*2);
		T alpha = (mu*sqrt(mu2+sigma2*4.0)/sigma2 + mu2/sigma2 + 2.0) / 2.0;
		T L = gammaPriorPoissonLikelihood()(k, alpha, beta);

		return L;
	}
};

struct LEff {
    template<typename T>
    T operator()(double k, T const & w_sum, T const & w2_sum) const {
        if(w_sum <= 0 || w2_sum < 0) {
            return(k==0?0:-std::numeric_limits<T>::infinity());
        }

        if(w2_sum == 0) {
            return poissonLikelihood()(k, w_sum, w2_sum);
        }

        T zero(0);
        if(w_sum == zero) {
            if(k == 0) {
                return zero;
            }
            else {
                return T(-std::numeric_limits<double>::infinity());
            }
        }

        T alpha = w_sum*w_sum/w2_sum + 1.0;
        T beta = w_sum/w2_sum;
        T L = gammaPriorPoissonLikelihood()(k, alpha, beta);

        return L;
    }
};


struct computeLMean {
    template<typename T>
    T operator()(unsigned int k, const std::vector<T>& wi) const{
        T mu = detail::get_mu(wi);
        T sigma2 = detail::get_sigma2(wi);
        return LMean()(k, mu, sigma2);
    }
};

struct computeLMode {
    template<typename T>
    T operator()(unsigned int k, const std::vector<T>& wi) const{
        T mu = detail::get_mu(wi);
        T sigma2 = detail::get_sigma2(wi);
        return LMode()(k, mu, sigma2);
    }
};

struct computeLEff {
    template<typename T>
    T operator()(double k, double mu, double sigma2) const{
        //T mu = detail::get_mu(wi);
        //T sigma2 = detail::get_sigma2(wi);
        return LEff()(k, mu, sigma2);
    }
};

} // namespace MCLLH

class PhotonPropagation {
    std::exception_ptr teptr = nullptr;
    std::string geometry_name_ = std::string("CCMGeometry");
    double smearing_mu_ = 1.0;
    double smearing_sigma_ = 1.0;
    double smearing_xi_ = 0.5;
    size_t n_convolution_chunks_ = 200.0;
    double portion_light_reflected_by_tpb_ = 1.0;
    double visible_absorption_length_ = 2000.0;
    size_t face_chunks_counter = 0;
    size_t side_chunks_counter = 0;
    // default values but can change with Set functions!
    double desired_chunk_width_ = 20.0; // use 5 for finer binning
    double desired_chunk_height_ = 20.0; // use 5 for finer binning
    size_t n_events_to_simulate_ = 1000;

    // place to store relevant information about our pmts!!!
    // pmt_x_loc, pmt_y_loc, pmt_z_loc, facing direction_x, facing_direction_y, facing_direction_z, coating flag, pmt facing area, pmt side area
    std::vector<std::vector<double>> pmt_parsed_information_;

    // similar place to store relevant information about our secondary locations
    // for our secondary locations, I think we can probably also pre-compute the yield and travel time for 1 photon from each location to each pmt
    // x_loc, y_loc, z_loc, facing_direction_x, facing_direction_y, facing_direction_z, TPB portion, facing area, side area
    std::vector<std::vector<double>> locations_to_check_information_;
    std::vector<std::vector<double>> locations_to_check_to_pmt_yield_;
    std::vector<std::vector<double>> locations_to_check_to_pmt_travel_time_;

    // place to store list of vertices to simulate
    std::vector<std::vector<double>> verticies_to_simuate_1275_;
    std::vector<std::vector<double>> verticies_to_simuate_;
    std::vector<std::vector<double>> verticies_to_simuate_511_;
    unsigned int coated_omtype = (unsigned int)10;
    unsigned int uncoated_omtype = (unsigned int)20;

    // defining some geomtry things for modelling the detector...not the most elegant
    double pmt_radius = 10.16; //radius in cm^2
    double pmt_facing_area = M_PI * std::pow(pmt_radius, 2);
    double pmt_side_area_factor = 0.217549;
    double pmt_side_area = pmt_facing_area * pmt_side_area_factor;
    double cylinder_max_x = 96.0;
    double cylinder_min_x = - cylinder_max_x;
    double cylinder_max_y = 96.0;
    double cylinder_min_y = - cylinder_max_y;
    double cylinder_max_z = 58.0;
    double cylinder_min_z = - cylinder_max_z;
    double cylinder_radius = cylinder_max_x;
    double cylinder_circumference = M_PI * 2 * cylinder_max_x;
    double cylinder_height = cylinder_max_z * 2;
    double chunk_side_area_factor = 0.1; // this number is a guess..maybe model it one day

    // now some geometry things for throwing source events
    double rod_diameter = 1.0;
    double source_diameter = 0.8;
    double rod_width = (rod_diameter - source_diameter)/2;
    double source_inset = -0.25;
    double decay_constant = 10.0;
    double source_rod_lower_end_cap = - source_inset;
    double detector_lower_end_cap = cylinder_min_z;
    double detector_radius = cylinder_radius;
    double pos_rad = source_diameter / 2;
    size_t total_events_that_escaped = 0;
    double source_z_offset_ = 0.0;

    // pmt noise rate
    double noise_photons = 1.0;
    double noise_triggers = 5.0;
    double digitization_time = 16 * std::pow(10, 3); //16 usec in nsec
    double noise_rate = noise_photons / (noise_triggers * digitization_time); // units of photons/nsec
    double noise_rate_per_time_bin = 2.0 * noise_rate; // 2nsec binning

    // some constants we use for the simulation
    double c = 2.998 * std::pow(10, 8); // speed of light in m/s
    double c_cm_per_nsec = c * std::pow(10, -7); // speed of light in cm/nsec
    double uv_index_of_refraction = 1.358;
    double vis_index_of_refraction = 1.23;
    double pmt_quantum_efficiency = 0.25;
    double full_acceptance = 4.0 * M_PI;
    size_t n_pmts_to_simulate = (size_t) 0;

    ctpl::thread_pool pool;
    size_t num_threads;
    std::vector<std::vector<std::vector<double>>> thread_verticies_;
    std::vector<std::vector<std::vector<double>>> thread_1275_verticies_;
    std::vector<std::vector<std::vector<double>>> thread_511_verticies_;
    size_t max_cached_vertices = (size_t) 2000;

    //std::deque<PhotonPropagationJob *> free_jobs;
    //std::deque<PhotonPropagationJob *> running_jobs;
    //std::deque<PhotonPropagationResult> results;

    I3Vector<I3Vector<double>> data_series_;
    size_t n_data_samples_;
    double max_data_value_;
    double time_of_max_data_value_;
    std::vector<double> times_of_data_points_;

public:
    PhotonPropagation();
    void SetData(I3Vector<I3Vector<double>> data_series);
    void SetDataSampleSize(size_t n_data_samples);
    void SetNThreads(size_t const & n_threads);
    size_t GetNFaceChunks();
    size_t GetNSideChunks();
    size_t GetNSimulatedEvents();
    void SetZOffset(double source_z_offset);
    void GetEventVertices(size_t const & n_events_to_simulate);
    void GetPMTInformation(I3FramePtr frame);
    void GetSecondaryLocs(double const & desired_chunk_width, double const & desired_chunk_height);
    //I3Vector<I3Vector<I3Vector<double>>> GetSimulation(double const & singlet_ratio_,
    //I3Vector<I3Vector<double>> GetSimulation(double const & singlet_ratio_,
    double GetSimulation(double const & singlet_ratio_,
                         double const & triplet_ratio_,
                         double const & singlet_tau_,
                         double const & triplet_tau_,
                         double const & recombination_tau_,
                         double const & TPB_ratio_,
                         double const & TPB_tau_,
                         double const & UV_absorption_length_,
                         double const & n_photons_produced_);
};
#endif
