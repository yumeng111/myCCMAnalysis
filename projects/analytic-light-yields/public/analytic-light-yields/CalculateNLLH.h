#ifndef CalculateNLLH_H
#define CalculateNLLH_H

#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/polygamma.hpp>

#include <tuple>
#include <vector>
#include <random>
#include <cmath>
#include <map>
#include <memory>
#include <type_traits>
#include <algorithm>

#include "icetray/I3Units.h"
#include "dataclasses/I3Position.h"
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Double.h"
#include "analytic-light-yields/GenerateExpectation.h"
#include "analytic-light-yields/autodiff.h"

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

struct gammaPriorPoissonLikelihoodDerivative {
	template<typename T>
	T operator()(double k, T const & w, T const & w2, T const & DwDTheta, T const & Dw2DTheta) {
		std::vector<T> items(5);
        items[0] = (((1 + ((w*w) / w2)) / w) + (2 * w * log(w/w2) / w2)) * DwDTheta + (-((1 + ((w*w) / w2)) / w2) - (w*w*log(w/w2) / (w2*w2))) * Dw2DTheta;
		items[1] = ((2 * w * boost::math::polygamma(0, 1 + k + (w * w) / w2)) / w2) * DwDTheta - ((w * w * boost::math::polygamma(0, 1 + k + (w * w) / w2)) / (w2 * w2)) * Dw2DTheta;
		items[2] = 0;
		items[3] = (((-1 - k - ((w * w) / w2)) / (w2 * (1 + (w / w2)))) - ((2 * w * detail::LogOnePlusX(w / w2)) / w2)) * DwDTheta +
                    (-((w * (-1 - k - ((w * w) / w2))) / (w2 * w2 * (1 + (w / w2)))) + ((w * w * detail::LogOnePlusX(w / w2)) / (w2 * w2))) * Dw2DTheta;
		items[4] = (- (2 * w * boost::math::polygamma(0, 1 + ((w * w) / w2))) / w2) * DwDTheta + ((w * w * boost::math::polygamma(0, 1 + ((w * w) / w2))) / (w2 * w2)) * Dw2DTheta;
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
            //std::cout << "w/w2 < 0" << std::endl;
            //std::cout << "w = " << w_sum << ", w2 = " << w2_sum << std::endl;
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
                //std::cout << "w = 0" << std::endl;
                return T(-std::numeric_limits<double>::infinity());
            }
        }

        T alpha = w_sum*w_sum/w2_sum + 1.0;
        T beta = w_sum/w2_sum;
        T L = gammaPriorPoissonLikelihood()(k, alpha, beta);

        return L;
    }
};

struct LEffDeriv {
    template<typename T>
    T operator()(double k, T const & w_sum, T const & w2_sum, T const & Dw_sumDTheta, T const & Dw2_sumDTheta) const {
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

        T L = gammaPriorPoissonLikelihoodDerivative()(k, w_sum, w2_sum, Dw_sumDTheta, Dw2_sumDTheta);

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

struct SinglePMTInfo {
    CCMPMTKey key;
    std::vector<double> data;
    double peak_time;
    double peak_value;
    double start_time;
    double max_time;
    size_t event_start_bin;
    double pre_event_average;
    std::vector<double> data_times;
};

class CalculateNLLH {
public:
    static constexpr size_t n_params = 12;
    typedef double Underlying;

    typedef phys_tools::autodiff::FD<n_params, Underlying> AD;
    typedef std::array<Underlying, n_params> Grad;

    
    static constexpr size_t n_params_cppmin = 18;
    typedef double Underlying_cppmin;

    typedef phys_tools::autodiff::FD<n_params_cppmin, Underlying_cppmin> AD_cppmin;
    typedef std::array<Underlying_cppmin, n_params_cppmin> Grad_cppmin;

    I3MapPMTKeyVectorDouble debug_all_data;
    I3MapPMTKeyVectorDouble debug_all_pred;
    I3MapPMTKeyVectorDouble debug_fit_data;
    I3MapPMTKeyVectorDouble debug_fit_pred;
    std::vector<double> debug_pred_times;
    std::vector<double> debug_data_times;

    double event_start_threshold = 10.0;

private:
    std::map<CCMPMTKey, SinglePMTInfo> data;
    std::vector<CCMPMTKey> keys_to_fit;
    I3FramePtr geo_frame;

    size_t n_sodium_events;
    double portion_light_reflected_by_tpb;
    double desired_chunk_width;
    double desired_chunk_height;

    boost::shared_ptr<GenerateExpectation> gen_expectation = nullptr;
    double n_data_events;

    size_t max_bins;

public:
    CalculateNLLH();
    CalculateNLLH(I3FramePtr data_frame, I3FramePtr geo_frame, size_t max_bins=100, size_t n_sodium_events=20, double portion_light_reflected_by_tpb=1.0, double desired_chunk_width=20, double desired_chunk_height=20, std::vector<CCMPMTKey> keys_to_fit=std::vector<CCMPMTKey>());
    void SetKeys(std::vector<CCMPMTKey> keys);
    void SetGeo(I3FramePtr geo_frame);
    void SetData(I3FramePtr data_frame);
    void SetGenExpectation(boost::shared_ptr<GenerateExpectation> gen_expectation);

    double GetNDataEvents() const { return n_data_events; };
    I3MapPMTKeyVectorDoublePtr GetData() const;
    boost::shared_ptr<GenerateExpectation> GetGenExpectation() const { return gen_expectation; };

    I3MapPMTKeyVectorDouble GetAllDataForDebug() {return debug_all_data;};
    I3MapPMTKeyVectorDouble GetAllPredForDebug() {return debug_all_pred;};
    std::vector<double> GetPredTimesForDebug() {return debug_pred_times;};
    std::vector<double> GetDataTimesForDebug() {return debug_data_times;};
    I3MapPMTKeyVectorDouble GetFitDataForDebug() {return debug_fit_data;};
    I3MapPMTKeyVectorDouble GetFitPredForDebug() {return debug_fit_pred;};

    std::vector<double> GetLightProfileDebug() {return gen_expectation->GetLightProfileDebug();};
    std::vector<double> GetLightProfileTimesDebug() {return gen_expectation->GetLightProfileTimesDebug();};
    std::vector<double> GetLightProfileTOffsetGradDebug() {return gen_expectation->GetLightProfileTOffsetGradDebug();};

    template<typename T>
    T ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, T const_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale, 
            double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type);

    AD GetNLLH(CCMPMTKey key, AnalyticLightYieldGenerator const & params, const double & late_pulse_mu, const double & late_pulse_sigma, const double & late_pulse_scale);
    double GetNLLHValue(CCMPMTKey key, AnalyticLightYieldGenerator const & params, const double & late_pulse_mu, const double & late_pulse_sigma, const double & late_pulse_scale);
    std::vector<double> GetNLLHDerivative(CCMPMTKey key, AnalyticLightYieldGenerator const & params, const double & late_pulse_mu, const double & late_pulse_sigma, const double & late_pulse_scale);

    template<typename T> T Interpolation(double this_time, std::vector<T> data_times, std::vector<T> data);
    template<typename T> T LatePulseGaussian(double this_time, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale);
};

template<typename T> T CalculateNLLH::Interpolation(double this_time, std::vector<T> data_times, std::vector<T> data){
        
    // let's find closest times in data_times
    size_t closest_data_idx = 0;
    T data_time_diff = abs(data_times.at(0) - this_time);
    for (size_t data_it = 1; data_it < data_times.size(); data_it++){
        if (abs(data_times.at(data_it) - this_time) < data_time_diff){
            closest_data_idx = data_it;
            data_time_diff = abs(data_times.at(data_it) - this_time);
        } 
    }
    T closest_data_point;
    if (this_time == data_times.at(closest_data_idx)){
        // this is the case where our llh grid point overlaps with a data grid point
        closest_data_point = data.at(closest_data_idx);
    }
    else{
        // case where our llh grid point is not in data_times! 
        // we need to interpolate!
        T closest_data_time_below;
        T closest_data_value_below;
        T closest_data_time_above;
        T closest_data_value_above;

        if (data_times.at(closest_data_idx) < this_time){
            closest_data_time_below = data_times.at(closest_data_idx);
            closest_data_value_below = data.at(closest_data_idx);

            closest_data_time_above = data_times.at(closest_data_idx + 1);
            closest_data_value_above = data.at(closest_data_idx + 1);
        }
        else {
            closest_data_time_below = data_times.at(closest_data_idx - 1);
            closest_data_value_below = data.at(closest_data_idx - 1);

            closest_data_time_above = data_times.at(closest_data_idx);
            closest_data_value_above = data.at(closest_data_idx);
        }

        // now interpolate!!!
        closest_data_point = closest_data_value_below + (this_time - closest_data_time_below) * ((closest_data_value_above - closest_data_value_below)/ (closest_data_time_above - closest_data_time_below));
    }

    return closest_data_point;

}

template<typename T> T CalculateNLLH::LatePulseGaussian(double this_time, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale){
    return late_pulse_scale * (1 / (late_pulse_sigma * std::sqrt(2.0 * M_PI))) * exp(-0.5 * ((this_time - late_pulse_mu) * (this_time - late_pulse_mu)) / (late_pulse_sigma * late_pulse_sigma));   
}

template<typename T>
T CalculateNLLH::ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, T const_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale, 
            double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type) {
    // let's grab our data
    SinglePMTInfo const & pmt_data = data[key];
    double start_time = pmt_data.start_time;
    double max_time = pmt_data.max_time;
    double peak_time = pmt_data.peak_time;
    size_t event_start_bin = pmt_data.event_start_bin;
    double pre_event_average = pmt_data.pre_event_average;
    double data_peak_value = pmt_data.peak_value;
    std::vector<double> data_times = pmt_data.data_times;

    // let's grab our expectation
    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> pred = gen_expectation->GetExpectation(key, start_time, max_time,
            peak_time, Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_time_offset, uv_absorption, z_offset, n_sodium_events, light_profile_type);
    
    // unpack into yields, yields^2, and times
    boost::shared_ptr<std::vector<T>> pred_yields;
    boost::shared_ptr<std::vector<T>> pred_yields_squared;
    boost::shared_ptr<std::vector<T>> pred_times;
    pred_yields = std::get<0>(pred);
    pred_yields_squared = std::get<1>(pred);
    pred_times = std::get<2>(pred);
    T minimum_pred_time = pred_times->at(0);

    T total_nllh(0.0);

    size_t n_bins = pmt_data.data.size();

    // let's make a llh grid -- this combines both the times we have data values at along with the times we have prediction values at
    double starting_llh_grid_time = data_times.front();
    double ending_llh_grid_time = data_times.back();
    
    // to make things easier, let's make sure we have a vector of doubles containing pred times
    std::vector<double> pred_times_double;
    for (size_t j = 0; j < pred_times->size(); j++){
        if constexpr (std::is_same<T, double>::value) {
            pred_times_double.push_back(pred_times->at(j) + light_time_offset);
            //pred_times_double.push_back(pred_times->at(j));
        } else {
            pred_times_double.push_back(pred_times->at(j).value() + light_time_offset.value());
        }
    }


    std::vector<double> llh_grid;
    for (size_t i = 0; i < data_times.size(); i++){
        llh_grid.push_back(data_times.at(i));
    }
    for (size_t i = 0; i < pred_times_double.size(); i++){
        if (pred_times_double.at(i) >= starting_llh_grid_time and pred_times_double.at(i) <= ending_llh_grid_time){
            llh_grid.push_back(pred_times_double.at(i));
        }
    }

    std::vector<T> pred_times_offset;
    for (size_t j = 0; j < pred_times->size(); j++){
        pred_times_offset.push_back(pred_times->at(j) + light_time_offset);   
    }
    
    // now sort our llh grid!
    std::sort(llh_grid.begin(), llh_grid.end());

    // let's add empty vector to debug_data and debug_pred
    debug_all_data[key] = std::vector<double> (llh_grid.size(), 0.0);
    debug_all_pred[key] = std::vector<double> (llh_grid.size(), 0.0);
    debug_fit_data[key] = std::vector<double> (llh_grid.size(), 0.0);
    debug_fit_pred[key] = std::vector<double> (llh_grid.size(), 0.0);
    debug_pred_times.clear();
    debug_data_times.clear();

    // set some times at the beginning of the wf to ignore for llh calculation!
    double ignore_start = 0.0;
    double ignore_end = 5.5;

    // let's try looping over our llh grid
    for (size_t i = 0; i < llh_grid.size(); i++){
        double this_time = llh_grid.at(i);

        // grab our data
        double k = Interpolation<double>(this_time, data_times, pmt_data.data);
        
        // now call our interpolation function for pred yields
        T pred_yields_this_time = Interpolation<T>(this_time, pred_times_offset, *pred_yields);
        T pred_yields_squared_this_time = Interpolation<T>(this_time, pred_times_offset, *pred_yields_squared);

        // and grab our late pulse gaussian fit
        T relative_late_pulse_height = late_pulse_scale * normalization;
        T late_pulse_gauss = LatePulseGaussian<T>(this_time, late_pulse_mu, late_pulse_sigma, relative_late_pulse_height);

        // now put together to make mu and sigma2!
        T mu = pred_yields_this_time * normalization + pre_event_average + late_pulse_gauss;
        T sigma_squared = pred_yields_squared_this_time * (normalization * normalization) + (pre_event_average * pre_event_average) + (late_pulse_gauss * late_pulse_gauss);

        // now check the time before computing the llh
        if (this_time < ignore_start or this_time > ignore_end){
            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
            
            if constexpr (std::is_same<T, AD>::value){
                //std::cout << "late_pulse_scale = " << late_pulse_scale << std::endl;
                //std::cout << "mu = " << mu << std::endl;
                //std::cout << "at time = " << this_time << " nllh = " << MCLLH::LEff()(k, mu, sigma_squared) << std::endl;
                
                Grad_cppmin grad;
                MCLLH::LEff()(k, mu, sigma_squared).copyGradient(grad.data());
                debug_fit_data[key].at(i) = grad.at(11);
                //std::cout << "late pulse scale gradient = " << grad.at(14) << std::endl;
            }

                
        }
        if constexpr (std::is_same<T, double>::value) {
            debug_all_data[key].at(i) = k;
            debug_all_pred[key].at(i) = mu;
            debug_fit_pred[key].at(i) = late_pulse_gauss;
            debug_data_times.push_back(this_time);
        }
    }
    //std::cout << "total LP scale grad = " << total_LP_scale_grad << std::endl;

    //// now loop over the time bins in our data 
    //size_t pred_idx;
    //for (size_t i=0; i<n_bins; ++i) {
    //    bool in_bump_region = false;
    //    if (i >= bump_start_bin and i  <= bump_end_bin){
    //        in_bump_region = true;
    //    }
    //    double k = pmt_data.data.at(i);
    //    double time_in_data = data_times.at(i);

    //    // let's find closest time bin in prediction
    //    T closest_pred_time = time_in_data + light_time_offset;
    //    if constexpr (std::is_same<T, double>::value) {
    //        pred_idx = (size_t) ((closest_pred_time - minimum_pred_time) / 2.0);
    //    }
    //    if constexpr (std::is_same<T, AD>::value) {
    //        pred_idx = (size_t) ((closest_pred_time.value() - minimum_pred_time.value()) / 2.0);
    //    }
    //    
    //    T mu = pred_yields->at(pred_idx) * normalization + pre_event_average;
    //    T sigma_squared = pred_yields_squared->at(pred_idx) * (normalization * normalization) + (pre_event_average * pre_event_average);
    //    
    //    // we want to compute the llh everywhere except first 3 bins of event and the bump region
    //    if ((i < event_start_bin) or (i > (event_start_bin + 3))){
    //        if (in_bump_region == false){
    //            total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
    //            if constexpr (std::is_same<T, double>::value) {
    //                debug_fit_data[key].at(i) = k;
    //                debug_fit_pred[key].at(i) = mu;
    //            }
    //        }
    //    }
    //    //total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
    //    if constexpr (std::is_same<T, double>::value) {
    //        debug_all_data[key].at(i) = k;
    //        debug_all_pred[key].at(i) = mu;
    //        debug_pred_times.push_back(closest_pred_time);
    //        debug_data_times.push_back(time_in_data);
    //    }
    //    // let's also save the gradient of time offset at every time bin
    //    // over-writing debug_fit_pred to contain gradient of t offset in every bin fyi
    //    //if constexpr (std::is_same<T, AD>::value) {
    //    //    Grad grad;
    //    //    MCLLH::LEff()(k, mu, sigma_squared).copyGradient(grad.data());
    //    //    debug_fit_pred[key].at(i) = grad.at(7);
    //    //}

    //}
    return total_nllh;
}


#endif // CalculateNLLH_H
