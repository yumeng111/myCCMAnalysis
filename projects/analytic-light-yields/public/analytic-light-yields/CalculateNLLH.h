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
#include <icetray/I3Frame.h>
#include "dataclasses/I3Map.h"
#include "dataclasses/I3Vector.h"
#include "dataclasses/I3Double.h"
#include "analytic-light-yields/GenerateExpectation.h"
#include <analytic-light-yields/YieldsPerPMT.h>
#include "analytic-light-yields/autodiff.h"
#include "simclasses/CCMMCPE.h"
#include "simclasses/PhotonSummary.h"

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

template<typename T>
class CalculateNLLH {
public:

    std::vector<double> data_vector;
    std::vector<double> pred_vector;
    std::vector<double> times_vector;

    double event_start_threshold = 10.0;

private:
    std::map<CCMPMTKey, SinglePMTInfo> data;
    I3VectorCCMPMTKey keys_to_fit;
    size_t max_bins;
    size_t n_data_events;

public:
    CalculateNLLH();
    CalculateNLLH(I3FramePtr data_frame, size_t max_bins=100, I3VectorCCMPMTKey keys_to_fit=I3VectorCCMPMTKey());
    void SetKeys(I3VectorCCMPMTKey keys);
    void SetData(I3FramePtr data_frame);
    void SetGenExpectation(boost::shared_ptr<GenerateExpectation<T>> gen_expectation);

    double GetNDataEvents() const { return n_data_events; };
    I3MapPMTKeyVectorDoublePtr GetData() const;
    boost::shared_ptr<GenerateExpectation<T>> GetGenExpectation() const { return gen_expectation; };

    std::vector<double> GetDataVector() {return data_vector;};
    std::vector<double> GetPredVector() {return pred_vector;};
    std::vector<double> GetTimesVector() {return times_vector;};

    std::vector<double> GetLightProfileDebug() {return gen_expectation->GetLightProfileDebug();};
    std::vector<double> GetLightProfileTimesDebug() {return gen_expectation->GetLightProfileTimesDebug();};
    std::vector<double> GetLightProfileTOffsetGradDebug() {return gen_expectation->GetLightProfileTOffsetGradDebug();};

    T ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB, T normalization, T light_time_offset,
            T const_offset, T late_pulse_mu, T late_pulse_sigma, T late_pulse_scale, T pmt_efficiency,
            std::vector<T> uv_absorption, T rayl, T photons_per_mev, T z_offset, size_t n_sodium_events,
            AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, bool UseG4Yields);

    boost::shared_ptr<GenerateExpectation<T>> gen_expectation = nullptr;
};

template<typename T>
T CalculateNLLH<T>::ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
                                T normalization, T light_time_offset, T const_offset, T late_pulse_mu,
                                T late_pulse_sigma, T late_pulse_scale, T pmt_efficiency, std::vector<T> uv_absorption,
                                T rayl, T photons_per_mev, T z_offset, size_t n_sodium_events,
                                AnalyticLightYieldGenerator::LArLightProfileType light_profile_type, bool fix_z_rayl) {

    //std::cout << "In CalculateNLLH<T>::ComputeNLLH : Rs = " << Rs << ", tau_s = " << tau_s << ", tau_TPB = " << tau_TPB <<
    //    ", norm = " << normalization << ", pmt eff = " << pmt_efficiency << ", uv abs = " << uv_absorption << ", photons/mev = "
    //    << photons_per_mev << ", rayl = " << rayl << std::endl;
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
    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> pred = gen_expectation->GetExpectation(key, start_time, max_time,
                                                            peak_time, Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, light_time_offset, late_pulse_mu, late_pulse_sigma, late_pulse_scale,
                                                            uv_absorption, rayl, z_offset, n_sodium_events, light_profile_type, fix_z_rayl);
    // unpack into yields, yields^2, and times
    boost::shared_ptr<std::vector<T>> pred_yields;
    boost::shared_ptr<std::vector<T>> pred_yields_squared;
    pred_yields = std::get<0>(pred);
    pred_yields_squared = std::get<1>(pred);

    T total_nllh = 0.0;

    //size_t n_bins = pmt_data.data.size();

    //// let's make a llh grid -- this combines both the times we have data values at along with the times we have prediction values at
    //double starting_llh_grid_time = data_times.front();
    //double ending_llh_grid_time = data_times.back();
    //
    //// to make things easier, let's make sure we have a vector of doubles containing pred times
    //std::vector<double> pred_times_double;
    //for (size_t j = 0; j < pred_times->size(); j++){
    //    if constexpr (std::is_same<T, double>::value) {
    //        pred_times_double.push_back(pred_times->at(j) + light_time_offset);
    //    } else {
    //        pred_times_double.push_back(pred_times->at(j).value() + light_time_offset.value());
    //    }
    //}


    //std::vector<double> llh_grid;
    //for (size_t i = 0; i < data_times.size(); i++){
    //    llh_grid.push_back(data_times.at(i));
    //}
    //for (size_t i = 0; i < pred_times_double.size(); i++){
    //    if (pred_times_double.at(i) >= starting_llh_grid_time and pred_times_double.at(i) <= ending_llh_grid_time){
    //        llh_grid.push_back(pred_times_double.at(i));
    //    }
    //}

    //std::vector<T> pred_times_offset;
    //for (size_t j = 0; j < pred_times->size(); j++){
    //    pred_times_offset.push_back(pred_times->at(j) + light_time_offset);
    //}
    //
    //// now sort our llh grid!
    //std::sort(llh_grid.begin(), llh_grid.end());

    //// let's empty vectors to save data, pred, and times
    //data_vector.clear();
    //pred_vector.clear();
    //times_vector.clear();

    //// set some times at the beginning of the wf to ignore for llh calculation!
    //double ignore_start = 0.0;
    //double ignore_end = 6.0;
    //double n_event_scaling = n_data_events / (double)n_sodium_events;

    //// let's try looping over our llh grid
    //for (size_t i = 0; i < llh_grid.size(); i++){
    //    double this_time = llh_grid.at(i);
    //
    //    // grab our data
    //    double k = InterpolationDouble(this_time, data_times, pmt_data.data);

    //    // now call our interpolation function for pred yields
    //    T pred_yields_this_time = Interpolation(this_time, pred_times_offset, *pred_yields);
    //    T pred_yields_squared_this_time = Interpolation(this_time, pred_times_offset, *pred_yields_squared);

    //    // now put together to make mu and sigma2!
    //    T adjusted_norm = normalization * pmt_efficiency;
    //    T mu = (pred_yields_this_time * adjusted_norm * n_event_scaling) + pre_event_average;
    //    T sigma_squared = (pred_yields_squared_this_time * (adjusted_norm * adjusted_norm) * (n_event_scaling * n_event_scaling) )+ (pre_event_average * pre_event_average);

    //    // now check the time before computing the llh
    //    if (this_time < ignore_start or this_time > ignore_end){
    //        T partial_nllh = MCLLH::LEff()(k, mu, sigma_squared);
    //        total_nllh += partial_nllh;
    //    }

    //    // and save some stuff for plotting
    //    data_vector.push_back(k);
    //    times_vector.push_back(this_time);
    //    if constexpr (std::is_same<T, double>::value) {
    //        pred_vector.push_back(mu);
    //    } else {
    //        pred_vector.push_back(mu.value());
    //    }
    //}
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

template<typename T> CalculateNLLH<T>::CalculateNLLH() : max_bins(100),  keys_to_fit(I3VectorCCMPMTKey()) {}

template<typename T> CalculateNLLH<T>::CalculateNLLH(I3FramePtr data_frame, size_t max_bins, I3VectorCCMPMTKey keys_to_fit) :
    max_bins(max_bins), keys_to_fit(keys_to_fit) {

    SetData(data_frame);
}

template<typename T> void CalculateNLLH<T>::SetKeys(I3VectorCCMPMTKey keys){
    keys_to_fit = keys;
    gen_expectation = boost::make_shared<GenerateExpectation<T>>(keys_to_fit);
}

template<typename T> void CalculateNLLH<T>::SetData(I3FramePtr data_frame) {
    bool restrict_keys = keys_to_fit.size() > 0;
    // grab our data out of the frame
    I3MapPMTKeyVectorDoubleConstPtr data_map = data_frame->Get<I3MapPMTKeyVectorDoubleConstPtr>("AccumulatedEventsMap");
    for(I3MapPMTKeyVectorDouble::const_iterator i = data_map->begin(); i != data_map->end(); i++) {
        if(restrict_keys and std::find(keys_to_fit.begin(), keys_to_fit.end(), i->first) == keys_to_fit.end()) {
            continue;
        }
        // based on the way I accumulated data, the event starts around the 50th bin
        // so to be safe, let's use the first 40 bins for out pre-event average
        double pre_event_values = 0.0;
        double pre_event_deriv = 0.0;
        double pre_event_bins = 0.0;
        double time_counter = 0.0;
        double prev_data = 0.0;
        for(double const & data_points: i->second) {
            if (time_counter < 80.0){
                pre_event_values += data_points;
                pre_event_bins += 1.0;
                pre_event_deriv += abs(data_points - prev_data);
            }
            prev_data = data_points;
            time_counter += 2.0;
        }

        SinglePMTInfo pmt_data;
        pmt_data.key = i->first;
        size_t peak_idx = std::distance(i->second.begin(), std::max_element(i->second.begin(), i->second.end()));
        pmt_data.peak_value = i->second.at(peak_idx);
        size_t start_idx = std::max(size_t(3), peak_idx) - 3;
        size_t min_idx = std::max(size_t(15), peak_idx) - 15;
        size_t max_idx = std::min(min_idx + max_bins, i->second.size());
        pmt_data.data = std::vector(i->second.begin() + min_idx, i->second.begin() + max_idx);
        pmt_data.start_time = (start_idx - min_idx) * 2.0;
        pmt_data.max_time = (std::max(int(min_idx), int(max_idx) - 1) - min_idx) * 2.0 + 30.0;
        pmt_data.peak_time = (peak_idx - min_idx) * 2.0;

        // let's also get the event start bin using derivs
        bool found_start = false;
        double deriv_threshold = (pre_event_deriv / pre_event_bins) * 5.0;
        pmt_data.data_times.push_back(0.0);
        for (size_t data_it = 1; data_it < pmt_data.data.size(); data_it++){
            double deriv = pmt_data.data.at(data_it) - pmt_data.data.at(data_it - 1);
            if (deriv > deriv_threshold and found_start == false){
                pmt_data.event_start_bin = data_it - 1;
                found_start = true;
            }
            pmt_data.data_times.push_back(pmt_data.data_times[data_it - 1] + 2.0);
        }

        // let's subtract off our event start time from pmt_data.data_times
        for (size_t i = 0; i < pmt_data.data_times.size(); i++){
            pmt_data.data_times.at(i) -= (pmt_data.event_start_bin * 2.0);
        }
        pmt_data.pre_event_average = pre_event_values / pre_event_bins;
        data[pmt_data.key] = pmt_data;
    }
    n_data_events = data_frame->Get<I3Double>("TotalEventsPastCuts").value;
}

template<typename T> I3MapPMTKeyVectorDoublePtr CalculateNLLH<T>::GetData() const {
    I3MapPMTKeyVectorDoublePtr data_map = boost::make_shared<I3MapPMTKeyVectorDouble> ();
    for (std::map<CCMPMTKey, SinglePMTInfo>::const_iterator i = data.begin(); i != data.end(); i++) {
        data_map->insert(std::make_pair(i->first, i->second.data));
    }
    return data_map;
}

template<typename T> void CalculateNLLH<T>::SetGenExpectation(boost::shared_ptr<GenerateExpectation<T>> gen_expectation) {
    this->gen_expectation = gen_expectation;
}


#endif // CalculateNLLH_H
