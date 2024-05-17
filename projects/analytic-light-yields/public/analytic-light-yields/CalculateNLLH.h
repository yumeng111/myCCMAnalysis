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
            std::cout << "w/w2 < 0" << std::endl;
            std::cout << "w = " << w_sum << ", w2 = " << w2_sum << std::endl;
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
                std::cout << "w = 0" << std::endl;
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
    double start_time;
    double max_time;
};

class CalculateNLLH {
public:
    static constexpr size_t n_params = 9;
    typedef double Underlying;

    typedef phys_tools::autodiff::FD<n_params, Underlying> AD;
    typedef std::array<Underlying, n_params> Grad;

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
    CalculateNLLH(I3FramePtr data_frame, I3FramePtr geo_frame, size_t max_bins=125, size_t n_sodium_events=2000, double portion_light_reflected_by_tpb=1.0, double desired_chunk_width=20, double desired_chunk_height=20, std::vector<CCMPMTKey> keys_to_fit=std::vector<CCMPMTKey>());
    void SetKeys(std::vector<CCMPMTKey> keys);
    void SetGeo(I3FramePtr geo_frame);
    void SetData(I3FramePtr data_frame);
    void SetGenExpectation(boost::shared_ptr<GenerateExpectation> gen_expectation);

    double GetNDataEvents() const { return n_data_events; };
    I3MapPMTKeyVectorDoublePtr GetData() const;
    boost::shared_ptr<GenerateExpectation> GetGenExpectation() const { return gen_expectation; };

    template<typename T>
    T ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, T const_offset, double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type);

    AD GetNLLH(CCMPMTKey key, AnalyticLightYieldGenerator const & params);
    double GetNLLHValue(CCMPMTKey key, AnalyticLightYieldGenerator const & params);
    std::vector<double> GetNLLHDerivative(CCMPMTKey key, AnalyticLightYieldGenerator const & params);
};

template<typename T>
T CalculateNLLH::ComputeNLLH(CCMPMTKey key, T Rs, T Rt, T tau_s, T tau_t, T tau_rec, T tau_TPB,
            T normalization, T light_time_offset, T const_offset, double uv_absorption, double z_offset, size_t n_sodium_events, AnalyticLightYieldGenerator::LArLightProfileType light_profile_type) {

    // let's grab our data
    SinglePMTInfo const & pmt_data = data[key];
    double start_time = pmt_data.start_time;
    double max_time = pmt_data.max_time;
    double peak_time = pmt_data.peak_time;
    // let's grab our expectation
    std::tuple<boost::shared_ptr<std::vector<T>>, boost::shared_ptr<std::vector<T>>> pred = gen_expectation->GetExpectation(key, start_time, max_time, peak_time, Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB,
            normalization, light_time_offset, uv_absorption, z_offset, n_sodium_events, light_profile_type);
    
    // unpack into yields and yields^2
    boost::shared_ptr<std::vector<T>> pred_yields;
    boost::shared_ptr<std::vector<T>> pred_yields_squared;
    pred_yields = std::get<0>(pred);
    pred_yields_squared = std::get<1>(pred);

    T total_nllh(0.0);

    size_t n_bins = pmt_data.data.size();
    
    // now loop over the time bins in our pred
    for (size_t i=0; i<n_bins; ++i) {
        double k = pmt_data.data.at(i);
        T mu = pred_yields->at(i) * normalization + const_offset;
        T sigma_squared = pred_yields_squared->at(i) * (normalization * normalization) + (const_offset * const_offset);
        total_nllh += MCLLH::LEff()(k, mu, sigma_squared);
    }
    return total_nllh;
}


#endif // CalculateNLLH_H
