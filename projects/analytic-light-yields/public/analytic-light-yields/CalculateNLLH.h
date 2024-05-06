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
#include <analytic-light-yields/GenerateExpectation.h>

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
        items[0] = ((2 * log(w * w) / w2) + (log(1 + ((w * w) / w2)) / w2)) * DwDTheta -
                   ((log(w * w * w) / (w2 * w2 * w2)) + (log(w) * (1 + (w * w) / w2 * w2) / (w2 * w2))) * Dw2DTheta;
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

class CalculateNLLH {
    std::shared_ptr<GenerateExpectation> gen_expectation = nullptr;
    I3MapPMTKeyVectorDouble data;
    double n_data_events;
    bool grabbed_data = false;

public:
    CalculateNLLH();
    void GrabData(I3FramePtr data_frame);
    double GetNLLH(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
                   std::vector<double> nuisance_params, I3FramePtr data_frame, size_t time_bin_offset);
    I3Vector<double> GetNLLHDerivative(AnalyticLightYieldGenerator analytic_light_yield_setup, I3FramePtr geo_frame,
                                       std::vector<double> nuisance_params, I3FramePtr data, size_t time_bin_offset);
};
#endif
