#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions.hpp>

#include <string>
#include <vector>
#include <iostream>

#include "icetray/ctpl.h"
#include "icetray/open.h"
#include "icetray/I3Frame.h"
#include "icetray/I3Units.h"
#include "icetray/I3TrayInfo.h"
#include "icetray/I3Module.h"
#include "icetray/I3Logging.h"
#include "icetray/CCMPMTKey.h"
#include "dataio/I3File.h"
#include "dataclasses/geometry/CCMGeometry.h"
#include "dataclasses/I3Map.h"
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

#include "analytic-light-yields/lbfgsb.h"
#include "analytic-light-yields/G4YieldsPerPMT.h"
#include "analytic-light-yields/CalculateNLLH.h"
#include "analytic-light-yields/StandAloneUVAbsMinimizer.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"

struct UVAbsLikelihoodFunctor {

    static constexpr int NewDerivativeDimension = 2; // max number of dimensions :
                                                     // uv absoprtion and scaling 

    std::shared_ptr<G4YieldsPerPMT> yields_constructor;
    std::vector<CCMPMTKey> all_keys; // list of PMTs in fit
    std::map<CCMPMTKey, double> data; // data!!!
    

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        if constexpr (std::is_same<T, double>::value) {
            std::cout << "uv abs = " << x[0] << ", scaling = " << x[1] << std::endl; 
        } else {
            std::cout << "uv abs = " << x[0].value() << ", scaling = " << x[1].value() << std::endl; 
        }
        
        // first grab the summed yields for given uv abs and scaling
        std::tuple<std::map<CCMPMTKey, T>, std::map<CCMPMTKey, T>> yields = yields_constructor->GetSummedYieldsMap<T>(all_keys, x[0], x[1]); 
        
        std::map<CCMPMTKey, T> summed_yields;
        std::map<CCMPMTKey, T> summed_yields_squared;
        summed_yields = std::get<0>(yields);
        summed_yields_squared = std::get<1>(yields);

        // now loop over each key, calculate LLH 
        T total_llh = 0;
        for (size_t pmt_it = 0; pmt_it < all_keys.size(); pmt_it ++){
            // loop over each PMT
            CCMPMTKey key = all_keys.at(pmt_it);
            total_llh += MCLLH::LEff()(data.at(key), summed_yields.at(key), summed_yields_squared.at(key));
        }
        return total_llh; 
    }
};

template<typename FuncType>
class BFGS_Function : public phys_tools::lbfgsb::BFGS_FunctionBase {
private:
    FuncType func;
public:
    BFGS_Function(FuncType f) : func(f) {}
    
    // This evaluates the NLLH value
    virtual double evalF(std::vector<double> x) const {
        return(-func.template evaluateLikelihood<double>(x)); // note negation!
    }

    // This evaluates the NLLH value and the NLLH gradient
    virtual std::pair<double,std::vector<double>> evalFG(std::vector<double> x) const {
        const size_t size=x.size();
        using GradType = phys_tools::autodiff::FD<FuncType::NewDerivativeDimension>;
        std::vector<GradType> params(size);
        for(size_t i=0; i<size; i++)                                                                                                                                                         
            params[i] = GradType(x[i],i);
        GradType result=func.template evaluateLikelihood<GradType>(params);
        std::vector<double> grad(size);
        for(unsigned int i=0; i<size; i++)
            grad[i] = -result.derivative(i); // note negation!
        return(std::make_pair(-result.value(), grad)); // note negation!
    }
};

typedef UVAbsLikelihoodFunctor LikelihoodType;

StandAloneUVAbsMinimizer::StandAloneUVAbsMinimizer() {}

double StandAloneUVAbsMinimizer::GrabScalingSeed(CCMPMTKey key, double data,  double uv_abs, std::shared_ptr<G4YieldsPerPMT> yields_constructor){
   
    // grab yields
    double baseline_scaling = 1.0;
    std::tuple<std::map<CCMPMTKey, double>, std::map<CCMPMTKey, double>> yields = yields_constructor->GetSummedYieldsMap<double>({key}, uv_abs, baseline_scaling);
    std::map<CCMPMTKey, double> summed_yields = static_cast<I3MapPMTKeyDouble>(std::get<0>(yields));  
       
    // now get norm scaling to make summed_yields == data
    baseline_scaling *= (data / summed_yields.at(key)); 

    return baseline_scaling; 
}


std::vector<double> StandAloneUVAbsMinimizer::MinimizeUVAbsorption(std::vector<CCMPMTKey> keys_to_fit, I3MapPMTKeyDouble data){

    std::vector<std::string> paramter_names = {"UV Abs", "scaling"};
    
    // now set up our likelihood object
    LikelihoodType likelihood;
    
    // make our G4YieldsPerPMT constructors 
    std::shared_ptr<G4YieldsPerPMT> yields_constructor = std::make_shared<G4YieldsPerPMT>();
    
    // now save to likelihood object
    likelihood.yields_constructor = yields_constructor;
    likelihood.all_keys = keys_to_fit;
    likelihood.data = data;

    // now set up our minimizer
    phys_tools::lbfgsb::LBFGSB_Driver minimizer;
    
    double uv_abs_seed = 55.0;
    size_t reference_pmt_idx = (size_t) (keys_to_fit.size() / 2);
    CCMPMTKey reference_pmt = keys_to_fit.at(reference_pmt_idx);
    double scaling_seed = GrabScalingSeed(reference_pmt, data.at(reference_pmt), uv_abs_seed, yields_constructor); 
    std::cout << "for uv abs = " << uv_abs_seed << " scaling seed = " << scaling_seed << std::endl;

    minimizer.addParameter(uv_abs_seed, 1e-3, 5.0, 100.0); // uv abs
    minimizer.addParameter(scaling_seed, 1e-3, scaling_seed * 1e-5, scaling_seed * 1e5); // scaling
    minimizer.setHistorySize(20);
    
    bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
 
    std::vector<double> data_to_return;
    
    if(succeeded) {
        std::cout << "joint fit converged!" << std::endl;
        std::cout << "minimization finished with " << minimizer.errorMessage() << " error message" << std::endl; 
        
        // save things for returning!
        double value = minimizer.minimumValue();                                                                                                                                                 
        std::vector<double> params = minimizer.minimumPosition();
        
        std::cout << "Function value at minimum: " << value << " after " << minimizer.numberOfEvaluations() << " function evaluations" << std::endl;
        std::cout << "Parameters: " << std::endl;
        data_to_return.push_back(value);
        data_to_return.push_back((double) minimizer.numberOfEvaluations());
        for(size_t i=0; i<UVAbsLikelihoodFunctor::NewDerivativeDimension; ++i) {
            data_to_return.push_back(params.at(i));
            std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
        }
    
        // one last thing -- take this best fit point and grab our pred 
        std::tuple<std::map<CCMPMTKey, double>, std::map<CCMPMTKey, double>> yields = yields_constructor->GetSummedYieldsMap<double>(keys_to_fit, data_to_return.at(2), data_to_return.at(3));
        best_fit_summed_yields = static_cast<I3MapPMTKeyDouble>(std::get<0>(yields));  
    
    } else {
        std::cout << "minimizer did not converge :( error message = " << minimizer.errorMessage() << std::endl;
    }
    
    return data_to_return;

}


void StandAloneUVAbsMinimizer::MakeSpatialDistros(std::vector<CCMPMTKey> keys_to_fit, std::vector<double> uv_abs, double scaling){

    // make our G4YieldsPerPMT constructors 
    std::shared_ptr<G4YieldsPerPMT> yields_constructor = std::make_shared<G4YieldsPerPMT>();
    
    for(size_t i = 0; i < uv_abs.size(); i++){
        std::tuple<std::map<CCMPMTKey, double>, std::map<CCMPMTKey, double>> yields = yields_constructor->GetSummedYieldsMap<double>(keys_to_fit, uv_abs.at(i), scaling);
        best_fit_summed_yields_multiple_uv_abs.push_back(static_cast<I3MapPMTKeyDouble>(std::get<0>(yields)));  
    }

}





