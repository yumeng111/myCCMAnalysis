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
#include "analytic-light-yields/CalculateNLLH.h"
#include "analytic-light-yields/AllPMTcppMinimizer.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"

struct NewLikelihoodFunctor {

    static constexpr int NewDerivativeDimension = 7; // max number of dimensions :
                                                     // Rs, tau_s, tau_TPB
                                                     // norm, norm, norm, norm -- 1 / data set!
                                                     
                                                     // things that are fixed for all PMTs:
                                                     // Rt, tau_t, tau_rec, const_offset 
                                                     
                                                     // things that are fixed for each PMT : 
                                                     // time offset (1 / PMT / data set), LPmu, LPsigma, and LPscale
    
    std::vector<std::shared_ptr<CalculateNLLH>> llh_constructor; // list of constructors
    I3VectorCCMPMTKey all_keys; // list of PMTs in fit
    I3MapPMTKeyDouble LPmu; // map between CCMPMTKey and LPmu
    I3MapPMTKeyDouble LPsigma; // map between CCMPMTKey and LPsigma
    I3MapPMTKeyDouble LPscale; // map between CCMPMTKey and LPscale
    std::vector<I3MapPMTKeyDouble> time_offsets; // map between CCMPMTKey and time offset (1 map/data set) 
    double uv_absorption = 55.0;
    double Rt = 0.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    std::vector<double> z_offset; // list of z offsets
    size_t n_sodium_events;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        T total_llh = 0;
        for (size_t data_it = 0; data_it < llh_constructor.size(); data_it ++){
            // loop over each data set we are fitting to
            for (size_t pmt_it = 0; pmt_it < all_keys.size(); pmt_it ++){
                // loop over each PMT
                CCMPMTKey key = all_keys.at(pmt_it);
                total_llh += llh_constructor.at(data_it)->ComputeNLLH<T>(key, x[0], Rt, x[1], tau_t, tau_rec, x[2], x[3 + data_it], // normalization for data set!!!
                                 time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), uv_absorption,
                                 z_offset.at(data_it), n_sodium_events, light_profile_type);
            }
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

typedef NewLikelihoodFunctor LikelihoodType;

AllPMTcppMinimizer::AllPMTcppMinimizer() {}

std::vector<double> AllPMTcppMinimizer::MultiplePMTMinimization(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble PMT_efficiencies,
                                                                I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                                                I3MapPMTKeyDouble time_offset1, I3MapPMTKeyDouble time_offset2, I3MapPMTKeyDouble time_offset3, I3MapPMTKeyDouble time_offset4,
                                                                std::vector<std::string> data_file_names,
                                                                std::vector<double> z_offsets, size_t n_sodium_events) {        

    std::vector<I3MapPMTKeyDouble> time_offsets = {time_offset1, time_offset2, time_offset3, time_offset4}; 
    std::vector<std::string> paramter_names = {"Rs", "tau_s", "tau_TPB", "norm1", "norm2", "norm3", "norm4"};

    double seeds[4] = {0.34, 6.7, 6.7, 1e5};
    double grad_scales[4] = {1e-3, 1e-3, 1e-3, 1e-6};
    double mins[4] = {0.2, 2.0, 1.0, 10.0};
    double maxs[4] = {0.5, 16.0, 9.0, 1e10};
    
    // let's initialize our constructor
    // grab our geometry frame
    std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
    dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
    I3FramePtr geo_frame = geometry_file.pop_frame();
    
    // now grab data frame(s) to make our constructors
    std::vector<std::shared_ptr<CalculateNLLH>> all_constructors;
    for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
        std::string this_data_fname = data_file_names.at(data_it);
        dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
        I3FramePtr this_data_frame = this_data_file.pop_frame();

        // and now make constructor
        std::shared_ptr<CalculateNLLH> this_llh_constructor = std::make_shared<CalculateNLLH>();
        this_llh_constructor->SetKeys(keys_to_fit);
        this_llh_constructor->SetPMTEff(PMT_efficiencies);
        this_llh_constructor->SetData(this_data_frame);
        this_llh_constructor->SetGeo(geo_frame);
    
        // now save!
        all_constructors.push_back(this_llh_constructor);
    }

    // now set up our likelihood object
    LikelihoodType likelihood;
    likelihood.n_sodium_events = n_sodium_events;
    likelihood.z_offset = z_offsets; 
    likelihood.llh_constructor = all_constructors;
    likelihood.all_keys = keys_to_fit;
    likelihood.LPmu = LPmu;
    likelihood.LPsigma = LPsigma;
    likelihood.LPscale = LPscale;
    likelihood.time_offsets = time_offsets;

    // now set up our minimizer
    phys_tools::lbfgsb::LBFGSB_Driver minimizer;
        
    minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]); // Rs
    minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]); // tau_s
    minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]); // tau_TPB
    // this is where we add our normalization parameter...need to add one for each data set
    for (size_t n = 0; n < 4; n++){
        minimizer.addParameter(seeds[3], grad_scales[3], mins[3], maxs[3]); // norm
    } 
    minimizer.setHistorySize(20);

    // fix parameter idx of guys we are not minimizing
    size_t data_sets_to_minimize = data_file_names.size();

    // data_sets_to_minimize should either be 1 or 4
    // but let's try to make this general
    size_t data_sets_included = 0;
    for (size_t n = 0; n < 4; n++){
        data_sets_included += 1;
        if (data_sets_included > data_sets_to_minimize){
            minimizer.fixParameter(2 + data_sets_included); // norms
        }
    } 
    
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
        for(size_t i=0; i<NewLikelihoodFunctor::NewDerivativeDimension; ++i) {
            data_to_return.push_back(params.at(i));
            std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
        }
    } else {
        std::cout << "oops! fit did not converge :( error message = " << minimizer.errorMessage() << std::endl;
    } 
    
    return data_to_return;

}



