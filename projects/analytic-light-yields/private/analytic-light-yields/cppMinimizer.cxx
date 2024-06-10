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
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"

#include "analytic-light-yields/lbfgsb.h"
#include "analytic-light-yields/CalculateNLLH.h"
#include "analytic-light-yields/cppMinimizer.h"
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"

struct LikelihoodFunctor {

    static constexpr int DerivativeDimension = 18; // max number of dimensions :
                                                   // Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB
                                                   // norm, norm, norm, norm -- 1 / data set!
                                                   // time offset, const offset
                                                   // LP mu, LP sigma,
                                                   // LP scale, LP scale, LP scale, LP scale -- 1 / data set!
    std::vector<std::shared_ptr<CalculateNLLH>> llh_constructor;
    CCMPMTKey key_to_fit;
    double uv_absorption = 55.0;
    std::vector<double> z_offset;
    size_t n_sodium_events;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
    std::vector<double> time_offsets;

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        T total_llh = 0;
        for (size_t data_it = 0; data_it < llh_constructor.size(); data_it ++){
            T time_offset;
            if (time_offsets.size() > 0){
                time_offset = time_offsets.at(data_it);
            } else {
                time_offset = x[10];
            }
            T partial_nllh = llh_constructor.at(data_it)->ComputeNLLH<T>(key_to_fit, x[0], x[1], x[2], x[3], x[4], x[5], x[6 + data_it], // normalization for data set!!!
                             time_offset, x[11], x[12], x[13], x[14 + data_it], uv_absorption, z_offset.at(data_it), n_sodium_events, light_profile_type);
            total_llh += partial_nllh;
            //std::cout << "partial nllh = " << partial_nllh << std::endl; 
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
        using GradType = phys_tools::autodiff::FD<FuncType::DerivativeDimension>;
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


typedef LikelihoodFunctor LikelihoodType;

cppMinimizer::cppMinimizer() {}

std::vector<double> cppMinimizer::OnePMTOneDataSetMinimization(CCMPMTKey this_key, std::vector<std::string> data_file_names, std::vector<double> z_offsets, double norm_seed, size_t n_sodium_events) {

    // let's try scanning over t offset
    double t_offset_start = 25.0;
    double t_offset_end = 35.0;
    double t_offset_range = t_offset_end - t_offset_start;
    size_t n_t_offsets = 20;

    std::vector<double> func_val;
    std::vector<std::vector<double>> params;
    std::vector<int> nevals;
    std::vector<std::string> paramter_names = {"Rs", "Rt", "tau_s", "tau_t", "tau_rec", "tau_TPB", "norm1", "norm2", "norm3", "norm4",  "time offset", "const offset",
                                               "late pulse mu", "late pulse sigma", "late pulse scale1", "late pulse scale2", "late pulse scale3", "late pulse scale4"};

    // free paramters are : Rs, tau_s, tau_TPB, normalization, late pulse mu, late pulse sigma, and late pulse scale
    // fixed paramters are : Rt, tau_t, tau_rec, time offset, and const_offset
    double seeds[7] = {0.34, 8.2, 3.0, norm_seed, 58.0, 8.0, 1.0};
    double grad_scales[7] = {1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6, 1e-6};
    double mins[7] = {0.2, 6.0, 1.0, norm_seed / 5, 50.0, 2.0, 0.0};
    double maxs[7] = {0.5, 16.0, 9.0, norm_seed * 1e5, 68.0, 15.0, 30.0};

    std::cout << "norm seed = " << norm_seed << ", lower bound = " << mins[3] << " and upper bound = " << maxs[3] << std::endl;
    
    // let's initialize our constructor
    // grab our geomtry frame
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
        this_llh_constructor->SetKeys({this_key});
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
    likelihood.key_to_fit = this_key;

    for (size_t toff_it = 0; toff_it < (n_t_offsets + 1); toff_it++){
        double this_t_offset = t_offset_start + ((double)toff_it / (double)n_t_offsets) * t_offset_range;
        
        phys_tools::lbfgsb::LBFGSB_Driver minimizer;

        minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]); // Rs
        minimizer.addParameter(0.0);                                        // Rt
        minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]); // tau_s
        minimizer.addParameter(743.0);                                      // tau_t
        minimizer.addParameter(0.0);                                        // tau_rec
        minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]); // tau_TPB
        // this is where we add our normalization parameter...need to add one for each data set
        for (size_t n = 0; n < 4; n++){
            minimizer.addParameter(seeds[3], grad_scales[3], mins[3], maxs[3]); // norm
        } 
        minimizer.addParameter(this_t_offset);                              // time offset
        minimizer.addParameter(0.0);                                        // const offset
        minimizer.addParameter(seeds[4], grad_scales[4], mins[4], maxs[4]); // late pulse mu
        minimizer.addParameter(seeds[5], grad_scales[5], mins[5], maxs[5]); // late pulse sigma
        // this is where we add our LP scale parameter...need to add one for each data set
        for (size_t n = 0; n < 4; n++){
            minimizer.addParameter(seeds[6], grad_scales[6], mins[6], maxs[6]); // late pulse scale 
        } 
        minimizer.setHistorySize(20);

        // fix parameter idx of guys we are not minimizing
        minimizer.fixParameter(1); // Rt
        minimizer.fixParameter(3); // tau_t
        minimizer.fixParameter(4); // tau_rec
        minimizer.fixParameter(7); // norm #2 
        minimizer.fixParameter(8); // norm #3
        minimizer.fixParameter(9); // norm #4 
        minimizer.fixParameter(10); // time offset
        minimizer.fixParameter(11); // const offset
        minimizer.fixParameter(15); // LP scale #2 
        minimizer.fixParameter(16); // LP scale #3 
        minimizer.fixParameter(17); // LP scale #4 

        bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));

        if(succeeded) {
            func_val.push_back(minimizer.minimumValue());
            params.push_back(minimizer.minimumPosition());
            nevals.push_back(minimizer.numberOfEvaluations());
            std::cout << "function minimum at " << minimizer.minimumValue() << " and minimization message = " << minimizer.errorMessage() << std::endl; 
        } 

    }

    // now let's find the smallest func_val!
    size_t smallest_idx = std::distance(func_val.begin(), std::min_element(func_val.begin(), func_val.end()));
    std::cout << "Best evaluation!!! Function value at minimum: " << func_val.at(smallest_idx) << " after " << nevals.at(smallest_idx) << " function evaluations" << std::endl;
    std::cout << "Parameters: " << std::endl;
    std::vector<double> data_to_return;
    data_to_return.push_back(func_val.at(smallest_idx));
    for(size_t i=0; i<LikelihoodFunctor::DerivativeDimension; ++i) {
        data_to_return.push_back(params.at(smallest_idx).at(i));
        if (i == 1 or i == 3 or i == 4 or i == 7 or i == 8 or i == 9 or i == 11 or i == 15 or i == 16 or i == 17){
            continue;
        }
        std::cout << paramter_names.at(i) << " = " << params.at(smallest_idx).at(i) << std::endl;
    }
    return data_to_return;
}

std::vector<double> cppMinimizer::OnePMTMultipleDataSetMinimization(CCMPMTKey this_key, std::vector<std::string> data_file_names, std::vector<double> z_offsets,
                            std::vector<double> time_offsets, std::vector<double> LPscale_seeds, std::vector<double> norm_seeds, size_t n_sodium_events) {

    std::vector<std::string> paramter_names = {"Rs", "Rt", "tau_s", "tau_t", "tau_rec", "tau_TPB", "norm1", "norm2", "norm3", "norm4",
                                               "time offset", "const offset", "late pulse mu", "late pulse sigma", "late pulse scale1",
                                               "late pulse scale2", "late pulse scale3", "late pulse scale4"};

    size_t n_data_sets = data_file_names.size();

    // grab our geomtry frame
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
        this_llh_constructor->SetKeys({this_key});
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
    likelihood.key_to_fit = this_key;
    likelihood.time_offsets = time_offsets;

    phys_tools::lbfgsb::LBFGSB_Driver minimizer;

    // free paramters are : Rs, tau_s, tau_TPB, normalization, late pulse mu, late pulse sigma, and late pulse scale
    // fixed paramters are : Rt, tau_t, tau_rec, time offset, and const_offset
    double seeds[7] = {0.34, 8.2, 3.0, 1e5, 60.0, 4.0, 1.0};
    double grad_scales[7] = {1e-3, 1e-3, 1e-3, 1e-6, 1e-3, 1e-3, 1e-3};
    double mins[7] = {0.2, 6.0, 1.0, 1.0, 50.0, 2.0, 1e-1};
    double maxs[7] = {0.5, 16.0, 9.0, 1e10, 68.0, 15.0, 30.0};
    
    minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]); // Rs
    minimizer.addParameter(0.0);                                        // Rt
    minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]); // tau_s
    minimizer.addParameter(743.0);                                      // tau_t
    minimizer.addParameter(0.0);                                        // tau_rec
    minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]); // tau_TPB
    // this is where we add our normalization parameter...need to add one for each data set
    for (size_t n = 0; n < n_data_sets; n++){
        minimizer.addParameter(norm_seeds.at(n), grad_scales[3], norm_seeds.at(n) / 1e2, norm_seeds.at(n) * 1e5); // norm
    } 
    minimizer.addParameter(0.0);                                        // time offset
    minimizer.addParameter(0.0);                                        // const offset
    minimizer.addParameter(seeds[4], grad_scales[4], mins[4], maxs[4]); // late pulse mu
    minimizer.addParameter(seeds[5], grad_scales[5], mins[5], maxs[5]); // late pulse sigma
    minimizer.addParameter(0.0);                                        // late pulse scale
    minimizer.setHistorySize(20);

    // fix parameter idx of guys we are not minimizing
    minimizer.fixParameter(1); // Rt
    minimizer.fixParameter(3); // tau_t
    minimizer.fixParameter(4); // tau_rec
    minimizer.fixParameter(10); // time offset
    minimizer.fixParameter(11); // const offset
    minimizer.fixParameter(14); // LP scale 

    bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));

    if (succeeded){
        std::cout << "joint fit converged!" << std::endl;
        std::cout << "minimization finished with " << minimizer.errorMessage() << " error message" << std::endl; 
    }

    double value = minimizer.minimumValue();                                                                                                                                                 
    std::vector<double> params = minimizer.minimumPosition();
    
    // now let's find the smallest func_val!
    std::cout << "Function value at minimum: " << value << " after " << minimizer.numberOfEvaluations() << " function evaluations" << std::endl;
    std::cout << "Parameters: " << std::endl;
    std::vector<double> data_to_return;
    data_to_return.push_back(value);
    for(size_t i=0; i<LikelihoodFunctor::DerivativeDimension; ++i) {
        data_to_return.push_back(params.at(i));
        if (i == 1 or i == 3 or i == 4 or i == 10 or i == 11 or i == 14){
            continue;
        }
        std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
    }
    return data_to_return;
}


