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

struct ZigZagPrior{
    private:
        double point;
        double scale;
        double m;
    public:
        ZigZagPrior(double point, double scale, bool small):
        point(point), scale(scale), m(int(small)*2.0 - 1.0) {}

        template<typename DataType>
        DataType operator()(DataType x) const{
            return log((tanh(scale*(point-x)*m)+1.)/2.+1e-18)-exp(-scale*(point-x)*m);
        }
};

struct NewLikelihoodFunctor {

    static constexpr int NewDerivativeDimension = 7 + 200 + 1 + 1; // max number of dimensions :
                                                                   // Rs, tau_s, tau_TPB
                                                                   // norm, norm, norm, norm -- 1 / data set!
                                                                   // 200 PMT efficiency terms (many will be fixed)                                                 
                                                                   // finally -- uv absorption
                                                                   // jk, one more variable -- photons per mev!                                                               

                                                                   // things that are fixed for all PMTs:
                                                                   // Rt, tau_t, tau_rec, const_offset 
                                                                   
                                                                   // things that are fixed for each PMT : 
                                                                   // time offset (1 / PMT / data set), LPmu, LPsigma, and LPscale
    
    typedef double Underlying;
    typedef phys_tools::autodiff::FD<NewDerivativeDimension, Underlying> AD;
    
    std::vector<std::shared_ptr<CalculateNLLH<AD>>> llh_constructorAD;
    std::vector<std::shared_ptr<CalculateNLLH<double>>> llh_constructorDouble;
    
    I3VectorCCMPMTKey all_keys; // list of PMTs in fit
    I3MapPMTKeyDouble LPmu; // map between CCMPMTKey and LPmu
    I3MapPMTKeyDouble LPsigma; // map between CCMPMTKey and LPsigma
    I3MapPMTKeyDouble LPscale; // map between CCMPMTKey and LPscale
    std::vector<I3MapPMTKeyDouble> time_offsets; // map between CCMPMTKey and time offset (1 map/data set) 
    double Rt = 0.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    std::vector<double> z_offset; // list of z offsets
    size_t n_sodium_events;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
    ZigZagPrior prior = ZigZagPrior(3.0, 6.0, false); // we want (tau_s - tau_TPB) > 3.0 with 600% scale
    ZigZagPrior uv_abs_prior = ZigZagPrior(65.0, 6.0, true); // we want uv abs < 65.0 with 600% scale
    bool fitting_uv_abs = true;
    bool apply_uv_abs_prior = false;

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        T total_llh = 0;
        if constexpr (std::is_same<T, double>::value) {
            std::cout << "Rs = " << x[0] << ", tau_s = " << x[1] << ", tau_TPB = " << x[2] << ", norm = " << x[3] << ", " << x[4] << ", " << x[5] << ", " << x[6]
                << ", uv abs = " << x[207] << ", and photons/mev = " << x[208] << std::endl;
        } else {
            std::cout <<"Rs = " << x[0].value() << ", tau_s = " << x[1].value() << ", tau_TPB = " << x[2].value() << ", norm = " << x[3].value() << ", " << x[4].value()
                << ", " << x[5].value() << ", " << x[6].value() << ", uv abs = " << x[207].value() << ", and photons/mev = " << x[208].value() << std::endl;
        }
        for (size_t data_it = 0; data_it < llh_constructorAD.size(); data_it ++){
            // loop over each data set we are fitting to
            T tau_s = x[1];
            T tau_TPB = x[2];        
            for (size_t pmt_it = 0; pmt_it < all_keys.size(); pmt_it ++){
                // loop over each PMT
                CCMPMTKey key = all_keys.at(pmt_it);
                if constexpr (std::is_same<T, double>::value) {
                    total_llh += llh_constructorDouble.at(data_it)->ComputeNLLH(key, x[0], Rt, x[1], tau_t, tau_rec, x[2], x[3 + data_it], // normalization for data set!!!
                                    time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), x[7 + pmt_it], // pmt eff!!
                                    x[207], x[208], fitting_uv_abs, z_offset.at(data_it), n_sodium_events, light_profile_type);
                } else {
                    total_llh += llh_constructorAD.at(data_it)->ComputeNLLH(key, x[0], Rt, x[1], tau_t, tau_rec, x[2], x[3 + data_it], // normalization for data set!!!
                                    time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), x[7 + pmt_it], // pmt eff!
                                    x[207], x[208], fitting_uv_abs, z_offset.at(data_it), n_sodium_events, light_profile_type);
                }
                total_llh += prior(tau_s - tau_TPB); 
                if (apply_uv_abs_prior){
                    total_llh += uv_abs_prior(x[207]); 
                }
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
    
    double seeds[4] = {0.34, 8.2, 3.0, 500.0};
    double grad_scales[4] = {1e-3, 1e-3, 1e-3, 1e-3};
    double mins[4] = {0.2, 2.0, 1.0, 10.0};
    double maxs[4] = {0.5, 16.0, 9.0, 1e10};

    // let's initialize our constructor
    // grab our geomtry frame
    std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
    dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
    I3FramePtr geo_frame = geometry_file.pop_frame();
    
    // now set up our likelihood object
    LikelihoodType likelihood;
    
    // make our AD object
    typedef double Underlying;
    typedef phys_tools::autodiff::FD<likelihood.NewDerivativeDimension, Underlying> AD;
    
    // now grab data frame(s) to make our constructors
    std::vector<std::shared_ptr<CalculateNLLH<AD>>> all_constructorsAD;
    std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDouble;
    for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
        std::string this_data_fname = data_file_names.at(data_it);
        dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
        I3FramePtr this_data_frame = this_data_file.pop_frame();

        // and now make constructor -- AD type
        std::shared_ptr<CalculateNLLH<AD>> this_llh_constructorAD = std::make_shared<CalculateNLLH<AD>>();
        this_llh_constructorAD->SetKeys(keys_to_fit);
        this_llh_constructorAD->SetData(this_data_frame);
        this_llh_constructorAD->SetGeo(geo_frame);
    
        // now save!
        all_constructorsAD.push_back(this_llh_constructorAD);
        
        // and now make constructor -- double type
        std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
        this_llh_constructorDouble->SetKeys(keys_to_fit);
        this_llh_constructorDouble->SetData(this_data_frame);
        this_llh_constructorDouble->SetGeo(geo_frame);
    
        // now save!
        all_constructorsDouble.push_back(this_llh_constructorDouble);
    }
    
    likelihood.llh_constructorAD = all_constructorsAD;
    likelihood.llh_constructorDouble = all_constructorsDouble;
    likelihood.n_sodium_events = n_sodium_events;
    likelihood.z_offset = z_offsets; 
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
    // now we are adding our pmt efficincy terms
    size_t n_pmts_fitting = keys_to_fit.size();
    for (size_t n = 0; n < 200; n++){
        double pmt_eff_seed = 1.0;
        if (n < n_pmts_fitting){
           pmt_eff_seed = 1.0 / PMT_efficiencies[keys_to_fit.at(n)]; 
        }
        minimizer.addParameter(pmt_eff_seed, 1e-4, 1.0, pmt_eff_seed * 1e2); // pmt eff
    } 
    minimizer.addParameter(55.0, 1e-5, 30.0, 100.0); // uv absorption!!!
    minimizer.addParameter(1.0); // photons per mev!!!
    minimizer.setHistorySize(20);

    // fix parameter idx of guys we are not minimizing
    size_t data_sets_to_minimize = data_file_names.size();
    size_t data_sets_included = 0;
    for (size_t n = 0; n < 4; n++){
        data_sets_included += 1;
        if (data_sets_included > data_sets_to_minimize){
            minimizer.fixParameter(2 + data_sets_included); // norms
        }
    }

    // now fix PMT efficincy params for PMTs we are not fitting to
    size_t n_pmts_included = 0;
    for (size_t n = 0; n < 200; n++){
        n_pmts_included += 1;
        if (n_pmts_included > n_pmts_fitting){
            minimizer.fixParameter(6 + n_pmts_included); // pmt eff
        }
    } 

    // fix photons per mev
    minimizer.fixParameter(208);
    
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
            if (i < 7) {
                std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
            } else {
                std::cout << params.at(i) << std::endl;
            }
        }
    
        // one last thing -- take this best fit point and grab our data, pred, and times
        double Rt = 0.0;
        double tau_t = 743.0;
        double tau_rec = 0.0;
        double const_offset = 0.0;
        AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
        for (size_t data_it = 0; data_it < all_constructorsDouble.size(); data_it ++){
            // loop over each data set we are fitting to
            I3MapPMTKeyVectorDouble this_data; 
            I3MapPMTKeyVectorDouble this_pred; 
            I3MapPMTKeyVectorDouble this_times; 
            
            for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++){
                // loop over each PMT
                CCMPMTKey key = keys_to_fit.at(pmt_it);
                double llh = all_constructorsDouble.at(data_it)->ComputeNLLH(key, params.at(0), Rt, params.at(1), tau_t, tau_rec, params.at(2), params.at(3 + data_it), // normalization for data set!!!
                                 time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), params.at(7 + pmt_it), params.at(207), params.at(208), true,
                                 z_offsets.at(data_it), n_sodium_events, light_profile_type);
                // now grab out data etc
                std::vector<double> this_pmt_data = all_constructorsDouble.at(data_it)->GetDataVector(); 
                std::vector<double> this_pmt_pred = all_constructorsDouble.at(data_it)->GetPredVector(); 
                std::vector<double> this_pmt_times = all_constructorsDouble.at(data_it)->GetTimesVector(); 
                // now save
                this_data[key] = this_pmt_data;
                this_pred[key] = this_pmt_pred;
                this_times[key] = this_pmt_times;
            }
            data.push_back(this_data);
            pred.push_back(this_pred);
            times.push_back(this_times);
        }
    
    } else {
        std::cout << "minimizer did not converge :( error message = " << minimizer.errorMessage() << std::endl;
    }
    
    return data_to_return;

}

void AllPMTcppMinimizer::ScanOverTOffsets(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble PMT_efficiencies, I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                          std::vector<std::string> data_file_names, std::vector<double> z_offsets, size_t n_sodium_events,
                                          double best_fit_Rs, double best_fit_tau_s, double best_fit_tau_other, std::vector<double> best_fit_norms) {        

    // let's initialize our constructor
    // grab our geometry frame
    std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
    dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
    I3FramePtr geo_frame = geometry_file.pop_frame();
    
    // now grab data frame(s) to make our constructors
    std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDouble;
    for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
        std::string this_data_fname = data_file_names.at(data_it);
        dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
        I3FramePtr this_data_frame = this_data_file.pop_frame();

        // and now make constructor -- double type
        std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
        this_llh_constructorDouble->SetKeys(keys_to_fit);
        this_llh_constructorDouble->SetData(this_data_frame);
        this_llh_constructorDouble->SetGeo(geo_frame);
    
        // now save!
        all_constructorsDouble.push_back(this_llh_constructorDouble);
    }

    // now loop over each data set and then loop over each pmt
    double uv_absorption = 55.0;
    double Rt = 0.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    double photons_per_mev = 1.0;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
    
    for (size_t data_it = 0; data_it < all_constructorsDouble.size(); data_it ++){
        // loop over each data set we are fitting to
        I3MapPMTKeyDouble this_data_set_best_time_offsets;

        for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++){
            // loop over each PMT
            CCMPMTKey key = keys_to_fit.at(pmt_it);
            
            // let's try scanning over t offset
            double t_offset_start = 15.0;
            double t_offset_end = 40.0;
            double t_offset_range = t_offset_end - t_offset_start;
            size_t n_t_offsets = 100;

            std::vector<double> func_val;
            std::vector<double> time_offsets;
            
            for (size_t toff_it = 0; toff_it < (n_t_offsets + 1); toff_it++){
                double this_t_offset = t_offset_start + ((double)toff_it / (double)n_t_offsets) * t_offset_range;
                double this_nllh = -1.0 * all_constructorsDouble.at(data_it)->ComputeNLLH(key, best_fit_Rs, Rt, best_fit_tau_s, tau_t, tau_rec, best_fit_tau_other, best_fit_norms.at(data_it), 
                             this_t_offset, const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), PMT_efficiencies.at(key), uv_absorption, photons_per_mev, false, z_offsets.at(data_it),
                             n_sodium_events, light_profile_type);
            
                // now let's save this nllh and corresponding time offset
                func_val.push_back(this_nllh);
                time_offsets.push_back(this_t_offset);
            }

            // now let's find the time offset that minimizes our nllh
            size_t smallest_idx = std::distance(func_val.begin(), std::min_element(func_val.begin(), func_val.end()));
            this_data_set_best_time_offsets[key] = time_offsets.at(smallest_idx);         
        }

        best_time_offsets.push_back(this_data_set_best_time_offsets);
    }
}

std::vector<double> AllPMTcppMinimizer::GrabNormSeed(CCMPMTKey key, double baseline_efficiency, double LPmu, double LPsigma, double LPscale, 
                                                     double uv_abs, std::vector<double> z_offsets, size_t n_sodium_events,
                                                     std::vector<I3MapPMTKeyDouble> time_offsets, std::vector<std::shared_ptr<CalculateNLLH<double>>> llh_constructorDouble){
    
    // for each data set, going to grab expectation for benchmark PMT
    // then we can get a reasonable looking norm
    double Rs = 0.3;
    double Rt = 0.0;
    double tau_s = 8.0;
    double tau_TPB = 4.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    double photons_per_mev = 1.0;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
   
    std::vector<double> norm_seeds; 
    for (size_t data_it = 0; data_it < llh_constructorDouble.size(); data_it ++){
        // loop over each data set we are fitting to
        double llh = llh_constructorDouble.at(data_it)->ComputeNLLH(key, Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, 1e10, // normalization for data set!!!
                         time_offsets.at(data_it).at(key), const_offset, LPmu, LPsigma, LPscale, baseline_efficiency,
                         uv_abs, photons_per_mev, true, z_offsets.at(data_it), n_sodium_events, light_profile_type);
        
        // now grab out data etc
        std::vector<double> this_pmt_data = llh_constructorDouble.at(data_it)->GetDataVector(); 
        std::vector<double> this_pmt_pred = llh_constructorDouble.at(data_it)->GetPredVector(); 
        
        // now let's get norm scaling
        size_t data_max_idx = std::distance(this_pmt_data.begin(), std::max_element(this_pmt_data.begin(), this_pmt_data.end()));
        size_t pred_max_idx = std::distance(this_pmt_pred.begin(), std::max_element(this_pmt_pred.begin(), this_pmt_pred.end()));
        
        double norm_scaling = (this_pmt_data.at(data_max_idx) / this_pmt_pred.at(pred_max_idx)) * 1e10;
        
        // now save
        norm_seeds.push_back(norm_scaling);
    }

    return norm_seeds;
}

std::vector<double> AllPMTcppMinimizer::ScanOverUVAbsorption(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble PMT_efficiencies,
                                                             I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                                             I3MapPMTKeyDouble time_offset1, I3MapPMTKeyDouble time_offset2, I3MapPMTKeyDouble time_offset3, I3MapPMTKeyDouble time_offset4,
                                                             std::vector<std::string> data_file_names, std::vector<double> z_offsets, size_t n_sodium_events) {        
    
    std::vector<I3MapPMTKeyDouble> time_offsets = {time_offset1, time_offset2, time_offset3, time_offset4}; 
    std::vector<std::string> paramter_names = {"Rs", "tau_s", "tau_TPB", "norm1", "norm2", "norm3", "norm4"};

    // for this fit -- we are looping over reasonable uv absorption values and minimizing at each
    // then returning best fit across this scan
    
    // set up our grid of uv absorptions
    double uv_abs_start = 50.0;
    double uv_abs_end = 60.0;
    double uv_abs_range = uv_abs_end - uv_abs_start;
    size_t n_uv_abs = 10;

    // place to store results of best fit
    std::vector<double> func_val;
    std::vector<double> n_evals;
    std::vector<std::vector<double>> all_params;
    std::vector<std::vector<I3MapPMTKeyVectorDouble>> all_data;
    std::vector<std::vector<I3MapPMTKeyVectorDouble>> all_pred;
    std::vector<std::vector<I3MapPMTKeyVectorDouble>> all_times;
    
    double seeds[4] = {0.34, 8.2, 3.0, 1000.0};
    double grad_scales[4] = {1e-3, 1e-3, 1e-3, 1e-3};
    double mins[4] = {0.2, 2.0, 1.0, 10.0};
    double maxs[4] = {0.5, 16.0, 9.0, 1e10};

    // now we need to make a new constructor for each uv abs
    // so time to loop over them
    for (size_t uvab_it = 0; uvab_it < (n_uv_abs + 1); uvab_it++){
        double this_uv_abs = uv_abs_start + ((double)uvab_it / (double)n_uv_abs) * uv_abs_range;
        bool this_uv_success = false;
        while (this_uv_success == false){        
            // let's initialize our constructor
            // grab our geomtry frame
            std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
            dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
            I3FramePtr geo_frame = geometry_file.pop_frame();
            
            // now set up our likelihood object
            LikelihoodType likelihood;
            
            // make our AD object
            typedef double Underlying;
            typedef phys_tools::autodiff::FD<likelihood.NewDerivativeDimension, Underlying> AD;
            
            // now grab data frame(s) to make our constructors
            std::vector<std::shared_ptr<CalculateNLLH<AD>>> all_constructorsAD;
            std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDouble;
            for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
                std::string this_data_fname = data_file_names.at(data_it);
                dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
                I3FramePtr this_data_frame = this_data_file.pop_frame();

                // and now make constructor -- AD type
                std::shared_ptr<CalculateNLLH<AD>> this_llh_constructorAD = std::make_shared<CalculateNLLH<AD>>();
                this_llh_constructorAD->SetKeys(keys_to_fit);
                this_llh_constructorAD->SetData(this_data_frame);
                this_llh_constructorAD->SetGeo(geo_frame);
            
                // now save!
                all_constructorsAD.push_back(this_llh_constructorAD);
                
                // and now make constructor -- double type
                std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
                this_llh_constructorDouble->SetKeys(keys_to_fit);
                this_llh_constructorDouble->SetData(this_data_frame);
                this_llh_constructorDouble->SetGeo(geo_frame);
            
                // now save!
                all_constructorsDouble.push_back(this_llh_constructorDouble);
            }
            
            likelihood.llh_constructorAD = all_constructorsAD;
            likelihood.llh_constructorDouble = all_constructorsDouble;
            likelihood.n_sodium_events = n_sodium_events;
            likelihood.z_offset = z_offsets; 
            likelihood.all_keys = keys_to_fit;
            likelihood.LPmu = LPmu;
            likelihood.LPsigma = LPsigma;
            likelihood.LPscale = LPscale;
            likelihood.time_offsets = time_offsets;
            likelihood.fitting_uv_abs = false;
           
            // now let's grab our norm seeds
            size_t reference_pmt_idx = (size_t) (keys_to_fit.size() / 2);
            CCMPMTKey reference_pmt = keys_to_fit.at(reference_pmt_idx);
            
            std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDoubleNormSeed;
            for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
                std::string this_data_fname = data_file_names.at(data_it);
                dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
                I3FramePtr this_data_frame = this_data_file.pop_frame();

                // and now make constructor -- double type
                std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
                this_llh_constructorDouble->SetKeys(I3VectorCCMPMTKey({reference_pmt}));
                this_llh_constructorDouble->SetData(this_data_frame);
                this_llh_constructorDouble->SetGeo(geo_frame);
            
                // now save!
                all_constructorsDoubleNormSeed.push_back(this_llh_constructorDouble);
            }
            std::vector<double> norm_seeds = GrabNormSeed(reference_pmt, 0.5, LPmu.at(reference_pmt), LPsigma.at(reference_pmt), LPscale.at(reference_pmt),
                                                         this_uv_abs, z_offsets, n_sodium_events, time_offsets, all_constructorsDouble);

            std::cout << "for uv abs = " << this_uv_abs << " norm seeds = " << norm_seeds << std::endl;
            // now set up our minimizer
            phys_tools::lbfgsb::LBFGSB_Driver minimizer;
                
            minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]); // Rs
            minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]); // tau_s
            minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]); // tau_TPB
            // this is where we add our normalization parameter...need to add one for each data set
            size_t total_data_sets = all_constructorsDouble.size();
            size_t data_sets_accounted_for = 0;
            for (size_t n = 0; n < 4; n++){
                data_sets_accounted_for += 1;
                if (data_sets_accounted_for <= total_data_sets) {
                    minimizer.addParameter(norm_seeds.at(n), grad_scales[3], mins[3], maxs[3]); // norm
                }
                else {
                    minimizer.addParameter(seeds[3], grad_scales[3], mins[3], maxs[3]); // norm
                }
            }

            // now we are adding our pmt efficincy terms
            size_t n_pmts_fitting = keys_to_fit.size();
            for (size_t n = 0; n < 200; n++){
                double pmt_eff_seed = 1.0;
                if (n < n_pmts_fitting){
                   pmt_eff_seed = 1.0 / PMT_efficiencies[keys_to_fit.at(n)]; 
                }
                //minimizer.addParameter(pmt_eff_seed, 1e-4, 1.0, pmt_eff_seed * 1e2); // pmt eff
                minimizer.addParameter(0.5, 1e-3, 0.1, 10.0); // pmt eff
            } 
            minimizer.addParameter(this_uv_abs, 1e-3, 30.0, 100.0); // uv absorption!!!
            minimizer.addParameter(1.0); // photons per mev!!!
            minimizer.setHistorySize(20);

            // fix parameter idx of guys we are not minimizing
            size_t data_sets_to_minimize = data_file_names.size();
            size_t data_sets_included = 0;
            for (size_t n = 0; n < 4; n++){
                data_sets_included += 1;
                if (data_sets_included > data_sets_to_minimize){
                    minimizer.fixParameter(2 + data_sets_included); // norms
                }
            }

            // now fix PMT efficincy params for PMTs we are not fitting to
            size_t n_pmts_included = 0;
            for (size_t n = 0; n < 200; n++){
                n_pmts_included += 1;
                if (n_pmts_included > n_pmts_fitting){
                    minimizer.fixParameter(6 + n_pmts_included); // pmt eff
                }
            }

            // finally, fix uv absorption since we are not minimizing but scanning over
            minimizer.fixParameter(207); // uv absorption!!! 
            minimizer.fixParameter(208); // photons per mev!!! 

            bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
            
            if(succeeded) {
                this_uv_success = true;
                std::cout << "fit converged for uv absorption = " << this_uv_abs << "!" << std::endl;
                std::cout << "minimization finished with " << minimizer.errorMessage() << " error message" << std::endl; 
                
                // save things for returning!
                double value = minimizer.minimumValue();                                                                                                                                                 
                std::vector<double> params = minimizer.minimumPosition();
                
                std::cout << "Function value at minimum: " << value << " after " << minimizer.numberOfEvaluations() << " function evaluations" << std::endl;
                std::cout << "Parameters: " << std::endl;
                func_val.push_back(value);
                n_evals.push_back((double) minimizer.numberOfEvaluations()); 
                std::vector<double> this_fit_params;
                for(size_t i=0; i<NewLikelihoodFunctor::NewDerivativeDimension; ++i) {
                    this_fit_params.push_back(params.at(i));
                    if (i < 7) {
                        std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
                    } else {
                        std::cout << params.at(i) << std::endl;
                    }
                }
                all_params.push_back(this_fit_params);
            
                // one last thing -- take this best fit point and grab our data, pred, and times
                double Rt = 0.0;
                double tau_t = 743.0;
                double tau_rec = 0.0;
                double const_offset = 0.0;
                AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
                
                std::vector<I3MapPMTKeyVectorDouble> this_fit_data;
                std::vector<I3MapPMTKeyVectorDouble> this_fit_pred;
                std::vector<I3MapPMTKeyVectorDouble> this_fit_times;
                for (size_t data_it = 0; data_it < all_constructorsDouble.size(); data_it ++){
                    // loop over each data set we are fitting to
                    I3MapPMTKeyVectorDouble this_data; 
                    I3MapPMTKeyVectorDouble this_pred; 
                    I3MapPMTKeyVectorDouble this_times; 
                    
                    for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++){
                        // loop over each PMT
                        CCMPMTKey key = keys_to_fit.at(pmt_it);
                        double llh = all_constructorsDouble.at(data_it)->ComputeNLLH(key, params.at(0), Rt, params.at(1), tau_t, tau_rec, params.at(2), params.at(3 + data_it), // normalization for data set!!!
                                         time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), params.at(7 + pmt_it), params.at(207),
                                         params.at(208), false, z_offsets.at(data_it), n_sodium_events, light_profile_type);
                        // now grab out data etc
                        std::vector<double> this_pmt_data = all_constructorsDouble.at(data_it)->GetDataVector(); 
                        std::vector<double> this_pmt_pred = all_constructorsDouble.at(data_it)->GetPredVector(); 
                        std::vector<double> this_pmt_times = all_constructorsDouble.at(data_it)->GetTimesVector(); 
                        // now save
                        this_data[key] = this_pmt_data;
                        this_pred[key] = this_pmt_pred;
                        this_times[key] = this_pmt_times;
                    }
                    this_fit_data.push_back(this_data);
                    this_fit_pred.push_back(this_pred);
                    this_fit_times.push_back(this_times);
                }
                all_data.push_back(this_fit_data);
                all_pred.push_back(this_fit_pred);
                all_times.push_back(this_fit_times);
            } else {
            std::cout << "minimizer did not converge for uv absorption = " << this_uv_abs << " :( error message = " << minimizer.errorMessage() << std::endl;
            }
        }
    }

    // now that we are done scanning over uv absorption, let's get the min!
    size_t smallest_idx = std::distance(func_val.begin(), std::min_element(func_val.begin(), func_val.end()));
    std::vector<double> data_to_return;
    data_to_return.push_back(func_val.at(smallest_idx));
    data_to_return.push_back(n_evals.at(smallest_idx));
    
    // now save best fit params
    std::vector<double> best_fit_params = all_params.at(smallest_idx);
    for (size_t i = 0; i < best_fit_params.size(); i++){
        data_to_return.push_back(best_fit_params.at(i));
    }

    // and now set data, pred and times equal to best fit
    data = all_data.at(smallest_idx);
    pred = all_pred.at(smallest_idx);
    times = all_times.at(smallest_idx);

    return data_to_return;

}

std::vector<double> AllPMTcppMinimizer::GrabPhotonsPerMeVSeed(CCMPMTKey key, double baseline_efficiency, double LPmu, double LPsigma, double LPscale, 
                                                              double uv_abs, std::vector<double> z_offsets, size_t n_sodium_events,
                                                              std::vector<I3MapPMTKeyDouble> time_offsets, std::vector<std::shared_ptr<CalculateNLLH<double>>> llh_constructorDouble){
    
    // for each data set, going to grab expectation for benchmark PMT
    // then we can get a reasonable looking norm
    double Rs = 0.3;
    double Rt = 0.0;
    double tau_s = 8.0;
    double tau_TPB = 4.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    double photons_per_mev = 1e10;
    double norm = 1.0;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
   
    std::vector<double> photons_per_mev_seeds; 
    for (size_t data_it = 0; data_it < llh_constructorDouble.size(); data_it ++){
        // loop over each data set we are fitting to
        double llh = llh_constructorDouble.at(data_it)->ComputeNLLH(key, Rs, Rt, tau_s, tau_t, tau_rec, tau_TPB, norm,
                         time_offsets.at(data_it).at(key), const_offset, LPmu, LPsigma, LPscale, baseline_efficiency,
                         uv_abs, photons_per_mev, true, z_offsets.at(data_it), n_sodium_events, light_profile_type);
        
        // now grab out data etc
        std::vector<double> this_pmt_data = llh_constructorDouble.at(data_it)->GetDataVector(); 
        std::vector<double> this_pmt_pred = llh_constructorDouble.at(data_it)->GetPredVector(); 
        
        // now let's get norm scaling
        size_t data_max_idx = std::distance(this_pmt_data.begin(), std::max_element(this_pmt_data.begin(), this_pmt_data.end()));
        size_t pred_max_idx = std::distance(this_pmt_pred.begin(), std::max_element(this_pmt_pred.begin(), this_pmt_pred.end()));
        
        double photons_per_mev_scaling = (this_pmt_data.at(data_max_idx) / this_pmt_pred.at(pred_max_idx)) * photons_per_mev;
        
        // now save
        photons_per_mev_seeds.push_back(photons_per_mev_scaling);
    }

    return photons_per_mev_seeds;
}


std::vector<double> AllPMTcppMinimizer::FitUVAbsorption(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                                        I3MapPMTKeyDouble time_offset1, I3MapPMTKeyDouble time_offset2, I3MapPMTKeyDouble time_offset3, I3MapPMTKeyDouble time_offset4,
                                                        std::vector<std::string> data_file_names, std::vector<double> z_offsets, size_t n_sodium_events, std::vector<double> best_fit_params,
                                                        bool apply_uv_abs_prior){

    std::vector<I3MapPMTKeyDouble> time_offsets = {time_offset1, time_offset2, time_offset3, time_offset4}; 
    std::vector<std::string> paramter_names = {"Rs", "tau_s", "tau_TPB", "norm1", "norm2", "norm3", "norm4"};

    bool uv_fit_success = false;
    std::vector<double> data_to_return;
    
    while (uv_fit_success == false){        
        
        // let's initialize our constructor
        // grab our geomtry frame
        std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
        dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
        I3FramePtr geo_frame = geometry_file.pop_frame();
        
        // now set up our likelihood object
        LikelihoodType likelihood;
        
        // make our AD object
        typedef double Underlying;
        typedef phys_tools::autodiff::FD<likelihood.NewDerivativeDimension, Underlying> AD;
        
        // now grab data frame(s) to make our constructors
        std::vector<std::shared_ptr<CalculateNLLH<AD>>> all_constructorsAD;
        std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDouble;
        for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
            std::string this_data_fname = data_file_names.at(data_it);
            dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
            I3FramePtr this_data_frame = this_data_file.pop_frame();

            // and now make constructor -- AD type
            std::shared_ptr<CalculateNLLH<AD>> this_llh_constructorAD = std::make_shared<CalculateNLLH<AD>>();
            this_llh_constructorAD->SetKeys(keys_to_fit);
            this_llh_constructorAD->SetData(this_data_frame);
            this_llh_constructorAD->SetGeo(geo_frame);
        
            // now save!
            all_constructorsAD.push_back(this_llh_constructorAD);
            
            // and now make constructor -- double type
            std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
            this_llh_constructorDouble->SetKeys(keys_to_fit);
            this_llh_constructorDouble->SetData(this_data_frame);
            this_llh_constructorDouble->SetGeo(geo_frame);
        
            // now save!
            all_constructorsDouble.push_back(this_llh_constructorDouble);
        }
        
        likelihood.llh_constructorAD = all_constructorsAD;
        likelihood.llh_constructorDouble = all_constructorsDouble;
        likelihood.n_sodium_events = n_sodium_events;
        likelihood.z_offset = z_offsets; 
        likelihood.all_keys = keys_to_fit;
        likelihood.LPmu = LPmu;
        likelihood.LPsigma = LPsigma;
        likelihood.LPscale = LPscale;
        likelihood.time_offsets = time_offsets;
        likelihood.fitting_uv_abs = true; 
        likelihood.apply_uv_abs_prior = apply_uv_abs_prior;
    
        // now set up our minimizer
        phys_tools::lbfgsb::LBFGSB_Driver minimizer;
            
        minimizer.addParameter(best_fit_params.at(0 + 2), 1e-3, 0.1, 0.5); // Rs
        minimizer.addParameter(best_fit_params.at(1 + 2), 1e-3, 2.0, 16.0); // tau_s
        minimizer.addParameter(best_fit_params.at(2 + 2), 1e-3, 1.0, 9.0); // tau_TPB
        
        // this is where we add our normalization parameter...need to add one for each data set
        size_t total_data_sets = all_constructorsDouble.size();
        size_t data_sets_accounted_for = 0;
        for (size_t n = 0; n < 4; n++){
            //data_sets_accounted_for += 1;
            //if (data_sets_accounted_for <= total_data_sets) {
            //    double norm_seed = best_fit_params.at(3 + n + 2);
            //    minimizer.addParameter(norm_seed, 1e-3, norm_seed * 1e-4, norm_seed * 1e4); // norm
            //}
            //else {
            //    minimizer.addParameter(1000.0, 1e-3, 10.0, 1e10); // norm
            //}
            minimizer.addParameter(1.0); // norm
        }
        
        // now we are adding our pmt efficincy terms
        for (size_t n = 0; n < 200; n++){
            //minimizer.addParameter(best_fit_params.at(7 + n + 2), 1e-3, 0.1, 10.0); // pmt eff
            minimizer.addParameter(0.5, 1e-3, 0.01, 10.0); // pmt eff
        } 
        double uv_abs_seed = 65.0;        
        minimizer.addParameter(uv_abs_seed, 1e-3, 20.0, 300.0); // uv absorption!!!
        
        // now let's grab our photons/mev seeds
        size_t reference_pmt_idx = (size_t) (keys_to_fit.size() / 2);
        CCMPMTKey reference_pmt = keys_to_fit.at(reference_pmt_idx);
        
        std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDoubleNormSeed;
        for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
            std::string this_data_fname = data_file_names.at(data_it);
            dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
            I3FramePtr this_data_frame = this_data_file.pop_frame();

            // and now make constructor -- double type
            std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
            this_llh_constructorDouble->SetKeys(I3VectorCCMPMTKey({reference_pmt}));
            this_llh_constructorDouble->SetData(this_data_frame);
            this_llh_constructorDouble->SetGeo(geo_frame);
        
            // now save!
            all_constructorsDoubleNormSeed.push_back(this_llh_constructorDouble);
        }
        std::vector<double> photons_per_mev_seeds = GrabPhotonsPerMeVSeed(reference_pmt, 0.5, LPmu.at(reference_pmt), LPsigma.at(reference_pmt), LPscale.at(reference_pmt),
                                                     uv_abs_seed, z_offsets, n_sodium_events, time_offsets, all_constructorsDouble);
        
        std::cout << "photons/mev seed = " << photons_per_mev_seeds << std::endl;
        
        // let's use the average photons/mev seed for our seed
        //double average_photons_per_mev = 0.0;
        //for (size_t p = 0; p < photons_per_mev_seeds.size(); p++){
        //    average_photons_per_mev += photons_per_mev_seeds.at(p);
        //}
        //average_photons_per_mev /= (double)photons_per_mev_seeds.size();

        minimizer.addParameter(photons_per_mev_seeds.at(0), 1e-3, photons_per_mev_seeds.at(0) * 1e-5, photons_per_mev_seeds.at(0) * 1e5); // photons per mev!!!
        minimizer.setHistorySize(20);
    
        // fix parameter idx of guys we are not minimizing
        size_t data_sets_to_minimize = data_file_names.size();
        size_t data_sets_included = 0;
        for (size_t n = 0; n < 4; n++){
            //data_sets_included += 1;
            //if (data_sets_included > data_sets_to_minimize){
            //    minimizer.fixParameter(2 + data_sets_included); // norms
            //}
            minimizer.fixParameter(3 + n); // norms
        }
        
        // now fix PMT efficincy params for PMTs we are not fitting to
        size_t n_pmts_included = 0;
        size_t n_pmts_fitting = keys_to_fit.size();
        for (size_t n = 0; n < 200; n++){
            //n_pmts_included += 1;
            //if (n_pmts_included > n_pmts_fitting){
            //    minimizer.fixParameter(6 + n_pmts_included); // pmt eff
            //}
            minimizer.fixParameter(7 + n); // pmt eff
        }
        
        bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
     
        if(succeeded) {
            uv_fit_success = true;
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
                if (i < 7) {
                    std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
                } else {
                    std::cout << params.at(i) << std::endl;
                }
            }
        
            // one last thing -- take this best fit point and grab our data, pred, and times
            double Rt = 0.0;
            double tau_t = 743.0;
            double tau_rec = 0.0;
            double const_offset = 0.0;
            AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
            for (size_t data_it = 0; data_it < all_constructorsDouble.size(); data_it ++){
                // loop over each data set we are fitting to
                I3MapPMTKeyVectorDouble this_data; 
                I3MapPMTKeyVectorDouble this_pred; 
                I3MapPMTKeyVectorDouble this_times; 
                
                for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++){
                    // loop over each PMT
                    CCMPMTKey key = keys_to_fit.at(pmt_it);
                    double llh = all_constructorsDouble.at(data_it)->ComputeNLLH(key, params.at(0), Rt, params.at(1), tau_t, tau_rec, params.at(2), params.at(3 + data_it), // normalization for data set!!!
                                     time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), params.at(7 + pmt_it), params.at(207),
                                     params.at(208), true, z_offsets.at(data_it), n_sodium_events, light_profile_type);
                    // now grab out data etc
                    std::vector<double> this_pmt_data = all_constructorsDouble.at(data_it)->GetDataVector(); 
                    std::vector<double> this_pmt_pred = all_constructorsDouble.at(data_it)->GetPredVector(); 
                    std::vector<double> this_pmt_times = all_constructorsDouble.at(data_it)->GetTimesVector(); 
                    // now save
                    this_data[key] = this_pmt_data;
                    this_pred[key] = this_pmt_pred;
                    this_times[key] = this_pmt_times;
                }
                data.push_back(this_data);
                pred.push_back(this_pred);
                times.push_back(this_times);
            }
        
        } else {
            std::cout << "minimizer did not converge :( error message = " << minimizer.errorMessage() << std::endl;
        }
        
    }    

    return data_to_return;

}


