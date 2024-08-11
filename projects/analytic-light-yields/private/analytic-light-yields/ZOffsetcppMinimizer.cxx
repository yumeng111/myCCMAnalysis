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
#include "analytic-light-yields/ZOffsetcppMinimizer.h"
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

struct ZOffsetFitLikelihoodFunctor {

    static constexpr int NewDerivativeDimension = 3 + 1 + 2 + 1 + 4 + 200; // max number of dimensions :
                                                                           // Rs, tau_s, tau_other
                                                                           // norm (for all data sets)
                                                                           // uv absorption 1, uv absorption 2
                                                                           // rayleigh scattering length!
                                                                           // z offset, z offset z offset, z offset
                                                                           // 200 PMT efficiencies

    typedef double Underlying;
    typedef phys_tools::autodiff::FD<NewDerivativeDimension, Underlying> AD;
    std::vector<std::shared_ptr<CalculateNLLH<AD>>> llh_constructorAD;

    I3VectorCCMPMTKey all_keys; // list of PMTs in fit
    I3MapPMTKeyDouble LPmu; // map between CCMPMTKey and LPmu
    I3MapPMTKeyDouble LPsigma; // map between CCMPMTKey and LPsigma
    I3MapPMTKeyDouble LPscale; // map between CCMPMTKey and LPscale
    std::vector<I3MapPMTKeyDouble> time_offsets; // map between CCMPMTKey and time offset (1 map/data set) 
    double Rt = 0.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    double photons_per_mev = 1.0;
    std::vector<size_t> n_sodium_events;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
    ZigZagPrior prior = ZigZagPrior(3.0, 6.0, false); // we want (tau_s - tau_TPB) > 3.0 with 600% scale
    bool use_g4_yields = true;
    bool use_tau_prior = true;

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        T total_llh = 0;
        T Rs = x[0];
        T tau_s = x[1];
        T tau_other = x[2];
        T norm = x[3];
        T uv_abs_1 = x[4];
        T uv_abs_2 = x[5];
        T rayl_length = x[6];

        std::cout << "Rs = " << Rs.value() << ", tau_s = " << tau_s.value() << ", tau_other = " << tau_other.value() << ", norm = " << norm.value()
            << ", uv abs 1 = " << uv_abs_1.value() << ", uv abs 2 = " << uv_abs_2.value() << ", rayl = " << rayl_length.value() << ", and z offsets = " << x[7].value()
            << ", " << x[8].value() << ", " << x[9].value() << ", " << x[10].value() << std::endl;

        for (size_t data_it = 0; data_it < llh_constructorAD.size(); data_it ++){
            // loop over each data set we are fitting to
            T this_z_offset = x[7 + data_it];

            for (size_t pmt_it = 0; pmt_it < all_keys.size(); pmt_it ++){
                // loop over each PMT
                CCMPMTKey key = all_keys.at(pmt_it);
                total_llh += llh_constructorAD.at(data_it)->ComputeNLLH(key, Rs, Rt, tau_s, tau_t, tau_rec, tau_other, norm,
                                time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), x[11 + pmt_it], // pmt eff!
                                {uv_abs_1, uv_abs_2}, rayl_length, photons_per_mev, this_z_offset, n_sodium_events.at(data_it), light_profile_type, use_g4_yields);
                if (use_tau_prior){
                    total_llh += prior(tau_s - tau_other); 
                }
            }
        }
        //std::cout << "total llh = " << total_llh << std::endl;
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
        return 0.0; // not using!!! plz dont need
        //return(-func.template evaluateLikelihood<double>(x)); // note negation!
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

typedef ZOffsetFitLikelihoodFunctor LikelihoodType;

ZOffsetcppMinimizer::ZOffsetcppMinimizer() {}

std::vector<double> ZOffsetcppMinimizer::GrabNormSeed(CCMPMTKey key, double baseline_efficiency, double LPmu, double LPsigma, double LPscale, 
                                                     std::vector<double> uv_abs_lengths, double rayl, std::vector<double> z_offsets, std::vector<size_t> n_sodium_events,
                                                     std::vector<I3MapPMTKeyDouble> time_offsets, std::vector<std::string> data_file_names){ 
    
    // for each data set, going to grab expectation for benchmark PMT
    // then we can get a reasonable looking norm
    
    // let's initialize our constructor
    // grab our geomtry frame
    std::string geometry_fname = "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/geo_run012490.i3.zst";
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
        this_llh_constructorDouble->SetKeys(I3VectorCCMPMTKey({key}));
        this_llh_constructorDouble->SetData(this_data_frame);

        // now save!
        all_constructorsDouble.push_back(this_llh_constructorDouble);
    }
   
    double Rs = 0.3;
    double Rt = 0.0;
    double tau_s = 8.0;
    double tau_other = 4.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    double photons_per_mev = 1.0;
    bool use_g4_yields = true;
    double norm_seed = 1e12;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 

    std::vector<double> norm_seeds; 
    for (size_t data_it = 0; data_it < all_constructorsDouble.size(); data_it ++){
        // loop over each data set we are fitting to
        double llh = all_constructorsDouble.at(data_it)->ComputeNLLH(key, Rs, Rt, tau_s, tau_t, tau_rec, tau_other, norm_seed, 
                         time_offsets.at(data_it).at(key), const_offset, LPmu, LPsigma, LPscale, baseline_efficiency,
                         uv_abs_lengths, rayl, photons_per_mev, z_offsets.at(data_it), n_sodium_events.at(data_it), light_profile_type, use_g4_yields);
        
        // now grab out data etc
        std::vector<double> this_pmt_data = all_constructorsDouble.at(data_it)->GetDataVector(); 
        std::vector<double> this_pmt_pred = all_constructorsDouble.at(data_it)->GetPredVector(); 
        
        // now let's get norm scaling
        size_t data_max_idx = std::distance(this_pmt_data.begin(), std::max_element(this_pmt_data.begin(), this_pmt_data.end()));
        size_t pred_max_idx = std::distance(this_pmt_pred.begin(), std::max_element(this_pmt_pred.begin(), this_pmt_pred.end()));
        
        double norm_scaling = (this_pmt_data.at(data_max_idx) / this_pmt_pred.at(pred_max_idx)) * norm_seed;
        
        // now save
        norm_seeds.push_back(norm_scaling);
    }

    return norm_seeds;
}


std::vector<double> ZOffsetcppMinimizer::FitParameters(I3VectorCCMPMTKey keys_to_fit, I3MapPMTKeyDouble LPmu, I3MapPMTKeyDouble LPsigma, I3MapPMTKeyDouble LPscale,
                                                       I3MapPMTKeyDouble time_offset1, I3MapPMTKeyDouble time_offset2, I3MapPMTKeyDouble time_offset3, I3MapPMTKeyDouble time_offset4,
                                                       std::vector<std::string> data_file_names, std::vector<double> z_offsets, std::vector<size_t> n_sodium_events,
                                                       I3MapPMTKeyDouble pmt_efficiency, std::vector<bool> fit_flags){

    std::cout << "fitting to " << data_file_names.size() << " data files" << std::endl;
    bool use_tau_prior = fit_flags.at(0);
    bool fix_z_offset = fit_flags.at(1);
    bool fix_pmt_eff = fit_flags.at(2);
    bool fix_second_abs_length = fit_flags.at(3);
    bool fix_rayl = fit_flags.at(4);

    std::vector<I3MapPMTKeyDouble> time_offsets = {time_offset1, time_offset2, time_offset3, time_offset4}; 
    std::vector<std::string> paramter_names = {"Rs", "tau_s", "tau_other", "norm", "uv_abs_1", "uv_abs_2", "rayl"};

    std::vector<double> data_to_return;

    // let's initialize our constructor
    // grab our geomtry frame
    std::string geometry_fname = "/lustre/scratch4/turquoise/darcyn/geant4_sodium_selected_events/geo_run012490.i3.zst";
    dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
    I3FramePtr geo_frame = geometry_file.pop_frame();

    // now set up our likelihood object
    LikelihoodType likelihood;
    typedef double Underlying;
    typedef phys_tools::autodiff::FD<ZOffsetFitLikelihoodFunctor::NewDerivativeDimension, Underlying> AD;

    // now grab data frame(s) to make our constructors
    std::vector<std::shared_ptr<CalculateNLLH<AD>>> all_constructorsAD;
    for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
        std::string this_data_fname = data_file_names.at(data_it);
        dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
        I3FramePtr this_data_frame = this_data_file.pop_frame();

        // and now make constructor -- AD type
        std::shared_ptr<CalculateNLLH<AD>> this_llh_constructorAD = std::make_shared<CalculateNLLH<AD>>();
        this_llh_constructorAD->SetKeys(keys_to_fit);
        this_llh_constructorAD->SetData(this_data_frame);

        // now save!
        all_constructorsAD.push_back(this_llh_constructorAD);
    }

    likelihood.llh_constructorAD = all_constructorsAD;
    likelihood.n_sodium_events = n_sodium_events;
    likelihood.all_keys = keys_to_fit;
    likelihood.LPmu = LPmu;
    likelihood.LPsigma = LPsigma;
    likelihood.LPscale = LPscale;
    likelihood.time_offsets = time_offsets;
    likelihood.use_tau_prior = use_tau_prior;

    // now set up our minimizer
    phys_tools::lbfgsb::LBFGSB_Driver minimizer;

    minimizer.addParameter(0.33, 1e-3, 0.1, 0.5); // Rs
    minimizer.addParameter(7.0, 1e-3, 1.0, 16.0); // tau_s
    minimizer.addParameter(3.0, 1e-3, 0.001, 9.0); // tau_TPB

    // let's grab the norm seed
    double uv_abs_1_seed = 45.0;
    double uv_abs_2_seed = 70.0;
    if (fix_second_abs_length){
        uv_abs_2_seed = 0.0;
    }
    double rayl_seed = 84.0;
    size_t reference_pmt_idx = (size_t) (keys_to_fit.size() / 2);
    CCMPMTKey reference_pmt = keys_to_fit.at(reference_pmt_idx);
    double baseline_eff = 1.0;
    if (pmt_efficiency.find(reference_pmt) != pmt_efficiency.end()){
        baseline_eff = pmt_efficiency.at(reference_pmt);
    }
    //std::vector<double> norm_seeds = GrabNormSeed(reference_pmt, baseline_eff, LPmu.at(reference_pmt), LPsigma.at(reference_pmt), LPscale.at(reference_pmt),
    //                                              {uv_abs_1_seed, uv_abs_2_seed}, rayl_seed, z_offsets, n_sodium_events, time_offsets, data_file_names);
    //double total_norm_seed = 0.0;
    //for (size_t n = 0; n < norm_seeds.size(); n++){
    //    total_norm_seed += norm_seeds.at(n);
    //}
    //double average_norm_seed = total_norm_seed / static_cast<double>(norm_seeds.size());
    
    //std::cout << "norm seeds = " << norm_seeds << " and adding average norm seed = " << average_norm_seed << std::endl;
    double average_norm_seed = 1.0;

    // now add normalization
    minimizer.addParameter(average_norm_seed, 1e-3, average_norm_seed * 1e-5, average_norm_seed * 1e5); // norm
   
    // now add uv abs
    minimizer.addParameter(uv_abs_1_seed, 1e-4, 5.0, 2800.0); // uv absorption 1!!!
    if (fix_second_abs_length){
        minimizer.addParameter(uv_abs_2_seed); // second uv abs length, not using, dont add bounds
    } else {
        minimizer.addParameter(uv_abs_2_seed, 1e-4, 5.0, 2800.0); // uv absorption 2!!! using!
    }
    
    // now add rayleigh scattering length
    minimizer.addParameter(rayl_seed, 1e-3, 55.01, 94.99);
    
    // now add z offsets
    for (size_t n = 0; n < 4; n++){
        if (n < data_file_names.size()){
            minimizer.addParameter(z_offsets.at(n), 1e-3, z_offsets.at(n) - 1.99, z_offsets.at(n) + 1.99);
        } else {
            minimizer.addParameter(0.0);
        }
    }
    
    // now we are adding our pmt efficincy terms
    for (size_t n = 0; n < 200; n++){
        if (n < keys_to_fit.size()){
            CCMPMTKey this_key = keys_to_fit.at(n);
            if (pmt_efficiency.find(this_key) != pmt_efficiency.end()){
                minimizer.addParameter(pmt_efficiency.at(this_key), 1e-3, pmt_efficiency.at(this_key) * 1e-2, pmt_efficiency.at(this_key) * 1e2); // pmt eff
            } 
            else {
                minimizer.addParameter(1.0, 1e-3, 0.001, 5.0); // pmt eff
            }
        }
        else {
            minimizer.addParameter(1.0, 1e-3, 0.001, 5.0); // pmt eff
        }
        //minimizer.addParameter(baseline_eff, 1e-3, 0.001, 5.0); // pmt eff
    } 
    
    minimizer.setHistorySize(20);

    // now fix some parameters
    // if we are only fitting for 1 uv absorption, fix the second one
    if (fix_second_abs_length){
        minimizer.fixParameter(5);
    }

    // if fix_rayl, fix rayl scattering length
    if (fix_rayl){
        minimizer.fixParameter(6);
    }

    // if fix_z_offset, fix ALL z offsets
    // otherwise only fix z offsets for data sets not in our fit
    if (fix_z_offset){
        for (size_t z = 0; z < 4; z++){
            minimizer.fixParameter(7 + z);
        }
    } else {
        for (size_t z = 0; z < 4; z++){
            if (z >= data_file_names.size()){
                minimizer.fixParameter(7 + z);
            }
        }
    }
    
    // same thing for pmt effiency
    if (fix_pmt_eff){
        for (size_t p = 0; p < 200; p++){
            minimizer.fixParameter(11 + p);    
        }
    } else{
        for (size_t p = 0; p < 200; p++){
            if (p >= keys_to_fit.size()){
                minimizer.fixParameter(11 + p);    
            }
        }
    }

    std::cout << "effective number of variables = " << minimizer.GeteffectivenVar() << std::endl;

    bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
 
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
        for(size_t i=0; i<ZOffsetFitLikelihoodFunctor::NewDerivativeDimension; ++i) {
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
        double photons_per_mev = 1.0;
        bool use_g4_yields = true;
        AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 
        for (size_t data_it = 0; data_it < all_constructorsAD.size(); data_it ++){
            // loop over each data set we are fitting to
            I3MapPMTKeyVectorDouble this_data; 
            I3MapPMTKeyVectorDouble this_pred; 
            I3MapPMTKeyVectorDouble this_times; 
            
            for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++){
                // loop over each PMT
                CCMPMTKey key = keys_to_fit.at(pmt_it);
                double llh = all_constructorsAD.at(data_it)->ComputeNLLH(key, params.at(0), Rt, params.at(1), tau_t, tau_rec, params.at(2), params.at(3), 
                                 time_offsets.at(data_it).at(key), const_offset, LPmu.at(key), LPsigma.at(key), LPscale.at(key), params.at(11 + pmt_it),
                                 {params.at(4), params.at(5)}, params.at(6), photons_per_mev, params.at(7 + data_it), n_sodium_events.at(data_it), light_profile_type, use_g4_yields).value();
                // now grab out data etc
                std::vector<double> this_pmt_data = all_constructorsAD.at(data_it)->GetDataVector(); 
                std::vector<double> this_pmt_pred = all_constructorsAD.at(data_it)->GetPredVector(); 
                std::vector<double> this_pmt_times = all_constructorsAD.at(data_it)->GetTimesVector(); 
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


