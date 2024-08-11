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

#include <fenv.h>

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

struct LikelihoodFunctor {

    static constexpr int DerivativeDimension = 8; // max number of dimensions :
                                                   // Rs, tau_s, tau_other
                                                   // time offset, norm
                                                   // LP mu, LP sigma, LP scale
    typedef double Underlying;
    typedef phys_tools::autodiff::FD<DerivativeDimension, Underlying> AD;

    std::vector<std::shared_ptr<CalculateNLLH<AD>>> llh_constructorAD;
    CCMPMTKey key_to_fit;
    std::vector<double> z_offset;
    size_t n_sodium_events;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified;
    std::vector<double> time_offsets;
    ZigZagPrior prior = ZigZagPrior(3.0, 3.0, false); // we want (tau_s - tau_TPB) > 3.0 with 600% scale
    bool fit_z_rayl = false;
    double nominal_rayl = 94.0;
    double Rt = 0.0;
    double tau_t = 743.0;
    double tau_rec = 0.0;
    double const_offset = 0.0;
    double photons_per_mev = 1.0;
    double pmt_eff = 1.0;
    double uv_abs_1 = 40.0;
    double uv_abs_2 = 0.0;
    std::vector<AD> all_uvs = {uv_abs_1, uv_abs_2};

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        T total_llh = 0;
        T Rs = x[0];
        T tau_s = x[1];
        T tau_other = x[2];
        T time_offset = x[3];
        T norm = x[4];
        T LPmu = x[5];
        T LPsigma = x[6];
        T LPscale = x[7];
        for (size_t data_it = 0; data_it < llh_constructorAD.size(); data_it ++){
            //T time_offset;
            //if (time_offsets.size() > 0){
            //    time_offset = time_offsets.at(data_it);
            //} else {
            //    time_offset = x[3];
            //}
            total_llh += llh_constructorAD.at(data_it)->ComputeNLLH(key_to_fit, Rs, Rt, tau_s, tau_t, tau_rec, tau_other, norm,
                            time_offset, const_offset, LPmu, LPsigma, LPscale, pmt_eff, all_uvs, nominal_rayl, photons_per_mev,
                            z_offset.at(data_it), n_sodium_events, light_profile_type, fit_z_rayl);

                // and now add our prior!
            total_llh += prior(tau_s - tau_other);
            //else {
            //    total_llh += std::numeric_limits<T>::infinity();
            //}

        }
        //std::cout << "Rs = " << Rs.value() << ", tau_s = " << tau_s.value() << ", tau_other = " << tau_other.value()
        //    << ", time offset = " << time_offset.value() << ", norm = " << norm.value() << ", LP mu = " << LPmu.value() << ", LPsigma = " << LPsigma.value()
        //    << ", LPscale = " << LPscale.value() << ", and total llh = " << total_llh << std::endl;

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
        //return(-func.template evaluateLikelihood<double>(x)); // note negation!
        return 0.0;
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

void cppMinimizer::OnePMTOneDataSetMinimization(std::vector<CCMPMTKey> all_keys_to_fit, std::vector<std::string> data_file_names, std::vector<double> z_offsets,
                                                double norm_seed, size_t n_sodium_events) {
    //feenableexcept(FE_INVALID);
    // let's initialize our constructor

    // make our AD object
    LikelihoodType likelihood;
    typedef double Underlying;
    typedef phys_tools::autodiff::FD<likelihood.DerivativeDimension, Underlying> AD;

    // now grab data frame(s) to make our constructors
    std::vector<std::shared_ptr<CalculateNLLH<AD>>> all_constructorsAD;
    for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
        std::string this_data_fname = data_file_names.at(data_it);
        dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
        I3FramePtr this_data_frame = this_data_file.pop_frame();

        // and now make constructor -- AD type
        std::shared_ptr<CalculateNLLH<AD>> this_llh_constructorAD = std::make_shared<CalculateNLLH<AD>>();
        this_llh_constructorAD->SetKeys(static_cast<I3VectorCCMPMTKey>(all_keys_to_fit));
        this_llh_constructorAD->SetData(this_data_frame);

        // now save!
        all_constructorsAD.push_back(this_llh_constructorAD);
    }

    for (size_t k = 0; k < all_keys_to_fit.size(); k++){
        CCMPMTKey this_key = all_keys_to_fit.at(k);
        std::cout << "fitting " << this_key << std::endl;
        // let's try scanning over t offset
        double t_offset_start = 12.0;
        double t_offset_end = 35.0;
        double t_offset_range = t_offset_end - t_offset_start;
        size_t n_t_offsets = 60;

        std::vector<double> func_val;
        std::vector<std::vector<double>> params;
        std::vector<int> nevals;

        // now set up our likelihood object
        LikelihoodType likelihood;

        likelihood.llh_constructorAD = all_constructorsAD;
        likelihood.n_sodium_events = n_sodium_events;
        likelihood.z_offset = z_offsets;
        likelihood.key_to_fit = this_key;

        for (size_t toff_it = 0; toff_it < (n_t_offsets + 1); toff_it++){
            double this_t_offset = t_offset_start + ((double)toff_it / (double)n_t_offsets) * t_offset_range;

            phys_tools::lbfgsb::LBFGSB_Driver minimizer;

            minimizer.addParameter(0.35, 1e-1, 0.1, 0.38); // Rs
            minimizer.addParameter(6.0, 1e-1, 1.0, 16.0); // tau_s
            minimizer.addParameter(1.0, 1e-1, 0.0001, 9.0); // tau_other
            minimizer.addParameter(this_t_offset); // time offset
            minimizer.addParameter(1.0, 1e-1, 0.001, 2.5e6); //norm
            minimizer.addParameter(50.0, 1e-1, 30.0, 70.0); // late pulse mu
            minimizer.addParameter(8.0, 1e-1, 3.0, 20.0); // late pulse sigma
            minimizer.addParameter(0.01, 1e-1, 1e-5, 0.1); // late pulse scale
            minimizer.setHistorySize(20);

            // fix parameter idx of guys we are not minimizing
            //minimizer.fixParameter(0); // Rs
            //minimizer.fixParameter(1); // tau_s
            //minimizer.fixParameter(2); // tau_other
            minimizer.fixParameter(3); // time offset

            bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));

            if(succeeded) {
                func_val.push_back(minimizer.minimumValue());
                params.push_back(minimizer.minimumPosition());
                nevals.push_back(minimizer.numberOfEvaluations());
                std::cout << "for time offset = " << this_t_offset << ", function minimum at " << minimizer.minimumValue() << " and minimization message = " << minimizer.errorMessage() << std::endl;
            } else{
                std::cout << "oops! for time offset = " << this_t_offset << ", minimizer says = " << minimizer.errorMessage() << std::endl;
            }


        }

        // now let's find the smallest func_val!
        size_t smallest_idx = std::distance(func_val.begin(), std::min_element(func_val.begin(), func_val.end()));
        std::cout << "Best evaluation!!! Function value at minimum: " << func_val.at(smallest_idx) << " after " << nevals.at(smallest_idx) << " function evaluations" << std::endl;
        std::cout << "Parameters: " << std::endl;
        std::vector<double> data_to_return;
        for(size_t i=0; i<LikelihoodFunctor::DerivativeDimension; ++i) {
            data_to_return.push_back(params.at(smallest_idx).at(i));
            std::cout << params.at(smallest_idx).at(i) << std::endl;
        }
        data_to_return.push_back(func_val.at(smallest_idx));

        // one last thing -- take this best fit point and grab our data, pred, and times
        AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified;
        double nominal_rayl = 94.0;
        double Rt = 0.0;
        double tau_t = 743.0;
        double tau_rec = 0.0;
        double const_offset = 0.0;
        double photons_per_mev = 1.0;
        double pmt_eff = 1.0;
        double uv_abs_1 = 40.0;
        double uv_abs_2 = 0.0;
        bool fit_z_rayl = false;
        std::vector<AD> all_uvs = {uv_abs_1, uv_abs_2};
        std::vector<CCMPMTKey> keys_to_fit = {this_key};
        double Rs = data_to_return.at(0);
        double tau_s = data_to_return.at(1);
        double tau_other = data_to_return.at(2);
        double time_offset = data_to_return.at(3);
        double norm = data_to_return.at(4);
        double LPmu = data_to_return.at(5);
        double LPsigma = data_to_return.at(6);
        double LPscale = data_to_return.at(7);
        for (size_t data_it = 0; data_it < all_constructorsAD.size(); data_it ++){
            // loop over each data set we are fitting to
            I3MapPMTKeyVectorDouble this_data;
            I3MapPMTKeyVectorDouble this_pred;
            I3MapPMTKeyVectorDouble this_times;

            for (size_t pmt_it = 0; pmt_it < keys_to_fit.size(); pmt_it ++){
                // loop over each PMT
                CCMPMTKey key = keys_to_fit.at(pmt_it);

                double llh = all_constructorsAD.at(data_it)->ComputeNLLH(key, Rs, Rt, tau_s, tau_t, tau_rec, tau_other, norm,
                                time_offset, const_offset, LPmu, LPsigma, LPscale, pmt_eff, all_uvs, nominal_rayl, photons_per_mev,
                                z_offsets.at(data_it), n_sodium_events, light_profile_type, fit_z_rayl).value();

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

        all_data_to_return.push_back(data_to_return);
    }
}

std::vector<double> cppMinimizer::OnePMTMultipleDataSetMinimization(CCMPMTKey this_key, std::vector<std::string> data_file_names, std::vector<double> z_offsets,
                            std::vector<double> time_offsets, std::vector<double> norm_seeds, size_t n_sodium_events, bool use_g4_yields) {

    //std::vector<std::string> paramter_names = {"Rs", "Rt", "tau_s", "tau_t", "tau_rec", "tau_TPB", "norm1", "norm2", "norm3", "norm4",
    //                                           "time offset", "const offset", "late pulse mu", "late pulse sigma", "late pulse scale", "uv absorption"};

    //size_t n_data_sets = data_file_names.size();

    //// grab our geomtry frame
    //std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
    //dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
    //I3FramePtr geo_frame = geometry_file.pop_frame();
    //
    //// now set up our likelihood object
    //LikelihoodType likelihood;
    //
    //// make our AD object
    //typedef double Underlying;
    //typedef phys_tools::autodiff::FD<likelihood.DerivativeDimension, Underlying> AD;
    //
    //// now grab data frame(s) to make our constructors
    //std::vector<std::shared_ptr<CalculateNLLH<AD>>> all_constructorsAD;
    //std::vector<std::shared_ptr<CalculateNLLH<double>>> all_constructorsDouble;
    //for (size_t data_it = 0; data_it < data_file_names.size(); data_it++){
    //    std::string this_data_fname = data_file_names.at(data_it);
    //    dataio::I3File this_data_file(this_data_fname, dataio::I3File::Mode::read);
    //    I3FramePtr this_data_frame = this_data_file.pop_frame();

    //    // and now make constructor -- AD type
    //    std::shared_ptr<CalculateNLLH<AD>> this_llh_constructorAD = std::make_shared<CalculateNLLH<AD>>();
    //    this_llh_constructorAD->SetKeys(I3VectorCCMPMTKey({this_key}));
    //    this_llh_constructorAD->SetData(this_data_frame);
    //    this_llh_constructorAD->SetGeo(geo_frame);
    //
    //    // now save!
    //    all_constructorsAD.push_back(this_llh_constructorAD);
    //
    //    // and now make constructor -- double type
    //    std::shared_ptr<CalculateNLLH<double>> this_llh_constructorDouble = std::make_shared<CalculateNLLH<double>>();
    //    this_llh_constructorDouble->SetKeys(I3VectorCCMPMTKey({this_key}));
    //    this_llh_constructorDouble->SetData(this_data_frame);
    //    this_llh_constructorDouble->SetGeo(geo_frame);
    //
    //    // now save!
    //    all_constructorsDouble.push_back(this_llh_constructorDouble);
    //}
    //
    //likelihood.llh_constructorAD = all_constructorsAD;
    //likelihood.llh_constructorDouble = all_constructorsDouble;
    //likelihood.n_sodium_events = n_sodium_events;
    //likelihood.z_offset = z_offsets;
    //likelihood.key_to_fit = this_key;
    //likelihood.time_offsets = time_offsets;
    //likelihood.use_g4_yields = use_g4_yields;
    //
    //phys_tools::lbfgsb::LBFGSB_Driver minimizer;

    //// free paramters are : Rs, tau_s, tau_TPB, normalization, late pulse mu, late pulse sigma, and late pulse scale
    //// fixed paramters are : Rt, tau_t, tau_rec, time offset, and const_offset
    //double seeds[7] = {0.34, 8.2, 2.0, 50.0, 55.0, 8.0, 3e-2};
    //double grad_scales[7] = {1e-3, 1e-3, 1e-3, 1e-6, 1e-3, 1e-3, 1e-3};
    //double mins[7] = {0.2, 2.0, 1.0, 1.0, 35.0, 2.0, 0.0};
    //double maxs[7] = {0.5, 16.0, 9.0, 1e10, 68.0, 20.0, 1e-1};
    //
    //minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]); // Rs
    //minimizer.addParameter(0.0);                                        // Rt
    //minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]); // tau_s
    //minimizer.addParameter(743.0);                                      // tau_t
    //minimizer.addParameter(0.0);                                        // tau_rec
    //minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]); // tau_TPB
    //// this is where we add our normalization parameter...need to add one for each data set
    //for (size_t n = 0; n < n_data_sets; n++){
    //    minimizer.addParameter(norm_seeds.at(n), grad_scales[3], norm_seeds.at(n) / 1e2, norm_seeds.at(n) * 1e5); // norm
    //}
    //minimizer.addParameter(0.0);                                        // time offset
    //minimizer.addParameter(0.0);                                        // const offset
    //minimizer.addParameter(seeds[4], grad_scales[4], mins[4], maxs[4]); // late pulse mu
    //minimizer.addParameter(seeds[5], grad_scales[5], mins[5], maxs[5]); // late pulse sigma
    //minimizer.addParameter(seeds[6], grad_scales[6], mins[6], maxs[6]); // late pulse scale
    //minimizer.addParameter(55.0, 1e-3, 40.0, 60.0);                     // uv absorption
    //minimizer.setHistorySize(20);

    //// fix parameter idx of guys we are not minimizing
    //minimizer.fixParameter(1); // Rt
    //minimizer.fixParameter(3); // tau_t
    //minimizer.fixParameter(4); // tau_rec
    //minimizer.fixParameter(10); // time offset
    //minimizer.fixParameter(11); // const offset
    //minimizer.fixParameter(15); // uv absorption

    //bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));

    //if (succeeded){
    //    std::cout << "joint fit converged!" << std::endl;
    //    std::cout << "minimization finished with " << minimizer.errorMessage() << " error message" << std::endl;
    //}

    //double value = minimizer.minimumValue();
    //std::vector<double> params = minimizer.minimumPosition();
    //
    //// now let's find the smallest func_val!
    //std::cout << "Function value at minimum: " << value << " after " << minimizer.numberOfEvaluations() << " function evaluations" << std::endl;
    //std::cout << "Parameters: " << std::endl;
    //std::vector<double> data_to_return;
    //data_to_return.push_back(value);
    //data_to_return.push_back((double) minimizer.numberOfEvaluations());
    //for(size_t i=0; i<LikelihoodFunctor::DerivativeDimension; ++i) {
    //    data_to_return.push_back(params.at(i));
    //    if (i == 1 or i == 3 or i == 4 or i == 10 or i == 11 or i == 15){
    //        continue;
    //    }
    //    std::cout << paramter_names.at(i) << " = " << params.at(i) << std::endl;
    //}
    std::vector<double> data_to_return;
    return data_to_return;
}


