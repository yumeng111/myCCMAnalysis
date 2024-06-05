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
#include "dataclasses/physics/AnalyticLightYieldGenerator.h"

struct LikelihoodFunctor {

    static constexpr int DerivativeDimension = 12;
    std::shared_ptr<CalculateNLLH> llh_constructor;
    CCMPMTKey key_to_fit;
    double uv_absorption = 55.0;
    double z_offset = 0.0;
    size_t n_sodium_events = 20;
    AnalyticLightYieldGenerator::LArLightProfileType light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified; 

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        return llh_constructor->ComputeNLLH<T>(key_to_fit, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], uv_absorption, z_offset, n_sodium_events, light_profile_type);
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

int main(int argc, char ** argv) {

    // let's initialize our constructor
    // grab our geomtry frame
    std::string geometry_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/geo_run012490.i3.zst";
    dataio::I3File geometry_file(geometry_fname, dataio::I3File::Mode::read);
    I3FramePtr geo_frame = geometry_file.pop_frame();
    
    // now grab a data frame
    std::string data_fname = "/Users/darcybrewuser/workspaces/CCM/notebooks/sodium_center_accumulated_events.i3.zst";
    dataio::I3File data_file(data_fname, dataio::I3File::Mode::read);
    I3FramePtr data_frame = data_file.pop_frame();

    std::shared_ptr<CalculateNLLH> llh_constructor = std::make_shared<CalculateNLLH>();
    // doing this for one PMT at the moment...
    CCMPMTKey this_key = CCMPMTKey(0,2,0);
    llh_constructor->SetKeys({this_key});
    llh_constructor->SetData(data_frame);
    llh_constructor->SetGeo(geo_frame);

    LikelihoodType likelihood;
    likelihood.llh_constructor = llh_constructor;
    likelihood.key_to_fit = this_key;

    phys_tools::lbfgsb::LBFGSB_Driver minimizer;

    // free paramters are : Rs, tau_s, tau_TPB, normalization, time offset, late pulse mu, late pulse sigma, and late pulse scale
    // fixed paramters are : Rt, tau_t, tau_rec, and const_offset
    double seeds[8] = {0.34, 8.2, 3.0, 1e5, 25.4, 60.0, 10.0, 1e3};
    double grad_scales[8] = {1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3};
    double mins[8] = {1e-3, 5.0, 1.0, 1.0, 20.0, 40.0, 2.0, 0.0};
    double maxs[8] = {1.0, 16.0, 9.0, 1e8, 30.0, 68.0, 30.0, 1e7};

    minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]); // Rs
    minimizer.addParameter(0.0);                                         // Rt
    minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]); // tau_s
    minimizer.addParameter(743.0);                                       // tau_t
    minimizer.addParameter(0.0);                                         // tau_rec
    minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]); // tau_TPB
    minimizer.addParameter(seeds[3], grad_scales[2], mins[2], maxs[2]); // norm
    minimizer.addParameter(seeds[4], grad_scales[4], mins[4], maxs[4]); // time offset
    minimizer.addParameter(0.0);                                         // const offset
    minimizer.addParameter(seeds[5], grad_scales[5], mins[5], maxs[5]); // late pulse mu
    minimizer.addParameter(seeds[6], grad_scales[6], mins[6], maxs[6]); // late pulse sigma
    minimizer.addParameter(seeds[7], grad_scales[7], mins[7], maxs[7]); // late pulse scale
    minimizer.setHistorySize(20);

    // fix parameter idx of guys we are not minimizing
    minimizer.fixParameter(1); // Rt
    minimizer.fixParameter(3); // tau_t
    minimizer.fixParameter(4); // tau_rec
    minimizer.fixParameter(8); // const offset

    bool succeeded = minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));

    if(succeeded) {
        std::cout << "Success!!" << std::endl;
    } else {
        std::cout << "Failure..." << std::endl;
    }

    double value = minimizer.minimumValue();
    std::vector<double> params = minimizer.minimumPosition();
    int nEval = minimizer.numberOfEvaluations();
    int nGrad = minimizer.numberOfEvaluations();
    std::string message = minimizer.errorMessage();


    std::cout << "Function value at minimum: " << value << std::endl;
    std::cout << "Parameters: (";
    for(size_t i=0; i<LikelihoodFunctor::DerivativeDimension; ++i) {
        if(i > 0) {
            std::cout << ", ";
        }
        std::cout << params[i];
    }
    std::cout << ")" << std::endl;
    std::cout << "Function evaluations: " << nEval << std::endl;
    std::cout << "Gradient evaluations: " << nGrad << std::endl;
    std::cout << "Message: \"" << message << "\"" << std::endl;

}
