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

    static constexpr int DerivativeDimension = 3;

    // This returns the LLH
    template<typename T>
    T evaluateLikelihood(std::vector<T> x) const {
        // set up our alygen
        AnalyticLightYieldGenerator alygen = AnalyticLightYieldGenerator();
        alygen.uv_absorption = 55.0;
        alygen.n_sodium_events = 20;
        alygen.z_offset = 0.0; 
        alygen.light_profile_type = AnalyticLightYieldGenerator::LArLightProfileType::Simplified;
        alygen.tau_t = 743.0;
        alygen.Rs = x[0];
        alygen.Rt = 1.0 - x[0];
        alygen.tau_s = x[1];
        alygen.normalization = 1e5; 
        alygen.tau_TPB = x[2];
        alygen.time_offset = 25.4; 
        
        return llh_constructor->GetNLLH(CCMPMTKey(0,2,0), alygen, 60.0, 10.0, 1e3); 
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
        for(size_t i=0; i<12; i++){
            if (i == 0 or i == 2 or i == 5){
                params[i] = GradType(x[i],i);
            }
        }
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
    llh_constructor->SetKeys({CCMPMTKey(0,2,0)});
    llh_constructor->SetData(data_frame);
    llh_constructor->SetGeo(geo_frame);

    LikelihoodType likelihood;

    phys_tools::lbfgsb::LBFGSB_Driver minimizer;

    double seeds[3]       = {0.34,    8.2,   2.5    };
    double grad_scales[3] = {1,    1,    1   };
    double mins[3]        = {1e-3, 5.0, 1.0};
    double maxs[3]        = {1.0,  16.0,  9.0};

    minimizer.addParameter(seeds[0], grad_scales[0], mins[0], maxs[0]);
    minimizer.addParameter(seeds[1], grad_scales[1], mins[1], maxs[1]);
    minimizer.addParameter(seeds[2], grad_scales[2], mins[2], maxs[2]);
    minimizer.setHistorySize(20);

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
