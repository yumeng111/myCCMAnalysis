
#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <nlopt.h>
#include <nlopt.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include <dataclasses/calibration/BaselineEstimate.h>

// struct to hold peak data when fitting
typedef struct {
    std::vector<double>::const_iterator v_start;
    std::vector<double>::const_iterator v_end;
    size_t t0;
} fit_data;


// struct to hold resulting best fit vals
struct BestFitParams{
    double c;
    double t0;
    double b1;
    double b2;
    double C;
    double T0;
    double B1;
    double B2;
};


class FitSummedPulses: public I3Module {
    std::string geometry_name_;
    bool geo_seen;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    static double GetPred(double & c, double & t0, double & b1, double & b2, double &t, double & C, double & T0, double & B1, double & B2);
    static std::vector<double> GetGradParams(double data, double pred, double & c, double & t0, double & b1, double & b2, double &t, double & C, double & T0, double & B1, double & B2);
    void FitPeak(std::vector<double> const & summed_pulses , std::pair<CCMPMTKey, BestFitParams> & best_fit_vals_per_pmt, int const & peak_position);

    static double PeakLossFunction(const std::vector<double> &x,
            std::vector<double> &grad,
            void * f_data);
    void AddTemplates(I3FramePtr frame);

    public:
    FitSummedPulses(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(FitSummedPulses);

FitSummedPulses::FitSummedPulses(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    }


void FitSummedPulses::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
}


void FitSummedPulses::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void FitSummedPulses::DAQ(I3FramePtr frame) {

    // let's read in our summed pulses
    I3Map<CCMPMTKey, std::vector<double>> const & SummedPulses = frame->Get<I3Map<CCMPMTKey, std::vector<double>> const> ("SummedPulses");
    I3Map<CCMPMTKey, std::vector<unsigned int>> const & counts = frame->Get<I3Map<CCMPMTKey, std::vector<unsigned int>> const> ("SummedPulsesCounts");
    I3Map<CCMPMTKey, int> const & peak_pos = frame->Get<I3Map<CCMPMTKey, int> const> ("SummedPulsesPeakPositions");

    // I3Map to store pmt key and best fit params
    I3Map<CCMPMTKey, BestFitParams> best_fit_vals;

    // loop over each pmt
    for(std::pair<CCMPMTKey const, std::vector<double>> const & it : SummedPulses){
        CCMPMTKey key = it.first;
        std::vector<double> summed_pulses = it.second;
        int peak_pos_per_key = peak_pos.at(key);

        std::pair<CCMPMTKey, BestFitParams> best_fit_vals_per_pmt;

        // now let's fit our summed pulse
        FitPeak(summed_pulses, best_fit_vals_per_pmt, peak_pos_per_key);
        //best_fit_vals.insert(std::make_pair<key, best_fit_vals_per_pmt>);
    }

    //frame->Put("8ParamSPETemplate", best_fit_vals);
    std::cout << "finished fitting SPEs!" << std::endl;
    PushFrame(frame);
}


double FitSummedPulses::GetPred(double & t, double & c, double & t0, double & b1, double & b2, double & C, double & T0, double & B1, double & B2){
    // we want to calculate the prediction value for a few values within this bin and return the average
    double total_pred;
    size_t n_samples_per_bin = 40;

    for(size_t idx = 0; idx<n_samples_per_bin; ++idx){
        double t_it = idx * 2.0 / double(n_samples_per_bin) + t;
        double pred = c / std::pow(std::exp(-(t_it- t0) / b1) + std::exp((t_it- t0) / b2) , 8) + C / std::pow(std::exp(-(t_it- T0) / B1) + std::exp((t_it- T0) / B2) , 8);
        total_pred += pred;
    }

    return total_pred/double(n_samples_per_bin);

}

std::vector<double> FitSummedPulses::GetGradParams(double data, double pred, double & c, double & t0, double & b1, double & b2, double &t, double & C, double & T0, double & B1, double & B2){

    std::vector<double> grad_vec (8);

    double total_grad0;
    double total_grad1;
    double total_grad2;
    double total_grad3;
    double total_grad4;
    double total_grad5;
    double total_grad6;
    double total_grad7;

    size_t n_samples_per_bin = 40;

    for(size_t idx = 0; idx<n_samples_per_bin; ++idx){
        double t_it = idx * 2.0 / double(n_samples_per_bin) + t;
        double prefactor = -2 * (data - pred) * pred;
        double prefactor2 = -8 * prefactor / (std::exp(-(t_it - t0)/b1) + std::exp( (t_it - t0)/b2));

        total_grad0 += prefactor / c;
        total_grad1 += prefactor2 * (std::exp(-(t_it - t0)/b1)/b1 - std::exp( (t_it - t0)/b2)/b2);
        total_grad2 += prefactor2 *  (t_it - t0)/std::pow(b1,2) * std::exp(-(t_it - t0)/b1);
        total_grad3 += prefactor2 * -(t_it - t0)/std::pow(b2,2) * std::exp( (t_it - t0)/b2);

        double Prefactor2 = -8 * prefactor / (std::exp(-(t_it - T0)/B1) + std::exp( (t_it - T0)/B2));

        total_grad4 += prefactor / C;
        total_grad5 += Prefactor2 * (std::exp(-(t_it - T0)/B1)/B1 - std::exp( (t_it - T0)/B2)/B2);
        total_grad6 += Prefactor2 *  (t_it - T0)/std::pow(B1,2) * std::exp(-(t_it - T0)/B1);
        total_grad7 += Prefactor2 * -(t_it - T0)/std::pow(B2,2) * std::exp( (t_it - T0)/B2);

    }

    double grad0 = total_grad0/double(n_samples_per_bin);
    double grad1 = total_grad1/double(n_samples_per_bin);
    double grad2 = total_grad2/double(n_samples_per_bin);
    double grad3 = total_grad3/double(n_samples_per_bin);
    double grad4 = total_grad4/double(n_samples_per_bin);
    double grad5 = total_grad5/double(n_samples_per_bin);
    double grad6 = total_grad6/double(n_samples_per_bin);
    double grad7 = total_grad7/double(n_samples_per_bin);

    grad_vec[0] = grad0;
    grad_vec[1] = grad1;
    grad_vec[2] = grad2;
    grad_vec[3] = grad3;
    grad_vec[4] = grad4;
    grad_vec[5] = grad5;
    grad_vec[6] = grad6;
    grad_vec[7] = grad7;

    return grad_vec;
}

double FitSummedPulses::PeakLossFunction(const std::vector<double> & x, std::vector<double> &grad, void * f_data) {

    // typecast to the our struct
    fit_data* d = (fit_data*) f_data;


    // initialize summation parameters
    double squared_residuals = 0;
    double pred, data;
    if (!grad.empty()) {
        for(size_t i = 0; i < grad.size(); ++i) grad[i] = 0;
    }

    // loop over each time bin, add to square residuals
    double c = x[0];
    double t0 = x[1];
    double b1 = x[2];
    double b2 = x[3];
    double C = x[4];
    double T0 = x[5];
    double B1 = x[6];
    double B2 = x[7];

    double t = d->t0;
    for(std::vector<double>::const_iterator it = d->v_start; it != d->v_end; ++it, t += 2) {
        pred = FitSummedPulses::GetPred(t, c, t0, b1, b2, C, T0, B1, B2);
        data = double(*it);

        squared_residuals += std::pow(data - pred, 2);
        // gradient calculated using SPE template function in I3DOMCalibration
        if(!grad.empty()) {
            std::vector<double> grad_vec = GetGradParams(data, pred, c, t0, b1, b2, t, C, T0, B1, B2);
            grad[0] += grad_vec[0];
            grad[1] += grad_vec[1];
            grad[2] += grad_vec[2];
            grad[3] += grad_vec[3];
            grad[4] += grad_vec[4];
            grad[5] += grad_vec[5];
            grad[6] += grad_vec[6];
            grad[7] += grad_vec[7];
        }
    }

    return squared_residuals;
}



void FitSummedPulses::FitPeak(std::vector<double> const & summed_pulses , std::pair<CCMPMTKey, BestFitParams> & best_fit_vals_per_pmt, int const & peak_position){
    size_t num_params = 8;
    //nlopt::opt opt(nlopt::LD_LBFGS, num_params);
    nlopt::opt opt(nlopt::LN_BOBYQA, num_params);

    size_t dont_fit_region = 1000;
    fit_data f_data;
    f_data.t0 = 0;
    f_data.v_start = summed_pulses.begin() + dont_fit_region;
    f_data.v_end = summed_pulses.end() - dont_fit_region;

    opt.set_min_objective(FitSummedPulses::PeakLossFunction, &f_data);

    // initial guesses
    std::vector<double> x(num_params);
    std::vector<double> lb(num_params);
    std::vector<double> ub(num_params);

    x[0] = 100;
    x[1] = peak_position * 2 - 10;
    x[2] = 0.5;
    x[3] = 20.0;
    x[4] = -0.5;
    x[5] = peak_position * 2 + 10;
    x[6] = 0.005;
    x[7] = 0.05;

    // guess some lower and upper bounds on each parameter
    lb[0] = 0.0; ub[0] = 10*x[0];
    lb[1] = peak_position * 2 - 40; ub[1] = peak_position * 2;
    lb[2] = 0.01*x[2]; ub[2] = 10*x[2];
    lb[3] = 0.01*x[3]; ub[3] = 10*x[3];
    lb[4] = -100.0; ub[4] = 0.0;
    lb[5] = peak_position * 2; ub[5] = peak_position * 2 + 100;
    lb[6] = 0.0; ub[6] = 100.0;
    lb[7] = 0.0; ub[7] = 100.0;

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    nlopt::result result;
    // now perform the minimization
    double minf;
    try {
        opt.set_xtol_rel(1e-12);
        opt.set_ftol_rel(1e-20);
        opt.set_maxtime(60);
        //opt.set_stopval(-HUGE_VAL);
        result = opt.optimize(x, minf);
    } catch(std::exception &e) {
    }

    //best_fit_vals_per_pmt.c = x[0];
    //best_fit_vals_per_pmt.t0 = x[1];
    //best_fit_vals_per_pmt.b1 = x[2];
    //best_fit_vals_per_pmt.b2 = x[3];
    //best_fit_vals_per_pmt.C = x[4];
    //best_fit_vals_per_pmt.T0 = x[5];
    //best_fit_vals_per_pmt.B1 = x[6];
    //best_fit_vals_per_pmt.B2 = x[7];

    std::cout << "nlopt result = " << result << std::endl;
    std::cout << "fit vals = " << x[0] << " , " << x[1] << " , " << x[2] << " , " << x[3] << " , " << x[4] << " , " << x[5] << " , " << x[6] << " , " << x[7] << std::endl;
}


void FitSummedPulses::Finish() {
}



