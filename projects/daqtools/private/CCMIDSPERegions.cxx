// this module identifies regions where we want to look for SPEs 
// then fits those regions using Nick's code

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
    std::vector<short unsigned int>::const_iterator v_start;
    std::vector<short unsigned int>::const_iterator v_end;
    int64_t t0;
    const double* baseline;
    int64_t peak_tmax;
    int64_t peak_amplitude;
} fit_data;

// struct to hold gradient vector values from 4 param fit 
struct FourParamGrad {
    double grad0;
    double grad1;
    double grad2;
    double grad3;
} ;


class CCMIDSPERegions: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    std::tuple<size_t, double, double> CheckForPulse(WaveformSmoother & smoother, double const & mode, size_t start_idx);
    void ProcessWaveform(WaveformSmoother & smoother,  CCMWaveformUInt16 const & waveform, double const & mode, I3Vector<SPETemplate> & template_per_channel, I3Vector<int64_t> & time_per_channel, I3Vector<int64_t> & length_per_channel, I3Vector<double> & amp_per_channel, I3Vector<short unsigned int> & SPE_wf_to_fit);
    void FindRegions(WaveformSmoother & smoother, std::vector<short unsigned int> const & samples, double const & mode, I3Vector<short unsigned int> & SPE_wf_to_fit);
    void GetPeakInfo(I3Vector<short unsigned int> & SPE_wf_to_fit, I3Vector<int> & window_amplitude, I3Vector<int64_t>& window_length, I3Vector<int64_t> &window_time);
    static double GetPred(double & c, double & t0, double & b1, double & b2, double &t);
    static std::pair<double, double> GetGradTwoParams(double & amp, double & data, double & pred, double & c, double & t0, double & b1, double & b2, double &t);
    static FourParamGrad GetGradFourParams(double data, double pred, double & c, double & t0, double & b1, double & b2, double &t);
    SPETemplate FitPeak(double peak_amplitude,
            int64_t length,
            int64_t time,
            std::vector<short unsigned int> const & waveform,
            double baseline,
            int num_params, 
            SPETemplate two_param_fit_vals = SPETemplate());

    int64_t fit_length = 100;
    int64_t pre_window;
    int min_amplitude;
    int max_amplitude;
    size_t n_daq_frames = 0;
    // int num_params;
    //static const int num_params = 2;

    double smoother_tau = 10.0 ;
    double smoother_delta_t = 2.0;


    double AmplitudeConstraint(const std::vector<double> &x,
            std::vector<double> &grad,
            void * data);

    static double PeakLossFunction(const std::vector<double> &x,
            std::vector<double> &grad,
            void * f_data);
    void AddTemplates(I3FramePtr frame);

    public:
    CCMIDSPERegions(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(CCMIDSPERegions);

CCMIDSPERegions::CCMIDSPERegions(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    }


void CCMIDSPERegions::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
}


void CCMIDSPERegions::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void CCMIDSPERegions::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    size_t size = waveforms->size();

    // a vector storing the SPE template and time for each channel
    boost::shared_ptr<I3Vector<I3Vector<SPETemplate>>> channel_templates(new I3Vector<I3Vector<SPETemplate>>(size));
    boost::shared_ptr<I3Vector<I3Vector<int64_t>>> channel_times(new I3Vector<I3Vector<int64_t>>(size));
    boost::shared_ptr<I3Vector<I3Vector<int64_t>>> channel_length(new I3Vector<I3Vector<int64_t>>(size));
    boost::shared_ptr<I3Vector<I3Vector<double>>> peak_amplitude(new I3Vector<I3Vector<double>>(size));
    boost::shared_ptr<I3Vector<I3Vector<short unsigned int>>> all_SPE_regions(new I3Vector<I3Vector<short unsigned int>>(size));

    // loop over each pmt 
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second;
        double mode = value.baseline * -1; //baseline mode is positive!!!
        uint32_t channel = pmt_channel_map_[key];

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);
        I3Vector<SPETemplate> template_per_channel; 
        I3Vector<int64_t> time_per_channel;
        I3Vector<int64_t> length_per_channel;
        I3Vector<double> amp_per_channel;
        I3Vector<short unsigned int> SPE_wf_to_fit(waveform.GetWaveform().size(), 0);
        WaveformSmoother smoother(waveform.GetWaveform().cbegin(), waveform.GetWaveform().cend(), smoother_delta_t, smoother_tau);

        ProcessWaveform(smoother, waveform, mode, template_per_channel, time_per_channel, length_per_channel, amp_per_channel, SPE_wf_to_fit);
        channel_templates->operator[](channel) = template_per_channel;
        channel_times->operator[](channel) = time_per_channel;
        channel_length->operator[](channel) = length_per_channel;
        all_SPE_regions->operator[](channel) = SPE_wf_to_fit;

    }

    frame->Put("SPETemplates", channel_templates);
    frame->Put("SPETemplatesTimes", channel_times);
    frame->Put("SPETemplatesLength", channel_length);
    frame->Put("SPEDataAmplitudes", peak_amplitude);
   // frame->Put("SPEPeakRegions", all_SPE_regions);
    std::cout << "finished fitting SPEs!" << std::endl;
    PushFrame(frame);
}

std::tuple<size_t, double, double> CCMIDSPERegions::CheckForPulse(WaveformSmoother & smoother, double const & mode, size_t start_idx){
    // 0 before start of pulse checking [checking for positive derivative]
    // 1 derivative is positive (in rising edge of pulse) [checking for negative derivative]
    // 2 derivative is negative (in falling edge of pulse) [checking for positive derivative]
    // 3 derivative is positive (recovering from droop) [done]

    int state = 1;
    double max_abs_derivative = std::abs(smoother.Derivative());
    double integral = 0;
    bool found_pulse = false;
    size_t final_idx = start_idx;
    size_t max_samples = 300;
    size_t N = std::min(smoother.Size(), start_idx + max_samples);
    double min_val_threshold = 20;
    double max_val_threshold = 40;
    double val_at_peak;
    bool at_turning_point = false;
    size_t turning_point_idx;
    double integral_threshold = 10.0;
    size_t min_length = 3;
    double inv_mode = -1*mode; //invert mode since smoother.Value() is negative
    double max_value = smoother.Value() - inv_mode;
    double max_unsmoothed_value = mode - smoother.RawValue();

    for(size_t i=start_idx; i<N; ++i) {
        max_value = std::max(max_value, smoother.Value() - inv_mode);
        max_unsmoothed_value = std::max(max_unsmoothed_value, mode - smoother.RawValue());
        max_abs_derivative = std::max(max_abs_derivative, std::abs(smoother.Derivative()));
        integral += (smoother.Value() - inv_mode);
        if(state % 2) {
            if(smoother.Derivative() < 0)
                state += 1;
        } else
            if(smoother.Derivative() > 0) {
                state += 1;
            }
        if(state >= 3) {
            found_pulse = true;
            final_idx = i;
            break;
        }
        smoother.Next();
    }
    smoother.Reset(start_idx);

    if(found_pulse
            // and (max_unsmoothed_value >= min_val_threshold)
            // and (max_unsmoothed_value <= max_val_threshold)
            // and (integral >= integral_threshold)
            // and (final_idx - start_idx >= min_length)
      ) {
        return std::make_tuple(final_idx, max_unsmoothed_value, integral);
    } else {
        return std::make_tuple(start_idx, max_unsmoothed_value, integral);
    }
}


void CCMIDSPERegions::FindRegions(WaveformSmoother & smoother, std::vector<short unsigned int> const & samples, double const & mode, I3Vector<short unsigned int> & SPE_wf_to_fit){
    // let's loop over derivs to find regions that might have SPEs
    // also good time to implement some cuts
    // cut on ADC counts above the baseline

    size_t N = smoother.Size(); 
    double deriv_threshold = 0.3;
    size_t pulse_first_index = 0;
    size_t pulse_last_index = 0;
    std::vector<size_t> pulse_tracker;
    double max_unsmoothed_value = 0;
    double integral = 0;
    double min_val_threshold = 20;
    double max_val_threshold = 40;
    double integral_threshold = 10.0;

    smoother.Reset();
    for(size_t i = 0; i < N; ++i){
        // looping over waveform 
        // need to look for peaks based on deriv > 0 then deriv < 0

        if(smoother.Derivative() > deriv_threshold){
            // derive is positive
            // let's call our check for pulse function
            std::tie(pulse_last_index, max_unsmoothed_value, integral) = CheckForPulse(smoother, mode, i);

            if(pulse_last_index > i) {
                pulse_first_index = size_t(std::max(ptrdiff_t(0), ptrdiff_t(i)));

                if(std::find(pulse_tracker.begin(), pulse_tracker.end(), pulse_last_index) != pulse_tracker.end()){
                    // we've already seen this pulse;
                }

                else{
                    // first time seeing this pulse!
                    pulse_tracker.push_back(pulse_last_index);
                    // let's add a cut here that continues if 100 bins before the start of our pulse and 100 bins after the end of our pulse are low charged
                    size_t pre_pulse_window;
                    size_t post_pulse_window;
                    size_t nearby_charge_window = 150;

                    // first checking the pre pulse window for activity
                    if (pulse_first_index < nearby_charge_window){
                        pre_pulse_window = pulse_first_index;
                    }
                    else{
                        pre_pulse_window = nearby_charge_window;
                    }
                    double max_val_pre_pulse = samples[pulse_first_index] - mode;
                    for (size_t pre_pulse_it = pulse_first_index - pre_pulse_window; pre_pulse_it < pulse_first_index; ++pre_pulse_it){
                        max_val_pre_pulse = std::min(max_val_pre_pulse, samples[pre_pulse_it] - mode);
                    }

                    // now checking the post pulse window for activity
                    if (pulse_last_index > (samples.size() - nearby_charge_window)){
                        post_pulse_window = samples.size() - pulse_last_index;
                    }
                    else{
                        post_pulse_window = nearby_charge_window;
                    }
                    double max_val_post_pulse = samples[pulse_last_index] - mode;
                    for (size_t post_pulse_it = pulse_last_index; post_pulse_it < pulse_last_index + post_pulse_window; ++post_pulse_it){
                        max_val_post_pulse = std::min(max_val_post_pulse, samples[post_pulse_it] - mode);
                    }
                    double surrounding_wf_threshold = 7; // counts above or below the baseline that we consider noise

                    if (max_val_pre_pulse < surrounding_wf_threshold and max_val_pre_pulse > -1*surrounding_wf_threshold
                        and max_val_post_pulse < surrounding_wf_threshold and max_val_post_pulse > -1*surrounding_wf_threshold) {
                        if((max_unsmoothed_value >= min_val_threshold) and (max_unsmoothed_value <= max_val_threshold) and (integral >= integral_threshold)){
                            // made it past our cuts on adc counts!
                            // let's save waveform to SPE_wf_to_fit
                            // std::cout << "found a good peak! starting bin = " << pulse_first_index << " , ending bin = " << pulse_last_index << std::endl;
                            for(size_t pulse_it = pulse_first_index; pulse_it <= pulse_last_index; ++pulse_it){
                                //std::cout << "wf in this peak = " << samples[pulse_it] << std::endl;
                                SPE_wf_to_fit[pulse_it] = samples[pulse_it];
                            }
                        }
                    }
                }

            }

        }
        smoother.Next();
    }
}


void CCMIDSPERegions::GetPeakInfo(I3Vector<short unsigned int>  & SPE_wf_to_fit, I3Vector<int>& window_amplitude, I3Vector<int64_t> &window_length, I3Vector<int64_t>& window_time){
    // need to fill in some vectors with peak info
    // looping over SPE_wf_to_fit that contains the wf in regions of interest

    bool in_peak_region = false;
    bool first_bin_of_peak = false;
    std::vector<short unsigned int> each_peak_wf;
    int64_t start_time;
    int64_t end_time;

    for(size_t i = 0; i < SPE_wf_to_fit.size(); ++i){

        if(SPE_wf_to_fit[i] != 0 and first_bin_of_peak == false){
            first_bin_of_peak = true;
            start_time = i*2; //start time in nsec
        }

        if(SPE_wf_to_fit[i] != 0){
            in_peak_region = true;
            each_peak_wf.push_back(SPE_wf_to_fit[i]);
        }

        if(in_peak_region == true and SPE_wf_to_fit[i] == 0){
            // at the end of peak region
            end_time = i*2; //end time in nsec

            // now let's loop over each_peak_wf and fill in some values
            int amp = each_peak_wf[0];
            int64_t length = each_peak_wf.size();
            //int64_t time = end_time - start_time; 
            int64_t time = start_time/2; 

            for(size_t j = 0; j < each_peak_wf.size(); ++j){
                if(each_peak_wf[j] < amp){
                    amp = each_peak_wf[j];
                } 
            }
            window_amplitude.push_back(amp);
            window_length.push_back(length);
            window_time.push_back(time);
            // clear vector that saves peak and reset variables 
            each_peak_wf.clear();
            in_peak_region = false;
            first_bin_of_peak = false;
        }
    }

}

double CCMIDSPERegions::GetPred(double & c, double & t0, double & b1, double & b2, double &t){
    // we want to calculate the prediction value for a few values within this bin and return the average
    double total_pred;
    size_t n_samples_per_bin = 40;

    for(size_t idx = 0; idx<n_samples_per_bin; ++idx){
        double t_it = idx * 2.0 / double(n_samples_per_bin) + t;
        double first_exp = std::exp(-(t_it- t0) / b1);
        double second_exp = std::exp((t_it - t0) / b2);
        double denom = std::pow(first_exp + second_exp , 8);
        double pred = c / denom;
        total_pred += pred;
    }

    return total_pred/double(n_samples_per_bin);

}

std::pair<double, double> CCMIDSPERegions::GetGradTwoParams(double & amp, double & data, double & pred, double & c, double & t0, double & b1, double & b2, double &t){

    double total_grad0;
    double total_grad1;
    size_t n_samples_per_bin = 40;

    for(size_t idx = 0; idx<n_samples_per_bin; ++idx){
        double t_it = idx * 2.0 / double(n_samples_per_bin) + t;
        double dwdb1 = - (8 * c * std::exp(-(t_it - t0)/b1) * (t_it - t0)) * (std::pow(b1,-2) * std::pow(std::exp(-(t_it - t0)/b1) + std::exp((t_it - t0)/b2), -9));
        double dwdc = (std::pow(std::exp(-(t_it - t0)/b1) + std::exp((t_it - t0)/b2), -8));
        double dcdb1 = ((8 * amp)/(b1 * std::pow((b1+b2), 2))) * std::pow( std::pow(b1/b2 , b2/(b1+b2)) + std::pow(b2/b1 , b1/(b1+b2)) , 7) * 
                        (std::pow(b1/b2 , b2/(b1+b2)) * b2 * (b1 + b2 - b1 * std::log(b1/b2)) - b1 * std::pow(b2/b1 , b1/(b1+b2))*(b1 + b2 - b2*std::log(b2/b1))); 
        double dwdt0 = -8 * c * ((std::exp(-(t_it - t0)/b1))/b1 - (std::exp((t_it - t0)/b2))/b2) * std::pow(std::exp(-(t_it - t0)/b1) + std::exp((t_it - t0)/b2) , -9);
        double dt0db1 = ((b2)/(b1+b2)) + ((b1 * b2 * std::log(b2/b1) )/(std::pow((b1 + b2), 2))) - ((b2 * std::log(b2/b1))/(b1+b2));

        double dwdb2 = (8 * c * std::exp((t_it - t0)/b2) * (t_it - t0)) * (std::pow(b2, -2) * std::pow(std::exp(-(t_it - t0)/b1) + std::exp((t_it - t0)/b2), -9));
        double dcdb2 = ((8 * amp)/(b2 * std::pow((b1+b2), 2))) * std::pow( std::pow(b1/b2 , b2/(b1+b2)) + std::pow(b2/b1 , b1/(b1+b2)) , 7) * 
                        (-std::pow(b1/b2 , b2/(b1+b2)) * b2 * (b1 + b2 - b1 * std::log(b1/b2)) + b1 * std::pow(b2/b1 , b1/(b1+b2))*(b1 + b2 - b2*std::log(b2/b1))); 
        double dt0db2 = ((-b1)/(b1+b2)) + ((b1 * b2 * std::log(b2/b1))/(std::pow((b1 + b2),2))) - ((b1 * std::log(b2/b1))/(b1+b2));

        double DwDb1 = dwdb1 + dwdc*dcdb1 + dwdt0*dt0db1;
        double DwDb2 = dwdb2 + dwdc*dcdb2 + dwdt0*dt0db2;

        total_grad0 += -2 * (data - pred) * DwDb1;
        total_grad1 += -2 * (data - pred) * DwDb2;
    }

    return std::make_pair(total_grad0/double(n_samples_per_bin) , total_grad1/double(n_samples_per_bin));
}

FourParamGrad CCMIDSPERegions::GetGradFourParams(double data, double pred, double & c, double & t0, double & b1, double & b2, double &t){
    FourParamGrad s;
    double total_grad0;
    double total_grad1;
    double total_grad2;
    double total_grad3;

    size_t n_samples_per_bin = 40;

    for(size_t idx = 0; idx<n_samples_per_bin; ++idx){
        double t_it = idx * 2.0 / double(n_samples_per_bin) + t;
        double prefactor = -2 * (data - pred) * pred;
        double prefactor2 = -8 * prefactor / (std::exp(-(t_it - t0)/b1) + std::exp( (t_it - t0)/b2));

        total_grad0 += prefactor / c;
        total_grad1 += prefactor2 * (std::exp(-(t_it - t0)/b1)/b1 - std::exp( (t_it - t0)/b2)/b2);
        total_grad2 += prefactor2 *  (t_it - t0)/std::pow(b1,2) * std::exp(-(t_it - t0)/b1);
        total_grad3 += prefactor2 * -(t_it - t0)/std::pow(b2,2) * std::exp( (t_it - t0)/b2);

    }
    s.grad0 = total_grad0/double(n_samples_per_bin);
    s.grad1 = total_grad1/double(n_samples_per_bin);
    s.grad2 = total_grad2/double(n_samples_per_bin);
    s.grad3 = total_grad3/double(n_samples_per_bin);

    return s;
}

double CCMIDSPERegions::PeakLossFunction(const std::vector<double> & xin, std::vector<double> &grad, void * f_data) {

    // typecast to the our struct
    fit_data* d = (fit_data*) f_data;
    int num_params = xin.size();

    if(not (num_params == 2 or num_params == 4)) {
        log_fatal("num_params should be 2 or 4!");
    }

    // figure out input based on constraint
    std::vector<double> x(4);
    if (num_params == 2) {
        x[2] = xin[0];
        x[3] = xin[1];
        x[1] = d->peak_tmax + x[2]*x[3] / (x[2] + x[3]) * std::log(x[2]/x[3]);
        x[0] = d->peak_amplitude * std::pow(std::pow(x[2]/x[3],x[3]/(x[2]+x[3]))
                +std::pow(x[3]/x[2],x[2]/(x[2]+x[3])),8);
    }
    else x = xin;

    // Get the SPE template for these parameters
    // SPETemplate pulse(x[0],x[1],x[2],x[3]);
    // CCMPMTCalibration::DroopedSPETemplate fit_func(pulse);

    // initialize summation parameters
    double squared_residuals = 0;
    double pred, data;
    if (!grad.empty()) {
        for(int i = 0; i < grad.size(); ++i) grad[i] = 0; 
    }

    //pred = fit_func(t); <- this doesnt work,,,not sure why
    // LETS CONVERT UNITS!!!
    double c = x[0];
    double t0 = x[1];
    double b1 = x[2];
    double b2 = x[3];
    double amp = d->peak_amplitude;
    double tmax = d->peak_tmax;

    // loop over each time bin, add to square residuals
    double t = d->t0;
    for(std::vector<short unsigned int>::const_iterator it = d->v_start; it != d->v_end; ++it, t += 2) {
        pred = CCMIDSPERegions::GetPred(c, t0, b1, b2, t);
        data = -(*(d->baseline) + double(*it));

        squared_residuals += std::pow(data - pred, 2);
        // gradient calculated using SPE template function in I3DOMCalibration
        if(!grad.empty()) {
            if(num_params == 4) {
                FourParamGrad s = GetGradFourParams(data, pred, c, t0, b1, b2, t);
                grad[0] += s.grad0;
                grad[1] += s.grad1;
                grad[2] += s.grad2;
                grad[3] += s.grad0;
            }
            else if(num_params == 2) {
                std::pair<double, double> two_param_grad = GetGradTwoParams(amp, data, pred, c, t0, b1, b2, t);
                grad[0] += two_param_grad.first;
                grad[1] += two_param_grad.second;
            }
        }
    }

    return squared_residuals;
}


double CCMIDSPERegions::AmplitudeConstraint(const std::vector<double> &x, std::vector<double> &grad, void * data) {
    double* peak_max = (double*) data;
    if (!grad.empty()) {
        grad[0] = 1;
        grad[1] = 0;
    }

}

SPETemplate CCMIDSPERegions::FitPeak(double peak_amplitude, int64_t length, int64_t time, std::vector<short unsigned int> const & waveform, double baseline, int num_params, SPETemplate two_param_fit_vals) {

    //nlopt::opt opt(nlopt::LD_LBFGS, num_params);
    nlopt::opt opt(nlopt::LN_BOBYQA, num_params);

    int64_t max_size = waveform.size();
    fit_data f_data;
    //f_data.t0 = std::max(long(0),time-pre_window);
    f_data.t0 = 0;
    f_data.baseline = &baseline;
    f_data.v_start = waveform.begin() + time;
    f_data.v_end = waveform.begin() + std::min(time + length, max_size);
    size_t fit_length = std::distance(f_data.v_start, f_data.v_end);

    int64_t peak_tmax = std::distance(f_data.v_start, std::min_element(f_data.v_start,f_data.v_end)) * 2;
    f_data.peak_tmax = peak_tmax;
    f_data.peak_amplitude = peak_amplitude;

    opt.set_min_objective(CCMIDSPERegions::PeakLossFunction, &f_data);

    // Initial guesses, x[0] and x[1] come from minimization of least squares 
    std::vector<double> x(num_params);
    std::vector<double> lb(num_params);
    std::vector<double> ub(num_params);

    if(num_params == 2) {
        double b1 = fit_length * 2;
        double b2 = b1;

        x[0] = b1;
        x[1] = b2;
        // guess some lower and upper bounds on each parameter
        lb[0] = 0.01 * b1;
        ub[0] = 10 * b1;
        lb[1] = 0.01 * b2;
        ub[1] = 10 * b2;
    } else if (num_params == 4) {
        x[0] = two_param_fit_vals.c;
        x[1] = two_param_fit_vals.x0;
        x[2] = two_param_fit_vals.b1;
        x[3] = two_param_fit_vals.b2;
        // guess some lower and upper bounds on each parameter
        lb[0] = 0.01*x[0]; ub[0] = 10*x[0];
        lb[1] = 0.0; ub[1] = fit_length * 2;
        lb[2] = 0.01*x[2]; ub[2] = 10*x[2];
        lb[3] = 0.01*x[3]; ub[3] = 10*x[3];
    }

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    nlopt::result result;
    // now perform the minimization
    double minf;
    try {
        opt.set_xtol_rel(1e-12);
        opt.set_ftol_rel(1e-20);
        opt.set_maxtime(400);
        //opt.set_stopval(-HUGE_VAL);
        result = opt.optimize(x, minf);
    } catch(std::exception &e) {
    }
    //std::cout << "nlopt result = " << result << std::endl;
    double c;
    double t0;
    double b1;
    double b2;
    if(num_params == 2) {
        b1 = x[0];
        b2 = x[1];
        c = peak_amplitude * std::pow(std::pow(b1 / b2, b2 / (b1 + b2))
                + std::pow(b2 / b1, b1 / (b1 + b2)),8);
        t0 = peak_tmax + b1 * b2 / (b1 + b2) * std::log(b1 / b2);
    } else if(num_params == 4) {
        c = x[0];
        t0 = x[1];
        b1 = x[2];
        b2 = x[3];
    }
   // std::cout << "fit vals = " << c << " , " << t0 << " , " << b1 << " , " << b2 << std::endl;
    return SPETemplate(c, t0, b1, b2);
}

void CCMIDSPERegions::ProcessWaveform(WaveformSmoother & smoother,  CCMWaveformUInt16 const & waveform, double const & mode, I3Vector<SPETemplate> & template_per_channel, I3Vector<int64_t> & time_per_channel, I3Vector<int64_t> & length_per_channel, I3Vector<double> & amp_per_channel, I3Vector<short unsigned int> & SPE_wf_to_fit){



    // get the vector of samples from the CCMWaveform object;
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();

    if (samples.size() == 0) {
        return;
    }

    // so now let's call a function to find regions of SPE
    // then we'll fit those regions

    FindRegions(smoother, samples, mode, SPE_wf_to_fit);

    // now let's define a few more quantities and fill using SPE_wf_to_fit
    I3Vector<int> window_amplitude;
    I3Vector<int64_t> window_length;
    I3Vector<int64_t> window_time;

    GetPeakInfo(SPE_wf_to_fit, window_amplitude, window_length, window_time);
    for(size_t ipk = 0; ipk < window_amplitude.size(); ++ipk){
        double peak_amplitude = -(window_amplitude[ipk] - mode);
        double inv_baseline = -1*mode;
        SPETemplate two_param_fit_vals = FitPeak(peak_amplitude, window_length[ipk], window_time[ipk], samples, inv_baseline, 2);
        template_per_channel.push_back(FitPeak(peak_amplitude, window_length[ipk], window_time[ipk], samples, inv_baseline, 4, two_param_fit_vals));
        time_per_channel.push_back(window_time[ipk]);
        length_per_channel.push_back(window_length[ipk]);
        amp_per_channel.push_back(peak_amplitude);
    }
}

void CCMIDSPERegions::Finish() {
}



