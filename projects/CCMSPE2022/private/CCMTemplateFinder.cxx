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


// struct to hold peak data when fitting 
typedef struct {
  std::vector<short unsigned int>::const_iterator v_start;
  std::vector<short unsigned int>::const_iterator v_end;
  int64_t t0;
  const double* baseline;
  int64_t peak_tmax;
  int64_t peak_amplitude;
} fit_data;


// Fits the SPE template shape to find a waveform template 
// Collects statistics about these regions

class CCMTemplateFinder : public I3Module {
    // Names for keys in the frame
    std::string geometry_name_;
    std::string daq_config_name_;
    std::string waveforms_name_;
    std::string peak_amplitude_name_;
    std::string peak_time_name_;
    std::string peak_length_name_;
    std::string baseline_estimate_name_;
    std::string template_output_name_;


    int64_t fit_length;
    int64_t pre_window;
    int min_amplitude;
    int max_amplitude;
    size_t n_daq_frames = 0;
    static const int num_params = 2;
    

    double AmplitudeConstraint(const std::vector<double> &x,
                               std::vector<double> &grad,
                               void * data);

    SPETemplate FitPeak(double const & peak_amplitude,
                        int64_t const & length,
                        int64_t const & time,
                        std::vector<short unsigned int> const & waveform,
                        double const & baseline);
    static double PeakLossFunction(const std::vector<double> &x,
                                   std::vector<double> &grad,
                                   void * f_data);
    void AddTemplates(I3FramePtr frame);

public:
    CCMTemplateFinder(const I3Context&);
    void Configure();
    void Process();
    void Finish();
};

I3_MODULE(CCMTemplateFinder);

CCMTemplateFinder::CCMTemplateFinder(const I3Context& context) : I3Module(context),
    geometry_name_(""), daq_config_name_(""), waveforms_name_(""), peak_amplitude_name_(""), peak_time_name_(""), peak_length_name_(""), baseline_estimate_name_(""), template_output_name_(""),
    fit_length(10), pre_window(3), min_amplitude(20), max_amplitude(40) {

    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMDAQConfigName", "Key for CCMDAQConfig", std::string(I3DefaultName<CCMAnalysis::Binary::CCMDAQConfig>::value()));
    AddParameter("CCMWaveformsName", "Key to output vector of CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("CCMPeakAmplitudesName", "Key to output vector of peak amplitudes", std::string("PeakAmplitudes"));
    AddParameter("CCMPeakLengthsName", "Key to output vector of peak lengths", std::string("PeakLengths"));
    AddParameter("CCMPeakTimesName", "Key to output vector of peak times", std::string("PeakTimes"));
    AddParameter("BaselineEstimateName", "Key for baseline information of frame", std::string("BaselineEstimates"));
    AddParameter("TemplateOutputName", "The output key of the baseline fit.", std::string("SPETemplate"));
}

void CCMTemplateFinder::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMDAQConfigName", daq_config_name_);
    GetParameter("CCMWaveformsName", waveforms_name_);
    GetParameter("CCMPeakAmplitudesName", peak_amplitude_name_);
    GetParameter("CCMPeakLengthsName", peak_length_name_);
    GetParameter("CCMPeakTimesName", peak_time_name_);
    GetParameter("BaselineEstimateName", baseline_estimate_name_);
    GetParameter("TemplateOutputName", template_output_name_);
}


double CCMTemplateFinder::PeakLossFunction(const std::vector<double> &xin,
                                           std::vector<double> &grad,
                                           void * f_data) {
  
  // typecast to the our struct
  fit_data* d = (fit_data*) f_data;

  // figure out input based on constraint
  std::vector<double> x(4);
  if (num_params==2) {
    x[2] = xin[0];
    x[3] = xin[1];
    x[1] = d->peak_tmax + x[2]*x[3] / (x[2] + x[3]) * std::log(x[2]/x[3]);
    x[0] = d->peak_amplitude * std::pow(std::pow(x[2]/x[3],x[3]/(x[2]+x[3]))
                                     +std::pow(x[3]/x[2],x[2]/(x[2]+x[3])),8);
  }
  else x = xin;
  
  // Get the SPE template for these parameters
  SPETemplate pulse(x[0],x[1],x[2],x[3]);
  CCMPMTCalibration::DroopedSPETemplate fit_func(pulse);
  
  // initialize summation parameters
  double squared_residuals = 0;
  double pred, data;
  if (!grad.empty()) {
    for(int i = 0; i < grad.size(); ++i) grad[i] = 0; 
  }
  int t = d->t0;
  
  // loop over each time bin, add to square residuals
  for(std::vector<short unsigned int>::const_iterator it = d->v_start; it != d->v_end; ++it, ++t) {
    pred = fit_func(t);
    data = -(*(d->baseline) + double(*it));
    squared_residuals += std::pow(data - pred, 2);
    // gradient calculated using SPE template function in I3DOMCalibration
    if (!grad.empty()) {
      double prefactor = -2 * (data - pred) * pred;
      double prefactor2 = -8 * prefactor / (std::exp(-(t - x[1])/x[2]) + 
                                            std::exp( (t - x[1])/x[3]));
      if (num_params==4) {
        grad[0] += prefactor / x[0];
        grad[1] += prefactor2 * (std::exp(-(t - x[1])/x[2])/x[2] - 
                                 std::exp( (t - x[1])/x[3])/x[3]);
        grad[2] += prefactor2 *  (t - x[1])/std::pow(x[2],2) * std::exp(-(t - x[1])/x[2]);
        grad[3] += prefactor2 * -(t - x[1])/std::pow(x[3],2) * std::exp( (t - x[1])/x[3]);
      }
      else if (num_params==2) {
        grad[0] += prefactor2 *  (t - x[1])/std::pow(x[2],2) * std::exp(-(t - x[1])/x[2]);
        grad[1] += prefactor2 * -(t - x[1])/std::pow(x[3],2) * std::exp( (t - x[1])/x[3]);
      }
    }
  }
  return squared_residuals;
}

double CCMTemplateFinder::AmplitudeConstraint(const std::vector<double> &x,
                                              std::vector<double> &grad,
                                              void * data) {
  double* peak_max = (double*) data;
  if (!grad.empty()) {
    grad[0] = 1;
    grad[1] = 0;
  }

}

SPETemplate CCMTemplateFinder::FitPeak(double const & peak_amplitude,
                                       int64_t const & length,
                                       int64_t const & time,
                                       std::vector<short unsigned int> const & waveform,
                                       double const & baseline) {
  
  nlopt::opt opt(nlopt::LD_MMA, num_params);
  fit_data f_data;
  f_data.t0 = std::max(long(0),time-pre_window);
  f_data.baseline = &baseline;
  f_data.v_start = waveform.begin() + f_data.t0;
  f_data.v_end = waveform.begin() + time + std::min(fit_length,length);
  int64_t peak_tmax = f_data.t0 + std::distance(f_data.v_start, std::min_element(f_data.v_start,f_data.v_end));
  f_data.peak_tmax = peak_tmax;
  f_data.peak_amplitude = peak_amplitude;
  
  opt.set_min_objective(CCMTemplateFinder::PeakLossFunction, &f_data);
  
  // Initial guesses, x[0] and x[1] come from minimization of least squares 
  std::vector<double> x(4);
  x[2] = double(time + std::min(fit_length,length) - f_data.t0)/2;
  x[3] = x[2];
  x[1] = peak_tmax + x[2]*x[3] / (x[2] + x[3]) * std::log(x[2]/x[3]);
  x[0] = peak_amplitude * std::pow(std::pow(x[2]/x[3],x[3]/(x[2]+x[3]))
                                  +std::pow(x[3]/x[2],x[2]/(x[2]+x[3])),8);
  
  // guess some lower and upper bounds on each parameter
  std::vector<double> lb(4);
  std::vector<double> ub(4);
  lb[0] = 0.01*x[0]; ub[0] = 10*x[0];
  lb[1] = f_data.t0; ub[1] = f_data.t0 + std::distance(f_data.v_start,f_data.v_end);
  lb[2] = 0.01*x[2]; ub[2] = 10*x[2];
  lb[3] = 0.01*x[3]; ub[3] = 10*x[3];

  if(num_params==2) {
    x = std::vector<double>(x.begin() + 2, x.end());
    lb = std::vector<double>(lb.begin() + 2, lb.end());
    ub = std::vector<double>(ub.begin() + 2, ub.end());
  }
  
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  // now perform the minimization
  double minf;
  try{
    opt.set_xtol_rel(1e-4);
    nlopt::result result = opt.optimize(x, minf);
  }
  catch(std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }
  std::vector<double> xp(4);
  if(num_params==2) {
    xp[2] = x[0];
    xp[3] = x[1];
    xp[1] = peak_tmax + xp[2]*xp[3] / (xp[2] + xp[3]) * std::log(xp[2]/xp[3]);
    xp[0] = peak_amplitude * std::pow(std::pow(xp[2]/xp[3],xp[3]/(xp[2]+xp[3]))
                                     +std::pow(xp[3]/xp[2],xp[2]/(xp[2]+xp[3])),8);
  }
  else xp = x; 
  
  return SPETemplate(xp[0],xp[1],xp[2],xp[3]);
}

void CCMTemplateFinder::AddTemplates(I3FramePtr frame) {
    
   
    // let's grab some stuff from the frame
    CCMWaveformUInt16Series const & waveforms = frame->Get<CCMWaveformUInt16Series>(waveforms_name_);
    size_t size = waveforms.size();
    I3Vector<I3Vector<int>> const & window_amplitudes = frame->Get<I3Vector<I3Vector<int>>>(peak_amplitude_name_);
    I3Vector<I3Vector<int64_t>> const & window_lengths = frame->Get<I3Vector<I3Vector<int64_t>>>(peak_length_name_);
    I3Vector<I3Vector<int64_t>> const & window_times = frame->Get<I3Vector<I3Vector<int64_t>>>(peak_time_name_);
    I3Vector<double> const & baseline_estimates = frame->Get<I3Vector<double>>(baseline_estimate_name_);
    
    // a vector storing the template for each channel
    boost::shared_ptr<I3Vector<I3Vector<SPETemplate>>> channel_templates(new I3Vector<I3Vector<SPETemplate>>(size));
    // a vector storing the time for each template
    boost::shared_ptr<I3Vector<I3Vector<int64_t>>> channel_times(new I3Vector<I3Vector<int64_t>>(size));
     
    for(size_t ich = 0; ich < waveforms.size(); ++ich) {
      I3Vector<SPETemplate> templates;
      I3Vector<int64_t> times;
      for(size_t ipk = 0; ipk < window_amplitudes[ich].size(); ++ipk) {
        double peak_amplitude = -(window_amplitudes[ich][ipk] + baseline_estimates[ich]);
        if(peak_amplitude < min_amplitude || peak_amplitude > max_amplitude) continue; 
        templates.push_back(FitPeak(peak_amplitude,
                                    window_lengths[ich][ipk],
                                    window_times[ich][ipk],
                                    waveforms[ich].GetWaveform(),
                                    baseline_estimates[ich]));
        times.push_back(window_times[ich][ipk]);
      }
      channel_templates->operator[](ich) = templates;
      channel_times->operator[](ich) = times;
    }
    frame->Put("SPETemplates", channel_templates);
    frame->Put("SPETemplatesTimes", channel_times);
}


void CCMTemplateFinder::Process() {
    I3FramePtr frame = PopFrame();

    if(frame->GetStop() == I3Frame::Geometry) {
        PushFrame(frame);
        return;
    }

    if(frame->GetStop() != I3Frame::DAQ) {
        PushFrame(frame);
        return;
    }

    n_daq_frames += 1;
    AddTemplates(frame);
    PushFrame(frame);

}

void CCMTemplateFinder::Finish() {
}
