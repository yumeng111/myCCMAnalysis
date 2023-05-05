#include <icetray/IcetrayFwd.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/I3PODHolder.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/I3Orientation.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"


class CCMNewBaselines : public I3Module {
    void AddBaselineStats(I3FramePtr frame);
    void ProcessWaveform(CCMWaveformUInt16 const & waveform, I3Vector<double> & baselines, I3Vector<double> & baselines_times, I3Vector<double>& derivs, I3Vector<double>& smoothed_wf);
    I3Vector<double> DoAllSmoothing(I3Vector<double> const & samples) const;
    I3Vector<double> ComputeDerivativeOfWaveform(I3Vector<double> const & smooth_wf) const;
    std::pair<I3Vector<double>, I3Vector<double>> FindZeroDeriv(I3Vector<double> const & derivs, I3Vector<double> const & smooth_wf) const;
    bool is_first_frame;
    double previous_frame_info;
public:
    CCMNewBaselines(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(CCMNewBaselines);

CCMNewBaselines::CCMNewBaselines(const I3Context& context) : I3Module(context) {
    // AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
}

void CCMNewBaselines::Configure() {
    // GetParameter("CCMGeometryName", geometry_name_);
}

void CCMNewBaselines::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms (one for each channel)
    // PYTHON: waveforms = frame.Get("CCMWaveforms")
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");

    // Place to store baselines
    //boost::shared_ptr<CCMWaveformUInt16Series> baselines = boost::make_shared<CCMWaveformUInt16Series>(waveforms->size());
    boost::shared_ptr<I3Vector<I3Vector<double>>> baselines (new I3Vector<I3Vector<double>>(waveforms->size()));
    boost::shared_ptr<I3Vector<I3Vector<double>>> baselines_times (new I3Vector<I3Vector<double>>(waveforms->size()));
    boost::shared_ptr<I3Vector<I3Vector<double>>> derivs (new I3Vector<I3Vector<double>>(waveforms->size()));
    boost::shared_ptr<I3Vector<I3Vector<double>>> smoothed_wf (new I3Vector<I3Vector<double>>(waveforms->size()));
    
    // loop over each channel / each waveform
    for(size_t i=0; i<waveforms->size(); ++i) {

        // get the CCMWaveform object from the vector
        // PYTHON: waveform = waveforms[i]
        CCMWaveformUInt16 const & waveform = waveforms->at(i);
        
	// pass the waveform to the function that computes the baselines
	ProcessWaveform(waveform, baselines->at(i), baselines_times->at(i), derivs->at(i), smoothed_wf->at(i));
    }
    
    frame->Put("Baselines", baselines);
    //frame->Put("BaselinesTimes", baselines_times);
    //frame->Put("Derivatives", derivs);
    //frame->Put("SmoothedWaveforms", smoothed_wf);
    PushFrame(frame);
    std::cout << "finished saving baselines" << std::endl;
}

I3Vector<double> CCMNewBaselines::DoAllSmoothing(I3Vector<double> const & samples) const {
    // Do the smoothing...
    // first thing, find inverse wf
    I3Vector<double> mult_factor(samples.size(), -1.0);
    I3Vector<double> inverse_wf(samples.size());

    for (int i = 0; i < samples.size(); i++) {
        inverse_wf[i] = samples[i] * mult_factor[i];
    }

    // next up, exponential smoothing filter
    double tau = 10.0;
    double delta_t = 2.0;

    I3Vector<double> exp_smooth_wf(inverse_wf.size());
    double alpha = exp((-1 * delta_t) / tau);
    double y_i = inverse_wf[0];

    for (int i = 0; i < inverse_wf.size(); i++) {
        y_i += (1.0 - alpha) * (inverse_wf[i] - y_i);
        exp_smooth_wf[i] = y_i;
    }

    // finally, smoothing filter
 
    I3Vector<double> smooth_wf(exp_smooth_wf.size());

    for (int i = 0; i < exp_smooth_wf.size(); i++) {
        if (i == 0) { // first element, special logic
            double smooth = (exp_smooth_wf[i] + exp_smooth_wf[i+1] + exp_smooth_wf[i+2])/3;
            smooth_wf[i] = smooth;
        }

        if (i == exp_smooth_wf.size() - 1) { // last element, special logic
            double smooth = (exp_smooth_wf[i-2] + exp_smooth_wf[i-1] + exp_smooth_wf[i])/3;
            smooth_wf[i] = smooth;
        }

        if (i > 0 && i < exp_smooth_wf.size() - 1) {
            double smooth = (exp_smooth_wf[i-1] + exp_smooth_wf[i] + exp_smooth_wf[i+1])/3; // smoothing over three values centered around i
            smooth_wf[i] = smooth;
        }
    }
    return smooth_wf; //return smoothed wf


}

I3Vector<double> CCMNewBaselines::ComputeDerivativeOfWaveform(I3Vector<double> const & smooth_wf) const {
    // compute the finite difference
    I3Vector<double> derivs(smooth_wf.size(), 0.0);

    for (size_t i = 0; i < smooth_wf.size(); ++i)
    {
        if (i == 0)
        {
            // first element, cannot find backwards deriv so let's just find the forwards deriv
            double rise = smooth_wf[i + 1] - smooth_wf[i];
            double run = 2.0;
            derivs[i] = rise / run;
        }
        else
        {
            double rise = smooth_wf[i] - smooth_wf[i - 1];
            double run = 2.0; // keeping the time in units of nsec
            derivs[i] = rise / run;
        }
    }

    return derivs;

}


std::pair<I3Vector<double>, I3Vector<double>> CCMNewBaselines::FindZeroDeriv(I3Vector<double> const & derivs, I3Vector<double> const & smooth_wf) const {
    // let's find region of zero deriv
    double sample_len=20;
    double tolerance = 0.30;

    I3Vector<double> baselines_for_plotting;
    I3Vector<double> baseline_times_for_plotting;


    int counter = 0;

    for (size_t i = 0; i < derivs.size() - 1; ++i)
    {
        if (derivs[i] < tolerance && derivs[i] > -1.0 * tolerance) // within tolerance at bin i
        {
            counter++;

            if (derivs[i+1] > tolerance || derivs[i+1] < -1.0 * tolerance) // within tolerance at bin i but NOT at bin i+1
            {
                if (counter >= sample_len)
                {
                    for (int j = counter - 1; j >= 0; --j)
                    {
                        baselines_for_plotting.push_back(smooth_wf[i-j]);
                        baseline_times_for_plotting.push_back(2.0 * (i-j));
                    }
                }

                counter = 0;
            }
        }
    }

    return std::make_pair(baselines_for_plotting, baseline_times_for_plotting);

    //so we have our vectors of times and sample for the baselines
    // let's interpolate between the baselines

    //I3Vector<double> & baselinedouble nsamples = smooth_wf.size();
    //I3Vector<double> final_baselines(smooth_wf.size());
    //double nsamples_counter = 0.0;
    //std::cout << "wf size = " << nsamples << std::endl;

    //for (size_t i = 0; i < baseline_times_for_plotting.size()-1; ++i){
    //    double current_time = baseline_times_for_plotting[i];
    //    double next_time = baseline_times_for_plotting[i+1];
    //    double nbins_to_add = 0; //offset that will be dertmined if baselines don't start at time equal zero
    //    double n_to_interpolate = 0;
    //    double x0 = 0;
    //    double y0 = 0;
    //    double x1 = 0;
    //    double y1 = 0;
    //    double x = 0;
    //    double y = 0;
    //    double bins_to_end = 0;

    //    //first let's check to see if the baselines start at time zero
    //    //if they don't, we add in the baselines 
    //    if (i==0 and current_time!= 0){
    //       nbins_to_add = current_time/2.0;
    //       for (size_t j = 0; j<nbins_to_add; ++j){
    //           final_baselines[i+j] = baselines_for_plotting[i];
    //           nsamples_counter += 1.0;
    //       }
    //    }
    //    
    //    //so the easy case, next time is one bin greater than current time
    //    
    //    if (next_time == (current_time + 2.0) and i!= 0){
    //       final_baselines[i+nbins_to_add] = baselines_for_plotting[i];
    //       nsamples_counter += 1.0;

    //    }

    //    //now the interpolation case
    //    if (next_time > (current_time + 2.0) and i!=0){
    //       n_to_interpolate = (next_time - current_time)/2;
    //       x0 = current_time;
    //       y0 = baselines_for_plotting[i];

    //       x1 = next_time; 
    //       y1 = baselines_for_plotting[i+1];

    //       for (size_t k = 0; k < n_to_interpolate; ++k){
    //         x = baselines_for_plotting[i+k];
    //         y = y0 + (x-x0)*(y1-y0)/(x1-x0);
    //         final_baselines[i+nbins_to_add+k] = y;
    //         nsamples_counter += 1.0;
    //       } 
    //    }

    //    //now some special logic for the last time bin
    //    //first the case where the baselines time vector reaches nsamples
    //    //second the case where it doesnt
    //    
    //    if (i == baseline_times_for_plotting.size()-2){
    //      if (next_time == 2*nsamples){
    //        final_baselines[nsamples-1] = baselines_for_plotting[-1];
    //        nsamples_counter += 1.0;
    //      }

    //      if (next_time < 2*nsamples){
    //         bins_to_end = (2*nsamples - next_time)/2.0;
    //         for (size_t l = 0; l < bins_to_end+1; ++l){
    //            final_baselines[nsamples-l-1] = baselines_for_plotting[-1];
    //            nsamples_counter += 1.0;
    //         }
    //      }
    //    
    //    }
    //}
    ////std::cout<< "all done w finding baselines!" << std::endl;
    ////std::cout << "final number of samples = " << nsamples << " and nbaselines = " << nsamples_counter << std::endl;
    //nsamples_counter = 0.0;
    //std::cout<<"found baselines"<<std::endl;
    //return final_baselines;
}

void CCMNewBaselines::ProcessWaveform(CCMWaveformUInt16 const & waveform, I3Vector<double>& baselines, I3Vector<double>& baselines_times, I3Vector<double>& derivs, I3Vector<double>& smoothed_wf) {
    // get the vector of samples from the CCMWaveform object;
    // PYTHON: samples = waveform.waveform
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();
    I3Vector<double> samples_double(samples.size());
   
    if (samples.size() == 0) {
    	I3Vector<double> interp_baselines_empty (8000, 0.0);
	baselines = interp_baselines_empty;
	std::cout << "oops! empty frame" << std::endl; 
    	return;
    }

    for (size_t i = 0; i < samples.size(); ++i) {
    samples_double[i] = static_cast<double>(samples[i]);
    } 
    
    //I3Vector<double> smoothed_samples(samples.size());
    //I3Vector<double> derivative_samples(samples.size());

    smoothed_wf = DoAllSmoothing(samples_double);
    derivs = ComputeDerivativeOfWaveform(smoothed_wf);
    
    std::pair<I3Vector<double>, I3Vector<double>> vecs = FindZeroDeriv(derivs, smoothed_wf);
    
    baselines = vecs.first;
    baselines_times = vecs.second; 

    // let's do some linear interpolation of baselines
    double baseline_size = baselines.size();
    double goal_size = samples_double.size();

    // first case where baselines is empty
    // fill baselines with 0

    if (baseline_size == 0){
       I3Vector<double> interp_baselines_empty (goal_size, 0.0);
       baselines = interp_baselines_empty;
    }

    // now that case where baseline_size < goal_size
    // time to interpolate

    if (baseline_size < goal_size && baseline_size > 0){
       I3Vector<double> interp_baselines (goal_size);
       // let's loop over the times
       
       // let's check to see at what time baselines start at
       double start_time = baselines_times[0];
       double beginning_offset = start_time/2.0;

       if (start_time != 0){
	  for (size_t i = 0; i < beginning_offset; ++i){
	      interp_baselines[i] = baselines[0];
	  }
       }

       // let's check to see at what time baselines end
       double end_time = baselines_times.back();
       double end_offset = (2.0*goal_size - end_time)/2;
       
       if (end_time < 2.0*goal_size){
	  for (size_t i = 0; i < end_offset; ++i){
	      interp_baselines[goal_size-i] = baselines.back();
	  }
       }
 
       // now let's do our normal interpolation
       
       for(size_t i = 0; i < baseline_size; ++i){

          if(i == baseline_size-1){
	    interp_baselines[goal_size-end_offset-1] = baselines[i]; 
	  }

          if(i < baseline_size-1){
            double current_time = baselines_times[i];
	    double next_time = baselines_times[i+1];

	    //so the easy case, next time is one bin greater than current time
    
            if (next_time == current_time + 2.0){
	       double new_time = current_time/2;
               interp_baselines[new_time] = baselines[i];
               } 

            //now the interpolation case
            if (next_time > current_time + 2.0){
               double n_to_interpolate = (next_time - current_time)/2;
               double x0 = current_time;
               double y0 = baselines[i];

               double x1 = next_time;
               double y1 = baselines[i+1];

               for (size_t k = 0; k < n_to_interpolate; ++k){
		 double new_time = current_time/2;
                 double x = current_time + 2*k; 
                 double y = y0 + (x-x0)*(y1-y0)/(x1-x0);
                 interp_baselines[new_time+k] = y;
               }
            }


	    }
       }

    baselines = interp_baselines;

    }

}

void CCMNewBaselines::Finish() {
}
