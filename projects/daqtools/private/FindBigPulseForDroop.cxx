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

class FindBigPulseForDroop: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    void Geometry(I3FramePtr frame);
    std::tuple<size_t, double, double> CheckForPulse(WaveformSmoother & smoother, double const & mode, size_t start_idx);
    void ProcessWaveform(WaveformSmoother & smoother,  CCMWaveformUInt16 const & waveform, double const & mode, std::vector<short unsigned int> & roi_per_pmt);
    void FindRegions(WaveformSmoother & smoother, std::vector<short unsigned int> const & samples, double const & mode, std::vector<short unsigned int> & roi_per_pmt);

    double smoother_tau = 10.0 ;
    double smoother_delta_t = 2.0;

    public:
    FindBigPulseForDroop(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(FindBigPulseForDroop);

FindBigPulseForDroop::FindBigPulseForDroop(const I3Context& context) : I3Module(context), 
    geometry_name_(""), geo_seen(false) {
        AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    }


void FindBigPulseForDroop::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
}


void FindBigPulseForDroop::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void FindBigPulseForDroop::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms, and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    size_t size = waveforms->size();

    // an I3Map between pmt key and a vector of vectors containing the waveform roi
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<short unsigned int>>> Droop_ROI = boost::make_shared<I3Map<CCMPMTKey, std::vector<short unsigned int>>>();

    // loop over each pmt
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        std::cout << "pmt key = " << key << std::endl;
        BaselineEstimate value = it.second;
        double mode = value.baseline; //mode is positive fyi!!!
        uint32_t channel = pmt_channel_map_[key];

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);
        std::vector<short unsigned int> roi_per_pmt (waveform.GetWaveform().size(), 0);

        WaveformSmoother smoother(waveform.GetWaveform().cbegin(), waveform.GetWaveform().cend(), smoother_delta_t, smoother_tau);
        ProcessWaveform(smoother, waveform, mode, roi_per_pmt);

        // now let's save!
        Droop_ROI->insert({key, roi_per_pmt});
    }

    frame->Put("DroopROI", Droop_ROI);
    PushFrame(frame);
}

std::tuple<size_t, double, double> FindBigPulseForDroop::CheckForPulse(WaveformSmoother & smoother, double const & mode, size_t start_idx){
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
      ) {
        return std::make_tuple(final_idx, max_unsmoothed_value, integral);
    } else {
        return std::make_tuple(start_idx, max_unsmoothed_value, integral);
    }
}


void FindBigPulseForDroop::FindRegions(WaveformSmoother & smoother, std::vector<short unsigned int> const & samples, double const & mode, std::vector<short unsigned int> & roi_per_pmt){
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
    double min_val_threshold = 100;
    //double max_val_threshold = 100;
    //double integral_threshold = 10.0;

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

                    size_t pre_pulse_start_window;
                    size_t pre_pulse_end_window;
                    size_t post_pulse_start_window;
                    size_t post_pulse_end_window;
                    size_t nearby_charge_window = 1000; // don't want charge in the surrounding 2usec!

                    // first checking the pre pulse window for activity
                    if (pulse_first_index < nearby_charge_window){
                        pre_pulse_start_window = 0;
                        pre_pulse_end_window = pulse_first_index;
                    }
                    else{
                        pre_pulse_start_window = pulse_first_index - nearby_charge_window;
                        pre_pulse_end_window = pulse_first_index;
                    }
                    double max_val_pre_pulse = samples[pre_pulse_start_window] - mode;

                    for (size_t pre_pulse_it = pre_pulse_start_window; pre_pulse_it < pre_pulse_end_window; ++pre_pulse_it){
                        max_val_pre_pulse = std::min(max_val_pre_pulse, samples[pre_pulse_it] - mode);
                    }

                    // now checking the post pulse window for activity
                    if (pulse_last_index > (samples.size() - nearby_charge_window)){
                        post_pulse_start_window = pulse_last_index;
                        post_pulse_end_window = samples.size();
                    }
                    else{
                        post_pulse_start_window = pulse_last_index;
                        post_pulse_end_window = pulse_last_index + nearby_charge_window;
                    }
                    double max_val_post_pulse = samples[post_pulse_start_window] - mode;
                    for (size_t post_pulse_it = post_pulse_start_window; post_pulse_it < post_pulse_end_window; ++post_pulse_it){
                        max_val_post_pulse = std::min(max_val_post_pulse, samples[post_pulse_it] - mode);
                    }
                    double surrounding_wf_threshold = 50; // ok if there are SPEs around but don't want any big pulses nearby

                    if (max_val_pre_pulse < surrounding_wf_threshold and max_val_pre_pulse > -1*surrounding_wf_threshold
                        and max_val_post_pulse < surrounding_wf_threshold and max_val_post_pulse > -1*surrounding_wf_threshold) {
                        // made it past our surrounding threshold range!
                        // now let's make sure this is a big pulse
                        if((max_unsmoothed_value >= min_val_threshold) ){
                            // made it past our cuts on adc counts!
                            // let's save waveform to roi_per_pmt
                            for(size_t pulse_it = pre_pulse_start_window; pulse_it <= post_pulse_end_window; ++pulse_it){
                                roi_per_pmt[pulse_it] = samples[pulse_it];
                            }
                        }
                    }
                }

            }

        }
        smoother.Next();
    }
}

void FindBigPulseForDroop::ProcessWaveform(WaveformSmoother & smoother,  CCMWaveformUInt16 const & waveform, double const & mode, std::vector<short unsigned int> & roi_per_pmt){

    // get the vector of samples from the CCMWaveform object;
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();

    if (samples.size() == 0) {
        return;
    }
    // so now let's call a function to find regions of big peaks
    // then we'll fit those regions
    FindRegions(smoother, samples, mode, roi_per_pmt);
}

void FindBigPulseForDroop::Finish() {
}



