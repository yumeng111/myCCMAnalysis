
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
#include <dataclasses/calibration/CCMCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

class  ElectronicsCorrection: public I3Module {
    bool geo_seen;
    bool calib_seen;
    std::string geometry_name_;
    std::string ccm_calibration_name_;
    std::string ccm_waveforms_name_;
    std::string baseline_estimates_name_;
    std::string output_name_;
    std::vector<double> default_droop_tau_;
    CCMPMTCalibration default_calib_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    CCMCalibration ccm_calibration_;

    // initialize droop correction time constant map
    double delta_t = 2.0;
    void Geometry(I3FramePtr frame);
    void Calibration(I3FramePtr frame);
    void BoxFilter(std::vector<double> const & samples, std::vector<double> & box_filter_results);
    void OutlierFilter(std::vector<double> const & samples, std::vector<double> & outlier_filter_results);
    void ECorrection(std::vector<uint16_t> const & samples, double const & mode, const CCMPMTCalibration& calibration,  std::vector<double> & electronics_correction_samples);
public:
     ElectronicsCorrection(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE( ElectronicsCorrection);

 ElectronicsCorrection:: ElectronicsCorrection(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false), calib_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("CCMCalibrationName", "Key for CCMCalibration", std::string("CCMCalibration"));
    AddParameter("CCMWaveformsName", "Key for input CCMWaveforms", std::string("CCMWaveforms"));
    AddParameter("BaselineEstimatesName", "Key for input BaselineEstimates", std::string("BaselineEstimates"));
    AddParameter("DefaultDroopTau", "The default droop time constant (in ns) to use if there is no entry in the calibratioin", std::vector((double)1200, (double)5000));
    AddParameter("OutputName", "Key to save output CCMWaveformDoubleSeries to", std::string("CCMCalibratedWaveforms"));
}


void  ElectronicsCorrection::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMCalibrationName", ccm_calibration_name_);
    GetParameter("CCMWaveformsName", ccm_waveforms_name_);
    GetParameter("BaselineEstimatesName", baseline_estimates_name_);
    GetParameter("DefaultDroopTau", default_droop_tau_);
    GetParameter("OutputName", output_name_);
    default_calib_.SetDroopTimeConstant(default_droop_tau_);
}


void  ElectronicsCorrection::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    PushFrame(frame);
}

void ElectronicsCorrection::Calibration(I3FramePtr frame) {
    if(not frame->Has(ccm_calibration_name_)) {
        log_fatal("Could not find CCMCalibration object with the key named \"%s\" in the Calibration frame.", ccm_calibration_name_.c_str());
    }
    ccm_calibration_ = frame->Get<CCMCalibration>(ccm_calibration_name_);
    calib_seen = true;
    PushFrame(frame);
}

void ElectronicsCorrection::OutlierFilter(std::vector<double> const & samples, std::vector<double> & outlier_filter_results) {

    // first let's find the average of the first 10 bins of our wf as the starting value
    double delta_tau = 20;
    double prev_tau = 2.0;
    double next_tau = 2.0;
    double starting_val = 0;
    double counter = 0;

    for(size_t it = 0; it < 10; ++it) {
        starting_val += samples[it];
        counter += 1;
    }

    double value = starting_val/counter;

    // now let's loop over the waveform
    for(size_t wf_it = 0; wf_it < samples.size(); ++wf_it) {
        double delta = samples[wf_it] - value;
        double e = std::fabs(delta) / delta_tau;

        if(wf_it > 0) {
            e += std::fabs(samples[wf_it] - samples[wf_it - 1]) / prev_tau;
        }

        if(wf_it + 1 < samples.size()) {
            e += std::fabs(samples[wf_it] - samples[wf_it + 1]) / next_tau;
        }
        delta *= std::exp(-e);
        value += delta;
        outlier_filter_results[wf_it] = value;
    }

}

void ElectronicsCorrection::BoxFilter(std::vector<double> const & samples, std::vector<double> & box_filter_results) {
    for(size_t i = 0; i < samples.size(); ++i) {
            // special logic for first element
        if(i == 0)
            box_filter_results[i] = (samples[i] + samples[i+1] + samples[i+2])/3;
        else if (i == samples.size() - 1)
            // special logic for last element
            box_filter_results[i] = (samples[i] + samples[i-1] + samples[i-2])/3;
        else
            box_filter_results[i] = (samples[i-1] + samples[i] + samples[i+1])/3;
    }

}

void ElectronicsCorrection::ECorrection(std::vector<uint16_t> const & samples, double const & mode, const CCMPMTCalibration& calibration,  std::vector<double> & electronics_correction_samples) {

    // let's start off by inverting samples and subtracting off the baseline
    std::vector<double> inv_wf_min_baseline(samples.size());

    for (size_t wf_it = 0; wf_it < samples.size(); ++wf_it) {
        inv_wf_min_baseline[wf_it] = -1 * (double(samples[wf_it]) + mode); // mode is negative and samples is positive!
    }

    // now let's run inv_wf_min_baseline through our box filter
    std::vector<double> box_filter_results(samples.size());
    BoxFilter(inv_wf_min_baseline, box_filter_results);

    // now droop correction time
    std::vector<double> tau = calibration.GetDroopTimeConstant();
    if(std::isnan(tau.at(0))) {
        tau = default_droop_tau_;
    }

    // now loop over tau
    
    double C;
    double B;
    double A;
    double S;
    double X;

    // place to store our results from droop correcting once
    std::vector<double> first_droop_correction_results(box_filter_results.size());

    for (size_t tau_it = 0; tau_it < 2; ++tau_it){
    
        // define our numbers on the first iteration
        C = std::exp(-delta_t/tau.at(tau_it));
        B = (1.0 - C);
        A = tau.at(tau_it) / delta_t * B;
        S = 0.0;

        if (tau_it == 0){
            X = double(box_filter_results[0]) / A + B * S;
            first_droop_correction_results[0] = X;
            
            for (size_t i = 1; i < box_filter_results.size(); ++i) {
                S = X + C * S;
                X = box_filter_results[i] / A + B * S;
                first_droop_correction_results[i] = X;
            }
        }
        
        if (tau_it == 1){
            X = double(first_droop_correction_results[0]) / A + B * S;
            electronics_correction_samples.push_back(X);
            
            for (size_t i = 1; i < box_filter_results.size(); ++i) {
                S = X + C * S;
                X = first_droop_correction_results[i] / A + B * S;
                electronics_correction_samples.push_back(X);
            }
        }

    
    }

    // now let's subtract off our outlier filter
    // place to store results of the outlier filter
    std::vector<double> outlier_filter_results(samples.size());
    OutlierFilter(electronics_correction_samples, outlier_filter_results);

    // now let's subtract off
    for(size_t it = 0; it < samples.size(); ++it) {
        electronics_correction_samples[it] -= outlier_filter_results[it];
    }
}


void  ElectronicsCorrection::DAQ(I3FramePtr frame) {
    if(not geo_seen) {
        log_fatal("No Geometry frame seen before DAQ frame! Did you forget to include a geometry file?");
    }
    if(not calib_seen) {
        log_fatal("No Calibration frame seen before DAQ frame! Did you forget to inlcude a calibration file?");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    CCMCalibration const & calibration = ccm_calibration_;

    boost::shared_ptr<CCMWaveformUInt16Series const> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>(ccm_waveforms_name_);
    boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const> baseline_mode = frame->Get<boost::shared_ptr<I3Map<CCMPMTKey, BaselineEstimate> const>>(baseline_estimates_name_);
    if(waveforms == nullptr) {
        log_fatal("No CCMWaveformUInt16Series under key name \"%s\"", ccm_waveforms_name_);
    }
    if(baseline_mode == nullptr) {
        log_fatal("No I3Map<CCMPMTKey, BaselineEstimate> under key name \"%s\"", baseline_estimates_name_);
    }

    size_t size = waveforms->size();

    // a vector storing the electronics-corrected wfs for each channel
    boost::shared_ptr<CCMWaveformDoubleSeries> electronics_corrected_wf = boost::make_shared<CCMWaveformDoubleSeries>(size);

    // loop over each pmt
    // We loop over pmt keys in baseline_estimates because a baseline estimate is required for the electronics correction
    // and the baseline_estimates object should have its pmt keys derived from the same geometry
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : *baseline_mode) {
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second;
        double mode = value.baseline; // baseline mode is negative fyi
        uint32_t channel = pmt_channel_map_.at(key);

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);

        // get the vector of samples from the CCMWaveform object;
        std::vector<uint16_t> const & samples = waveform.GetWaveform();

        // no work to be done if there is not waveform available
        if (samples.size() == 0) {
	        return;
        }

        // place to store our waveform after electronics correction
        std::vector<double> & electronics_correction_samples = electronics_corrected_wf->at(channel).GetWaveform();
        electronics_correction_samples.reserve(samples.size());

        // let's pass it to electronics correction module
        // this module first droop corrects
        // then substracts off the outlier filter
        if(calib == calibration.pmtCal.cend()) {
            ECorrection(samples, mode, default_calib_, electronics_correction_samples);
        } else {
            ECorrection(samples, mode, calib->second, electronics_correction_samples);
        }
    }

    frame->Put(output_name_, electronics_corrected_wf);
    PushFrame(frame);
}

void  ElectronicsCorrection::Finish() {
}


