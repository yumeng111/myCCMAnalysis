
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
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include "CCMAnalysis/CCMBinary/BinaryFormat.h"
#include "CCMAnalysis/CCMBinary/BinaryUtilities.h"
#include "icetray/robust_statistics.h"
#include "daqtools/WaveformSmoother.h"
#include "daqtools/WaveformAccumulator.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/CCMBCMSummary.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

class  ElectronicsCorrection: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string calibration_name_;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;

    // initialize droop correction time constant map
    I3Map<CCMPMTKey, double> droop_tau_map_;
    double delta_t = 2.0;
    void Geometry(I3FramePtr frame);
    void Calibration(I3FramePtr frame);
    void BoxFilter(std::vector<double> const & samples, std::vector<double> & box_filter_results);
    void OutlierFilter(std::vector<double> const & samples, std::vector<double> & outlier_filter_results);
    void ECorrection(std::vector<short unsigned int> const & samples, double const & mode, double const & tau,  CCMWaveformDouble & electronics_corrected_wf_per_channel);
public:
     ElectronicsCorrection(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE( ElectronicsCorrection);

 ElectronicsCorrection:: ElectronicsCorrection(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    //AddParameter("CCMCalibrationName", "Key for CCMCalibration", std::string(I3DefaultName<CCMCalibration>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void  ElectronicsCorrection::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("CCMCalibrationName", calibration_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void  ElectronicsCorrection::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    geo_seen = true;
    bool found_bcm_key = false;
    for(auto const & key_om : geo.pmt_geo) {
        if(key_om.second.omtype == CCMOMGeo::OMType::BeamCurrentMonitor) {
            bcm_key = key_om.first;
            found_bcm_key = true;
        }
    }
    if(not found_bcm_key) {
        log_fatal("CCMGeometry does not contain a channel corresponding to a BeamCurrentMonitor");
    }
    bcm_channel = geo.pmt_channel_map.at(bcm_key);
    cosmic_trigger_key = CCMTriggerKey(CCMTriggerKey::TriggerType::CosmicTrigger, 1);
    PushFrame(frame);
}

void ElectronicsCorrection::Calibration(I3FramePtr frame){
    if(not frame->Has(calibration_name_)) {
        log_fatal("Could not find CCMCalibration object with the key named \"%s\" in the Calibration frame.", calibration_name_);
    }

    // CCMCalibration const & calibration = frame->Get<CCMCalibration const>(calibration_name_);
    // ADD TAU TO MAP

}

void ElectronicsCorrection::OutlierFilter(std::vector<double> const & samples, std::vector<double> & outlier_filter_results){

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

void ElectronicsCorrection::BoxFilter(std::vector<double> const & samples, std::vector<double> & box_filter_results){

    for(size_t i = 0; i < samples.size(); ++i){

    if(i == 0){
        // special logic for first element
        box_filter_results[i] = (samples[i] + samples[i+1] + samples[i+2])/3;
    }

    if (i == samples.size() - 1){
        // special logic for last element
        box_filter_results[i] = (samples[i] + samples[i-1] + samples[i-2])/3;
    }

    else{
        box_filter_results[i] = (samples[i-1] + samples[i] + samples[i+1])/3;
    }

    }

}

void ElectronicsCorrection::ECorrection(std::vector<short unsigned int> const & samples, double const & mode, double const & tau,  CCMWaveformDouble & electronics_corrected_wf_per_channel){

    // first we are correcting for the droop!
    std::vector<double> & electronics_correction_samples = electronics_corrected_wf_per_channel.GetWaveform();
    electronics_correction_samples.resize(samples.size());

    // let's start off by inverting samples and subtracting off the baseline
    std::vector<double> inv_wf_min_baseline(samples.size());

    for (size_t wf_it = 0; wf_it < samples.size(); ++wf_it){
        inv_wf_min_baseline[wf_it] = -1 * (double(samples[wf_it]) + mode); // mode is negative and samples is positive!
    }

    // now let's run inv_wf_min_baseline through our box filter
    std::vector<double> box_filter_results(samples.size());
    BoxFilter(inv_wf_min_baseline, box_filter_results);

    // now droop correction time
    double C = std::exp(-delta_t/tau);
    double B = (1.0 - C);
    double A = tau / delta_t * B;
    double S = 0.0;
    double X = double(box_filter_results[0]) / A + B * S;
    electronics_correction_samples[0] = X;

    for (size_t i = 1; i < box_filter_results.size(); ++i){
        S = X + C * S;
        X = box_filter_results[i] / A + B * S;
        electronics_correction_samples[i] = X;
    }

    // now let's subtract off our outlier filter
    // place to store results of the outlier filter
    std::vector<double> outlier_filter_results(samples.size());
    OutlierFilter(electronics_correction_samples, outlier_filter_results);

    // now let's subtract off
    for(size_t it = 0; it < samples.size(); ++it){
        electronics_correction_samples[it] -= outlier_filter_results[it];
    }

}


void  ElectronicsCorrection::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    size_t size = waveforms->size();

    // a vector storing the electronics-corrected wfs for each channel
    boost::shared_ptr<CCMWaveformDoubleSeries> electronics_correced_wf = boost::make_shared<CCMWaveformDoubleSeries>(size);

    // loop over each pmt
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second;
        double mode = value.baseline; // baseline mode is negative fyi
        uint32_t channel = pmt_channel_map_[key];
        double tau = droop_tau_map_[key];

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);

        // get the vector of samples from the CCMWaveform object;
        std::vector<short unsigned int> const & samples = waveform.GetWaveform();

        if (samples.size() == 0) {
	        return;
        }

        // place to store our waveform after electronics correction
        CCMWaveformDouble & electronics_corrected_wf_per_channel = electronics_correced_wf->at(channel);

        // let's pass it to electronics correction module
        // this module first droop corrects
        // then substracts off the outlier filter
        ECorrection(samples, mode, tau, electronics_corrected_wf_per_channel);

        // now let's save our electronics corrected waveform!
        electronics_correced_wf->operator[](channel) = electronics_corrected_wf_per_channel;
    }

    frame->Put("ElectronicsCorrection", electronics_correced_wf);
    std::cout << "finished electronics corrections!" << std::endl;
    PushFrame(frame);
}

void  ElectronicsCorrection::Finish() {
}


