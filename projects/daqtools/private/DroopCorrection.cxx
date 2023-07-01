
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

class DroopCorrection: public I3Module {
    bool geo_seen;
    std::string geometry_name_;
    std::string nim_pulses_name_;
    CCMPMTKey bcm_key;
    size_t bcm_channel;
    CCMTriggerKey cosmic_trigger_key;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    // droop correction time constants:
    double tau = 4774.0;
    double delta_t = 2.0;
    void Geometry(I3FramePtr frame);
    void CorrectDroop(std::vector<double> const & samples, double const & mode, I3Vector<double> & droop_corrected_wf_per_channel);
    void ProcessWaveform(CCMWaveformUInt16 const & waveform, double const & mode, I3Vector<double> & droop_corrected_wf_per_channel);
public:
    DroopCorrection(const I3Context&);
    void Configure();
    void DAQ(I3FramePtr frame);
    void Finish();
};

I3_MODULE(DroopCorrection);

DroopCorrection::DroopCorrection(const I3Context& context) : I3Module(context),
    geometry_name_(""), geo_seen(false) {
    AddParameter("CCMGeometryName", "Key for CCMGeometry", std::string(I3DefaultName<CCMGeometry>::value()));
    AddParameter("NIMPulsesName", "Key for NIMLogicPulseSeriesMap", std::string("NIMPulses"));
}


void DroopCorrection::Configure() {
    GetParameter("CCMGeometryName", geometry_name_);
    GetParameter("NIMPulsesName", nim_pulses_name_);
}


void DroopCorrection::Geometry(I3FramePtr frame) {
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

void DroopCorrection::CorrectDroop(std::vector<double> const & samples, double const & mode, I3Vector<double> & droop_corrected_wf_per_channel){

    // this module takes our wf, called samples, and droop corrects and updates droop_corrected_wf_per_channel
    double C = std::exp(-delta_t/tau);
    double B = (1.0 - C);
    double A = tau / delta_t * B;
    double S = 0.0;
    double X = double(samples[0]) / A + B * S;
    droop_corrected_wf_per_channel[0] = X;

    for (size_t i = 1; i < samples.size(); ++i){
        S = X + C * S;
        X = samples[i] / A + B * S;
        droop_corrected_wf_per_channel[i] = X;
    }
}


void DroopCorrection::DAQ(I3FramePtr frame) {
    if(not frame->Has("CCMWaveforms")) {
        throw std::runtime_error("No waveforms!");
    }

    // ptr to vector of all waveforms, derivs, and baselines (one for each channel)
    boost::shared_ptr<const CCMWaveformUInt16Series> waveforms = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>("CCMWaveforms");
    I3Map<CCMPMTKey, BaselineEstimate> const & baseline_mode = frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");

    size_t size = waveforms->size();

    // a vector storing the droop-corrected wfs for each channel
    boost::shared_ptr<I3Vector<I3Vector<double>>> droop_correced_wf(new I3Vector<I3Vector<double>>(size));

    // loop over each pmt
    for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baseline_mode){
        CCMPMTKey key = it.first;
        BaselineEstimate value = it.second;
        double mode = value.baseline; // baseline mode is negative fyi
        uint32_t channel = pmt_channel_map_[key];

        CCMWaveformUInt16 const & waveform = waveforms->at(channel);

        I3Vector<double> droop_corrected_wf_per_channel(waveform.GetWaveform().size());

        ProcessWaveform(waveform, mode, droop_corrected_wf_per_channel);
        droop_correced_wf->operator[](channel) = droop_corrected_wf_per_channel;
    }

    frame->Put("DroopCorrectedWf", droop_correced_wf);
    std::cout << "finished droop correcting wf" << std::endl;
    PushFrame(frame);
}

void DroopCorrection::ProcessWaveform(CCMWaveformUInt16 const & waveform, double const & mode, I3Vector<double> & droop_corrected_wf_per_channel){

    // get the vector of samples from the CCMWaveform object;
    std::vector<short unsigned int> const & samples = waveform.GetWaveform();

    if (samples.size() == 0) {
	    return;
    }
    // let's subtract off the baseline and invert our waveform before passing it to the droop correction function
    std::vector<double> inv_wf (samples.size());

    for (size_t i = 0; i < samples.size(); ++i){
        inv_wf[i] = -1*(double(samples[i]) + mode);
    }

    // let's call our droop correction function and then we're done!
    CorrectDroop(inv_wf, mode, droop_corrected_wf_per_channel);
}

void DroopCorrection::Finish() {
}


