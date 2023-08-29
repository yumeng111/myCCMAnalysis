#include <cctype>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <cholmod.h>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/I3ConditionalModule.h>

#include "wavereform/CCMWavereform.h"
#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

#include <dataclasses/CCMPMTFunctions.h>
#include <dataclasses/I3TimeWindow.h>
#include <dataclasses/calibration/CCMCalibration.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/status/CCMDetectorStatus.h>
#include <dataclasses/status/CCMPMTStatus.h>


I3_MODULE(CCMWavereform);

CCMWavereform::CCMWavereform(const I3Context &ctx) : I3ConditionalModule(ctx)
{
	waveform_name_ = "CCMWaveforms";
	AddParameter("Waveforms", "Name of calibrated waveforms in the frame", waveform_name_);

    pulse_name_ = "WavedeformPulses";
	AddParameter("Pulses", "Name of the pulse series in the frame", pulse_name_);

	chi_name_ = "RecoPulseChi";
	AddParameter("Chi", "Where to store the chi^2 between waveform and pulses", chi_name_);

	chi_threshold_ = 2e3;
	AddParameter("ChiThreshold", "Flag any PMT with a waveform/pulse chi^2 greater than this.", chi_threshold_);

	flag_name_ = "LargeChiPMTs";
	AddParameter("Flag", "Where to store the list of bad PMTs", flag_name_);

    geometry_name_ = "CCMGeometry";
    AddParameter("CCMGeometryName", "Key for CCMGeometry", geometry_name_);
}


void CCMWavereform::Configure(){
    GetParameter("CCMGeometryName", geometry_name_);
	GetParameter("Waveforms", waveform_name_);
	GetParameter("Pulses", pulse_name_);
	GetParameter("Chi", chi_name_);
	GetParameter("ChiThreshold", chi_threshold_);
	GetParameter("Flag", flag_name_);
}

void CCMWavereform::Calibration(I3FramePtr frame){
	PushFrame(frame);
}

void CCMWavereform::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_);
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    PushFrame(frame);
}

void CCMWavereform::GetChi(double & chi, double & dof, std::vector<double> const & refolded_wf, std::vector<double> const & samples){
    chi = 0;
    dof = 0;
    // loop over each bin and compare refolded_wf to wf
    for(size_t wf_it = 0; wf_it < samples.size(); ++wf_it){
        if (samples[wf_it] < 7){
            chi += 0.0;
        }
        else {
            chi += (std::pow(samples[wf_it] - refolded_wf[wf_it], 2));
            dof += 1;
        }
    }
}

void CCMWavereform::GetRefoldedWf(std::vector<double> & refolded_wf, std::vector<CCMRecoPulse> const & pulses, CCMPMTCalibration const & calib){
    // time for how far we need to calculate our pulse out to
    // should replace with paramters from calibration one day
    double pulse_time;
    double pulse_charge;

    // now let's loop over pulses
    for (size_t pulse_it = 0; pulse_it < pulses.size(); ++pulse_it){
        pulse_time = pulses[pulse_it].GetTime();
        pulse_charge = pulses[pulse_it].GetCharge();
        // now let's figure out how many ADC counts this pulse is
        CCMSPETemplate channel_template = calib.GetSPETemplate();
        int start_bin = int((pulse_time - 20.0) / 2.0);
        int end_bin = int((pulse_time + 120.0) / 2.0);
        start_bin = std::max(int(0), start_bin);
        end_bin = std::min(int(refolded_wf.size()) - 1, end_bin+1);

        // let's loop over all the bins in this pulse
        int i = 0;
        for(int bin_it = start_bin; bin_it < end_bin; ++bin_it, ++i){
            double time = i * 2 - 20 + 2; // fitting to the left edge of the data
            double wf_value = pulse_charge * (channel_template.Evaluate(time));
            refolded_wf[bin_it] += wf_value;
        }
    }
}


void CCMWavereform::DAQ(I3FramePtr frame) {
    const CCMCalibration& calibration_ = frame->Get<CCMCalibration>("CCMCalibration");
    I3Map<CCMPMTKey, BaselineEstimate> const & baselines= frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");
	CCMRecoPulseSeriesMapConstPtr pulse_map = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulse_name_);
    boost::shared_ptr<const CCMWaveformUInt16Series> waveform_map = frame->Get<boost::shared_ptr<const CCMWaveformUInt16Series>>(waveform_name_);
    boost::shared_ptr<const CCMWaveformDoubleSeries> corrected_waveform_map = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>("CCMCalibratedWaveforms");

    // let's make places to store our chis, refolded wfs, and our bad pmts
    boost::shared_ptr<I3Map<CCMPMTKey, double>> chi_map = boost::make_shared<I3Map<CCMPMTKey, double>>();
    boost::shared_ptr<I3Map<CCMPMTKey, double>> dof_map = boost::make_shared<I3Map<CCMPMTKey, double>>();
    boost::shared_ptr<I3Map<CCMPMTKey, std::vector<double>>> refolded_wf_map = boost::make_shared<I3Map<CCMPMTKey, std::vector<double>>>();
    boost::shared_ptr<I3Map<CCMPMTKey, double>> screwy_pmts = boost::make_shared<I3Map<CCMPMTKey, double>>();

    for(std::pair<CCMPMTKey const, uint32_t> const & p : pmt_channel_map_) {
        CCMPMTKey key = p.first;

        if(calibration_.pmtCal.count(key) == 0){
            std::cout << "oops! no calibration for " << key << std::endl;
            continue;
        }
        if(pulse_map->count(key) == 0){
            std::cout << "oops! no pulses for " << key << std::endl;
            continue;
        }

        uint32_t channel = p.second;
        std::vector<CCMRecoPulse> pulse_series = pulse_map->at(key);

        // let's grab our baseline
        double baseline = 0;
        for(std::pair<CCMPMTKey const, BaselineEstimate> const & it : baselines){
            CCMPMTKey baseline_key = it.first;
            if(key == baseline_key){
                baseline = it.second.baseline;
            }
        }

        if(baseline == 0){
            std::cout << "oops! no baseine for key " << key << std::endl;
            return;
        }

        CCMWaveformUInt16 const & waveform = waveform_map->at(channel);
        std::vector<short unsigned int> const & samples = waveform.GetWaveform();

        CCMWaveformDouble const & electronics_corrected_wf = corrected_waveform_map->at(channel);
        std::vector<double> corrected_wf = electronics_corrected_wf.GetWaveform();

        // let's invert and subtract the baseline off our wf
        std::vector<double> inv_wf_minus_baseline(samples.size());
        for(size_t i = 0; i < samples.size(); ++i){
            inv_wf_minus_baseline[i] = -1 * ((double)samples[i] + baseline);
        }


        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration_.pmtCal.find(key);

		// first let's get our refolded wf
        std::vector<double> refolded_wf(waveform.GetWaveform().size(), 0.0);
        GetRefoldedWf(refolded_wf, pulse_series, calib->second);

        // now let's get our chi^2
        double chi = 0;
        double dof = 0;
		bool borked = false;


        GetChi(chi, dof, refolded_wf, corrected_wf);
        // std::cout << "for " << key << " chi/dof = " << chi/dof << std::endl;
        //GetChi(chi, refolded_wf, inv_wf_minus_baseline);
        // let's see if this is a big chi^2
        if (chi > chi_threshold_){
            borked = true;
        }

        // now let's save everything
        chi_map->emplace(key, chi);
        dof_map->emplace(key, dof);

        refolded_wf_map->emplace(key, refolded_wf);
    }

	frame->Put(chi_name_, chi_map);
	frame->Put("DegreesofFreedom", dof_map);
	frame->Put("RefoldedWaveforms", refolded_wf_map);

	PushFrame(frame);
}
