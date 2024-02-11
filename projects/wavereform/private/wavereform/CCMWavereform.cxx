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

	chi_name_ = "RecoPulse";
	AddParameter("Output", "Where to store the chi^2 between waveform and pulses", chi_name_);

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
	GetParameter("Output", chi_name_);
	GetParameter("ChiThreshold", chi_threshold_);
	GetParameter("Flag", flag_name_);
}

void CCMFillFWHM(double& start, double& stop, const std::vector<double>& data, double spacing, double min) {

    // Find the maximum
    double max = DBL_MIN;
    for (unsigned i = 0; i < data.size(); ++i) {
        if (data[i] > max) {
            max = data[i];
        }
    }

    double hm = max*0.5;
    unsigned hmbin = 0;
    for (unsigned i = 0; i < data.size(); ++i) {
        if (data[i] > hm) {
            hmbin = i;
            break;
        }
    }

    // Be conservative: choose the bin before the first rise above hm
    if (hmbin > 0) {
        --hmbin;
    }

    start = min + hmbin*spacing;

    for (int i = data.size() - 1; i >= 0; --i) {
        if (data[i] > hm) {
            hmbin = i;
            break;
        }
    }

    if (hmbin < data.size() - 1) {
        ++hmbin;
    }

    stop = min + hmbin*spacing;
    if (stop == start) {
        stop += spacing;
    }
}

void FillTemplate(CCMWaveformTemplate& wfTemplate, const CCMPMTCalibration& calibration, double const & start_time, int const & template_bins, double const & template_bin_spacing_) {
    CCMSPETemplate channel_template = calibration.GetSPETemplate();
    wfTemplate.digitizer_template.resize(template_bins);
    double total = 0.0;
    double wf_max = DBL_MIN;
    for(int i = 0; i < template_bins; i++) {
        double value = channel_template.Evaluate(start_time + i*template_bin_spacing_); // hard coding bin spacing...use calib one day
        wfTemplate.digitizer_template[i] = value;
        total += value;
        wf_max = std::max(wf_max, value);
    }

    double tail_fraction = 0.005;
    double target_amount = tail_fraction * total;
    double height_cut = tail_fraction * wf_max;

    double running_total = 0.0;
    size_t begin_idx = 0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        running_total += wfTemplate.digitizer_template[i];
        //if(running_total >= target_amount) {
        if(wfTemplate.digitizer_template[i] >= height_cut) {
            begin_idx = i;
            break;
        }
    }
    if(begin_idx > 0)
        begin_idx -= 1;

    target_amount = tail_fraction * total;
    running_total = 0.0;
    size_t end_idx = wfTemplate.digitizer_template.size();
    for(int i=wfTemplate.digitizer_template.size()-1; i>=0; --i) {
        running_total += wfTemplate.digitizer_template[i];
        //if(running_total >= target_amount) {
        if(wfTemplate.digitizer_template[i] >= height_cut) {
            end_idx = i;
            break;
        }
    }
    if(end_idx < wfTemplate.digitizer_template.size() - 1)
        end_idx += 1;

    wfTemplate.digitizer_template.erase(wfTemplate.digitizer_template.begin(), wfTemplate.digitizer_template.begin() + begin_idx);
    wfTemplate.start_time = start_time + template_bin_spacing_ * begin_idx;

    end_idx -= begin_idx;

    wfTemplate.digitizer_template.erase(wfTemplate.digitizer_template.begin() + end_idx, wfTemplate.digitizer_template.end());
    wfTemplate.end_time = wfTemplate.start_time + wfTemplate.digitizer_template.size() * template_bin_spacing_;

    wfTemplate.pulse_width = wfTemplate.end_time - wfTemplate.start_time;

    total = 0.0;
    for(size_t i=0; i<wfTemplate.digitizer_template.size(); ++i) {
        total += wfTemplate.digitizer_template[i];
    }

    wfTemplate.total_mass = total;

    CCMFillFWHM(wfTemplate.digitizerStart,
            wfTemplate.digitizerStop,
            wfTemplate.digitizer_template,
            template_bin_spacing_, wfTemplate.start_time);

    wfTemplate.filled = true;
}

void CCMWavereform::Calibration(I3FramePtr frame){
    calibration = frame->Get<CCMCalibration>("CCMCalibration");

    double start_time = 2 * -10;
    double end_time = 2 * 60;

    template_bin_spacing_ = 2.0 / spes_per_bin_;
    double range = end_time - start_time;

    for(size_t i=0; i<pmt_keys_.size(); ++i) {
        // let's get our wf and what not
        CCMPMTKey key = pmt_keys_[i];
        uint32_t channel = pmt_channel_map_[key];

        if(calibration.pmtCal.count(key) == 0) {
            // std::cout << "oops! no calibration for " << key << std::endl;
            continue;
        }

        std::map<CCMPMTKey, CCMPMTCalibration>::const_iterator calib = calibration.pmtCal.find(key);

        double placeholder = 1.0;

        if(not template_.at(key).filled) {
            int template_bins = (int)ceil(range / template_bin_spacing_);
            FillTemplate(template_.at(key), calib->second, start_time, template_bins, template_bin_spacing_);
            std::cout << key << " Template(" << template_.at(key).start_time << ", " << template_.at(key).end_time << ") FWHM(" << template_.at(key).digitizerStart << ", " << template_.at(key).digitizerStop << ")" << std::endl;
        }
    }
	PushFrame(frame);
}

void CCMWavereform::Geometry(I3FramePtr frame) {
    if(not frame->Has(geometry_name_)) {
        log_fatal("Could not find CCMGeometry object with the key named \"%s\" in the Geometry frame.", geometry_name_.c_str());
    }
    CCMGeometry const & geo = frame->Get<CCMGeometry const>(geometry_name_);
    pmt_channel_map_ = geo.pmt_channel_map;
    template_.clear();
    for(auto pmt_channel_pair : pmt_channel_map_) {
        template_[pmt_channel_pair.first] = CCMWaveformTemplate();
        template_[pmt_channel_pair.first].filled = false;
    }

    I3Map<CCMPMTKey, CCMOMGeo> const & pmt_geo_ = geo.pmt_geo;
    std::set<CCMPMTKey> allowed_pmt_keys(allowed_pmt_keys_.begin(), allowed_pmt_keys_.end());
    bool filter_pmts = allowed_pmt_keys.size() > 0;
    for(std::pair<CCMPMTKey const, CCMOMGeo> const & p : pmt_geo_) {
        if(filter_pmts and allowed_pmt_keys.find(p.first) == allowed_pmt_keys.end()) {
            continue;
        }
        if(p.second.omtype == CCMOMGeo::OMType::CCM8inCoated or p.second.omtype == CCMOMGeo::OMType::CCM8inUncoated or p.second.omtype == CCMOMGeo::OMType::CCM1in) {
            pmt_keys_.push_back(p.first);
        }
    }

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

void CCMWavereform::GetRefoldedWf(std::vector<double> & refolded_wf, std::vector<CCMRecoPulse> const & pulses, CCMPMTKey pmt_key){
    // time for how far we need to calculate our pulse out to
    // should replace with paramters from calibration one day
    double pulse_time;
    double pulse_charge;

    CCMWaveformTemplate & temp = template_[pmt_key];

    double templ_bin_spacing_inv = 1./template_bin_spacing_;
    std::vector<double> const & pulse_templ = temp.digitizer_template;
    size_t n_bins = pulse_templ.size();

    size_t nspes = pulses.size();
    int i = 0;
    int first_spe = 0;
    double last_t = 0;
    for(int i=0; i<refolded_wf.size(); ++i) {
        double bin_time = (i + 1) * 2.0;

        if (bin_time < last_t)
            first_spe = 0;
        last_t = bin_time;

        // The earliest pulse influencing this bin is PULSE_WIDTH in the past.
        // The template is defined up to (but not including) PULSE_WIDTH.
        while (first_spe < nspes && pulses[first_spe].GetTime() + temp.end_time <= bin_time) {
            first_spe++;
        }
        if (first_spe == nspes) {
            break;
        }

        // The last pulse for this bin is 2 ns in the future
        for (int j = first_spe; j < nspes; j++) {
            int templ_bin = int(((bin_time - pulses[j].GetTime()) - temp.start_time)*templ_bin_spacing_inv);
            if (templ_bin < 0)
                break;
            if (templ_bin >= n_bins)
                continue;

            refolded_wf.at(i) += pulse_templ.at(templ_bin) * pulses[j].GetCharge();
        }
    }
}


void CCMWavereform::DAQ(I3FramePtr frame) {
    const CCMCalibration& calibration_ = frame->Get<CCMCalibration>("CCMCalibration");
    I3Map<CCMPMTKey, BaselineEstimate> const & baselines= frame->Get<I3Map<CCMPMTKey, BaselineEstimate> const>("BaselineEstimates");
	CCMRecoPulseSeriesMapConstPtr pulse_map = frame->Get<CCMRecoPulseSeriesMapConstPtr>(pulse_name_);
    boost::shared_ptr<const CCMWaveformDoubleSeries> waveform_map = frame->Get<boost::shared_ptr<const CCMWaveformDoubleSeries>>(waveform_name_);

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

        CCMWaveformDouble const & electronics_corrected_wf = waveform_map->at(channel);
        std::vector<double> const & corrected_wf = electronics_corrected_wf.GetWaveform();

		// first let's get our refolded wf
        std::vector<double> refolded_wf(corrected_wf.size(), 0.0);
        GetRefoldedWf(refolded_wf, pulse_series, key);

        // now let's get our chi^2
        double chi = 0;
        double dof = 0;
		bool borked = false;

        GetChi(chi, dof, refolded_wf, corrected_wf);
        // let's see if this is a big chi^2
        if (chi > chi_threshold_){
            borked = true;
        }

        // now let's save everything
        chi_map->emplace(key, chi);
        dof_map->emplace(key, dof);

        refolded_wf_map->emplace(key, refolded_wf);
    }

	frame->Put(chi_name_ + "Chi2", chi_map);
	frame->Put(chi_name_ + "DegreesofFreedom", dof_map);
	frame->Put(chi_name_ + "RefoldedWaveforms", refolded_wf_map);

	PushFrame(frame);
}
