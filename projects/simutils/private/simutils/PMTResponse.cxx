#include <icetray/IcetrayFwd.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <set>
#include <tuple>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <filesystem>

#include <icetray/open.h>
#include <icetray/I3Frame.h>
#include <icetray/I3TrayInfo.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/I3DefaultName.h>

#include <dataclasses/I3Double.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/calibration/CCMSimulationCalibration.h>

#include <simclasses/CCMMCPE.h>
#include <simclasses/SampledRecoPulse.h>
#include <simclasses/DetectorResponseConfig.h>

#include <phys-services/I3GSLRandomService.h>

struct DiscreteDistribution {
    DiscreteDistribution() = default;
    std::vector<double> probs_;
    DiscreteDistribution(std::vector<double> const & probs) : probs_(probs) {
        double sum = 0.0;
        for(double const & p : probs_)
            sum += p;
        std::transform(probs_.begin(), probs_.end(), probs_.begin(), [=](double p) { return p / sum; });
        for(size_t i = 1; i < probs_.size(); ++i)
            probs_[i] += probs_[i - 1];

    }
    size_t operator()(I3RandomServicePtr rng) {
        double r = rng->Uniform(0.0, 1.0);
        size_t i = std::distance(probs_.begin(), std::lower_bound(probs_.begin(), probs_.end(), r));
        return i;
    }
};

class PMTResponse: public I3Module {
    double ratio_singlet_;
    double ratio_triplet_;
    double singlet_time_constant_;
    double triplet_time_constant_;
    double intermediate_time_constant_;
    double uv_absorption_a_;
    double uv_absorption_b_;
    double uv_absorption_d_;
    double uv_absorption_scaling_;

    double uv_absorption_min_wavelength_;
    double uv_absorption_max_wavelength_ = 200.0; // nm
    double uv_absorption_min_abs_length_;
    double uv_absorption_max_abs_length_;

    CCMSimulationPMTCalibrationMap pmt_cal_;
    boost::shared_ptr<const I3MapPMTKeyDouble> spe_lower_threshold_;
    boost::shared_ptr<const I3MapPMTKeyDouble> pmt_transit_times_;
    I3Map<CCMPMTKey, CCMTriggerKey> trigger_copy_map_;

    std::string input_hits_map_name_;
    std::string output_reco_pulse_name_;
    std::string output_time_offsets_name_;
    std::string detector_configuration_name_;
    std::string randomServiceName_;
    std::string output_true_event_time_name_;

    I3RandomServicePtr randomService_;

    double photon_sampling_factor_;
    bool simulated_enable_uv_absorption_;
    double simulated_uv_absorption_a_;
    double simulated_uv_absorption_b_;
    double simulated_uv_absorption_d_;
    double simulated_uv_absorption_scaling_;
    double simulated_uv_absorption_min_wavelength_;
    double simulated_uv_absorption_max_wavelength_ = 200.0; // nm
    double simulated_uv_absorption_min_abs_length_;
    double simulated_uv_absorption_max_abs_length_;

    double average_pmt_eff_;
    bool flat_eff_;
    bool remove_cherenkov_;
    bool weight_uv_absorption_;
    bool track_sampled_pe_;
    bool seen_cal_frame = false;
    bool check_for_detailed_photon_tracking = false;
    I3FramePtr calibration_frame_ = nullptr;
    std::vector<double> wavelength_qe_wavelength;
    std::vector<double> wavelength_qe_efficiency;

    CCMPulseTimeDistributionParameters default_main_pulse_ = {-0.45, 0.9, 0.7971466712310554};
    std::vector<CCMPulseTimeDistributionParameters> default_late_pulses_ = {
        {30.946544, 7.737966, 0.041169},
        {47.129206, 3.287489, 0.015341},
        {384.661399, 166.957997, 0.146343},
    };
    bool generate_late_pulses_ = true;

    std::map<CCMPMTKey, std::vector<CCMPulseTimeDistributionParameters>> pulse_types_;
    std::map<CCMPMTKey, DiscreteDistribution> pulse_selector_;

public:
    PMTResponse(const I3Context&);
    void Configure();
    void Geometry(I3FramePtr frame);
    void Simulation(I3FramePtr frame);
    void Calibration(I3FramePtr frame);
    double NormalDistributionCDF(double mu, double sigma, double x);
    double NormalDistributionInverseCDF(double p, double mu, double sigma);
    double ExponentialInverseCDF(double p, double tau);
    double IntermediateInverseCDF(double p, double tau);
    void DAQ(I3FramePtr frame);
    void Finish();

    static const std::vector<double> default_wavelength_qe_wavelength;
    static const std::vector<double> default_wavelength_qe_efficiency;

    static double UVAbsorptionLengthCM(double wavelength, double a, double b, double d, double scaling, double min_wavelength, double max_wavelength, double min_abs_length, double max_abs_length);
    static double UVAbsorptionScalingCM(double wavelength, double distance, double a, double b, double d, double scaling, double min_wavelength, double max_wavelength, double min_abs_length, double max_abs_length);
    double SimulatedUVAbsorptionCM(double wavelength) const;
    double SimulatedUVAbsorptionScalingCM(double wavelength, double distance) const;
    double ResponseUVAbsorptionCM(double wavelength) const;
    double ResponseUVAbsorptionScalingCM(double wavelength, double distance) const;

    double SampleGEV(double mu, double sigma);
    double SampleGEV(CCMPulseTimeDistributionParameters const & pulse) {
        return SampleGEV(pulse.mu, pulse.sigma);
    }

    int CheckForDetailedPhotonTracking(I3FramePtr frame) const;
    void UpdateCalibrationFrame();
};

int PMTResponse::CheckForDetailedPhotonTracking(I3FramePtr frame) const {
    boost::shared_ptr<CCMMCPESeriesMap const> mcpeseries_source = frame->Get<boost::shared_ptr<CCMMCPESeriesMap const>>(input_hits_map_name_);
    bool has_photon = false;
    for(CCMMCPESeriesMap::const_iterator it = mcpeseries_source->begin(); it != mcpeseries_source->end(); ++it)
        for(CCMMCPESeries::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            CCMMCPE const & pe = *it2;
            if(pe.distance_uv > 0.0 or pe.distance_visible > 0.0)
                return true;
            else
                has_photon = true;
        }
    if(not has_photon)
        return -1;
    return false;
}

void PMTResponse::UpdateCalibrationFrame() {
    if(weight_uv_absorption_)
        return;
    std::string sim_cal_name_ = "SimulationCalibration";
    CCMSimulationCalibration cal = calibration_frame_->Get<CCMSimulationCalibration>(sim_cal_name_);
    boost::shared_ptr<CCMSimulationCalibration> sim_calibration = boost::make_shared<CCMSimulationCalibration>(cal);
    calibration_frame_->Delete(sim_cal_name_);

    sim_calibration->uv_absorption_a = simulated_uv_absorption_a_; // [1/nm]
    sim_calibration->uv_absorption_b = simulated_uv_absorption_b_; // [nm]
    sim_calibration->uv_absorption_d = simulated_uv_absorption_d_; // [cm]
    sim_calibration->uv_absorption_scaling = simulated_uv_absorption_scaling_;

    calibration_frame_->Put(sim_cal_name_, sim_calibration);
}

double PMTResponse::UVAbsorptionLengthCM(double wavelength, double a, double b, double d, double scaling, double min_wavelength, double max_wavelength, double min_abs_length, double max_abs_length) {
    if(wavelength < min_wavelength) {
        return min_abs_length * I3Units::cm;
    } else if(wavelength > max_wavelength) {
        return max_abs_length * I3Units::cm;
    } else {
        double function_T = 1.0 - std::exp( - a * (wavelength - b));
        return (d / std::log(1.0 / function_T)) * scaling * I3Units::cm;
    }
}

double PMTResponse::UVAbsorptionScalingCM(double wavelength, double distance, double a, double b, double d, double scaling, double min_wavelength, double max_wavelength, double min_abs_length, double max_abs_length) {
    if(wavelength < min_wavelength) {
        return std::exp(-distance / (min_abs_length * I3Units::cm));
    } else if(wavelength > max_wavelength) {
        return std::exp(-distance / (max_abs_length * I3Units::cm));
    } else {
        return std::pow(1.0 - std::exp(a * (b - wavelength)), distance / (scaling * d * I3Units::cm));
    }
}

double PMTResponse::SimulatedUVAbsorptionCM(double wavelength) const {
    return UVAbsorptionLengthCM(wavelength, simulated_uv_absorption_a_, simulated_uv_absorption_b_, simulated_uv_absorption_d_, simulated_uv_absorption_scaling_, simulated_uv_absorption_min_wavelength_, simulated_uv_absorption_max_wavelength_, simulated_uv_absorption_min_abs_length_, simulated_uv_absorption_max_abs_length_);
}

double PMTResponse::SimulatedUVAbsorptionScalingCM(double wavelength, double distance) const {
    return UVAbsorptionScalingCM(wavelength, distance, simulated_uv_absorption_a_, simulated_uv_absorption_b_, simulated_uv_absorption_d_, simulated_uv_absorption_scaling_, simulated_uv_absorption_min_wavelength_, simulated_uv_absorption_max_wavelength_, simulated_uv_absorption_min_abs_length_, simulated_uv_absorption_max_abs_length_);
}

double PMTResponse::ResponseUVAbsorptionCM(double wavelength) const {
    return UVAbsorptionLengthCM(wavelength, uv_absorption_a_, uv_absorption_b_, uv_absorption_d_, uv_absorption_scaling_, uv_absorption_min_wavelength_, uv_absorption_max_wavelength_, uv_absorption_min_abs_length_, uv_absorption_max_abs_length_);
}

double PMTResponse::ResponseUVAbsorptionScalingCM(double wavelength, double distance) const {
    return UVAbsorptionScalingCM(wavelength, distance, uv_absorption_a_, uv_absorption_b_, uv_absorption_d_, uv_absorption_scaling_, uv_absorption_min_wavelength_, uv_absorption_max_wavelength_, uv_absorption_min_abs_length_, uv_absorption_max_abs_length_);
}

double PMTResponse::SampleGEV(double mu, double sigma) {
    double u = randomService_->Uniform(0.0, 1.0);
    return mu - sigma * std::log(-std::log(u));
}

const std::vector<double> PMTResponse::default_wavelength_qe_wavelength = {
    271.7297065293589, 272.5905302357397, 274.7689165854319, 275.774504708572, 277.42549012103206,
    279.60684224611464, 281.3715166369379, 284.02104543387406, 286.31476684771184, 288.36924963284537,
    291.20241580811705, 294.7489082161423, 297.3472537419332, 301.7761955951739, 308.1508675168173,
    315.0744829262529, 328.8322696687711, 348.0242036526301, 369.01527273218346, 389.99718683974567,
    410.969762875877, 431.01860223810456, 451.05865282722067, 471.091196339304, 488.38255813780705,
    508.41400305325146, 525.6989563713607, 537.484151815526, 548.3526497420297, 560.5942673320314,
    575.1251487359968, 589.6524047710537, 602.3532805137708, 614.1377435601768, 624.5502953993878,
    634.0448456460379, 643.0291587870497, 650.2913321812506, 657.0493219407376, 663.3466380871176,
    669.6431989483084, 675.4830136793763, 680.9333817872812, 685.7977420998144, 691.1783937881538,
    695.7505266681454, 700.5754003512864, 704.6935508854773, 709.4535259893923, 714.0548672262536,
    718.4585395392952, 722.9341917297394
};

const std::vector<double> PMTResponse::default_wavelength_qe_efficiency = {
    0.002398740835564059, 0.003121371056039485, 0.00425850887297264, 0.006001329657904212, 0.008024056054700794,
    0.011056341101904597, 0.015094429573531078, 0.020696143741601142, 0.02849836945068553, 0.03826555230720272,
    0.052991415285607905, 0.07240331867023701, 0.09610228291102357, 0.13207165321064374, 0.1815343556927802,
    0.23776843161775388, 0.30248367358133543, 0.3364736243973209, 0.34354861678310683, 0.3402008658768985,
    0.32653283426099217, 0.30174406822319666, 0.27076483199514884, 0.23694524445692558, 0.20457669204936288,
    0.17836806326025914, 0.15073783540716676, 0.11699878242414168, 0.08951160789158077, 0.06952998759784439,
    0.0555286953033716, 0.04381269414550079, 0.0343829721767211, 0.026621914607480143, 0.020384246137586923,
    0.015317782194202512, 0.011382145604947418, 0.00849320551185922, 0.006394900074833947, 0.0047434076893910345,
    0.0035095442415119683, 0.002591832591784181, 0.0019241990489623855, 0.001431051684337872, 0.001046402936813468,
    0.0007715369207947706, 0.0005675436529287897, 0.0004215655226197218, 0.0003103859991646788, 0.0002290412779971303,
    0.00016780030982719966, 0.00012576417985507728
};

I3_MODULE(PMTResponse);

PMTResponse::PMTResponse(const I3Context& context) : I3Module(context),
    input_hits_map_name_("PMTMCHitsMap"), output_reco_pulse_name_("MCRecoPulses"),
    output_time_offsets_name_("SimulatedBoardTimeOffsets"), detector_configuration_name_("DetectorResponseConfig"),
    output_true_event_time_name_("TrueEventTime"),
    flat_eff_(false), remove_cherenkov_(false), weight_uv_absorption_(false) {
    AddParameter("InputHitsMapName", "Name of the input hits map", input_hits_map_name_);
    AddParameter("OutputRecoPulseName", "Name of the output reco pulse series map", output_reco_pulse_name_);
    AddParameter("OutputTimeOffsetsName", "Name of the output time offsets map", output_time_offsets_name_);
    AddParameter("OutputTrueEventTimeName", "Name of the output true event time", output_true_event_time_name_);
    AddParameter("DetectorConfigurationName", "Name of the detector configuration in the Simulation frame", detector_configuration_name_);
    AddParameter("RandomServiceName", "Name of the random service in the context. If empty default random service will be used.", randomServiceName_);
    AddParameter("QEWavelengths", "Wavelengths for quantum efficiency", wavelength_qe_wavelength);
    AddParameter("QEValues", "Values for quantum efficiency", wavelength_qe_efficiency);
    AddParameter("RemoveCherenkov", "true removes cherenkov photons, false keeps them all", remove_cherenkov_);
    AddParameter("FlatEfficiency", "true to use flat efficiency", flat_eff_);
    AddParameter("WeightUVAbsorption", "true to add uv absorption via reweighting", weight_uv_absorption_);
}


void PMTResponse::Configure() {
    GetParameter("InputHitsMapName", input_hits_map_name_);
    GetParameter("OutputRecoPulseName", output_reco_pulse_name_);
    GetParameter("OutputTimeOffsetsName", output_time_offsets_name_);
    GetParameter("OutputTrueEventTimeName", output_true_event_time_name_);
    GetParameter("DetectorConfigurationName", detector_configuration_name_);
    GetParameter("RandomServiceName", randomServiceName_);
    if(randomServiceName_.empty()) {
        randomService_ = I3RandomServicePtr(new I3GSLRandomService(0));
        log_debug("+ Random service: I3GSLRandomService  (default)");
    } else {
        randomService_ = GetContext().Get<I3RandomServicePtr>(randomServiceName_);
        if(randomService_) log_debug("+ Random service: %s  (EXTERNAL)",  randomServiceName_.c_str());
        else log_fatal("No random service \"%s\" in context!", randomServiceName_.c_str());
    }
    GetParameter("QEWavelengths", wavelength_qe_wavelength);
    GetParameter("QEValues", wavelength_qe_efficiency);
    if(wavelength_qe_wavelength.size() != wavelength_qe_efficiency.size()) {
        log_fatal("Wavelengths and QE values must have the same size!");
    }
    if(wavelength_qe_wavelength.empty()) {
        wavelength_qe_wavelength = default_wavelength_qe_wavelength;
        wavelength_qe_efficiency = default_wavelength_qe_efficiency;
        log_debug("Using default wavelength qe values");
    }
    GetParameter("RemoveCherenkov", remove_cherenkov_);
    GetParameter("FlatEfficiency", flat_eff_);
    GetParameter("WeightUVAbsorption", weight_uv_absorption_);
}

void PMTResponse::Simulation(I3FramePtr frame) {
    if(not frame->Has(detector_configuration_name_)) {
        log_fatal("No detector configuration found in frame with name %s", detector_configuration_name_.c_str());
    }

    DetectorResponseConfig const & detector_config = frame->Get<DetectorResponseConfig>(detector_configuration_name_);
    photon_sampling_factor_ = detector_config.photon_sampling_factor_;
    simulated_enable_uv_absorption_ = detector_config.enable_uv_absorption_;
    simulated_uv_absorption_a_ = detector_config.uv_absorption_a_; // [1/nm]
    simulated_uv_absorption_b_ = detector_config.uv_absorption_b_; // [nm]
    simulated_uv_absorption_d_ = detector_config.uv_absorption_d_; // [cm]
    simulated_uv_absorption_scaling_ = detector_config.uv_absorption_scaling_;

    simulated_uv_absorption_min_wavelength_ = simulated_uv_absorption_b_ + 0.1;
    simulated_uv_absorption_max_wavelength_ = 200.0;
    double min_wavelength_function_T = 1.0 - std::exp( - simulated_uv_absorption_a_ * (simulated_uv_absorption_min_wavelength_ - simulated_uv_absorption_b_));
    simulated_uv_absorption_min_abs_length_ = (simulated_uv_absorption_d_ / std::log(1.0 / min_wavelength_function_T)); // units of cm!
    simulated_uv_absorption_min_abs_length_ *= simulated_uv_absorption_scaling_;
    double max_wavelength_function_T = 1.0 - std::exp( - simulated_uv_absorption_a_ * (simulated_uv_absorption_max_wavelength_ - simulated_uv_absorption_b_));
    simulated_uv_absorption_max_abs_length_ = (simulated_uv_absorption_d_ / std::log(1.0 / max_wavelength_function_T)); // units of cm!
    simulated_uv_absorption_max_abs_length_ *= simulated_uv_absorption_scaling_;

    PushFrame(frame);
}

void PMTResponse::Calibration(I3FramePtr frame) {
    // grab our simulation calibration out of the frame
    boost::shared_ptr<CCMSimulationCalibration const> sim_calibration = frame->Get<boost::shared_ptr<CCMSimulationCalibration const>>("SimulationCalibration");

    // and set our parameters using calibration frame
    ratio_singlet_ = sim_calibration->Rs;
    ratio_triplet_ = sim_calibration->Rt;
    singlet_time_constant_ = sim_calibration->tau_s;
    triplet_time_constant_ = sim_calibration->tau_t;
    intermediate_time_constant_ = sim_calibration->tau_other;
    pmt_cal_ = sim_calibration->pmt_calibration;
    uv_absorption_a_ = sim_calibration->uv_absorption_a; // [1/nm]
    uv_absorption_b_ = sim_calibration->uv_absorption_b; // [nm]
    uv_absorption_d_ = sim_calibration->uv_absorption_d; // [cm]
    uv_absorption_scaling_ = sim_calibration->uv_absorption_scaling;

    uv_absorption_min_wavelength_ = uv_absorption_b_ + 0.1;
    uv_absorption_max_wavelength_ = 200.0;
    double min_wavelength_function_T = 1.0 - std::exp( - uv_absorption_a_ * (uv_absorption_min_wavelength_ - uv_absorption_b_));
    uv_absorption_min_abs_length_ = (uv_absorption_d_ / std::log(1.0 / min_wavelength_function_T));
    uv_absorption_min_abs_length_ *= uv_absorption_scaling_;
    double max_wavelength_function_T = 1.0 - std::exp( - uv_absorption_a_ * (uv_absorption_max_wavelength_ - uv_absorption_b_));
    uv_absorption_max_abs_length_ = (uv_absorption_d_ / std::log(1.0 / max_wavelength_function_T));
    uv_absorption_max_abs_length_ *= uv_absorption_scaling_;

    spe_lower_threshold_ = frame->Get<boost::shared_ptr<I3MapPMTKeyDouble const>>("SPELowerThreshold");

    seen_cal_frame = true;

    // before we're done, let's grab the average pmt eff
    average_pmt_eff_ = 0.0;
    size_t total_keys = 0;
    for(const std::pair<CCMPMTKey const, CCMSimulationPMTCalibration>& entry : pmt_cal_) {
        average_pmt_eff_ += entry.second.pmt_efficiency;
        total_keys += 1;
    }

    average_pmt_eff_ /= (static_cast<double>(total_keys));

    pulse_types_.clear();
    for(std::pair<CCMPMTKey const, CCMTriggerKey> const & key : trigger_copy_map_) {
        auto it = pmt_cal_.find(key.first);
        std::vector<double> probs;
        std::vector<CCMPulseTimeDistributionParameters> pulses;
        if(generate_late_pulses_) {
            if(it != pmt_cal_.end()) {
                CCMSimulationPMTCalibration const & cal = it->second;
                std::vector<CCMPulseTimeDistributionParameters> const & late_pulse_params = cal.late_pulses;
                double sum = 0.0;
                for(CCMPulseTimeDistributionParameters const & pulse : late_pulse_params)
                    sum += pulse.fraction;
                double main_pulse_scale = 1.0 - sum;
                assert(main_pulse_scale >= 0.0);

                probs.push_back(main_pulse_scale);
                CCMPulseTimeDistributionParameters main_pulse = {cal.main_pulse_mu, cal.main_pulse_sigma, main_pulse_scale};
                pulses.emplace_back(main_pulse);

                for(CCMPulseTimeDistributionParameters const & pulse : late_pulse_params) {
                    probs.push_back(pulse.fraction);
                    pulses.emplace_back(pulse);
                }
            } else {
                // no late pulse params for this pmt, use default
                double sum = 0.0;
                for(CCMPulseTimeDistributionParameters const & pulse : default_late_pulses_)
                    sum += pulse.fraction;
                double main_pulse_scale = 1.0 - sum;
                assert(main_pulse_scale >= 0.0);

                probs.push_back(main_pulse_scale);
                CCMPulseTimeDistributionParameters main_pulse = {default_main_pulse_.mu, default_main_pulse_.sigma, main_pulse_scale};
                pulses.emplace_back(main_pulse);

                for(CCMPulseTimeDistributionParameters const & pulse : default_late_pulses_) {
                    probs.push_back(pulse.fraction);
                    pulses.emplace_back(pulse);
                }
            }
        } else {
            probs.push_back(1.0);
            if(it != pmt_cal_.end()) {
                CCMSimulationPMTCalibration const & cal = it->second;
                CCMPulseTimeDistributionParameters main_pulse = {cal.main_pulse_mu, cal.main_pulse_sigma, 1.0};
                pulses.emplace_back(main_pulse);
            } else {
                CCMPulseTimeDistributionParameters main_pulse = {default_main_pulse_.mu, default_main_pulse_.sigma, 1.0};
                pulses.emplace_back(main_pulse);
            }
        }
        pulse_types_[key.first] = pulses;
        pulse_selector_[key.first] = DiscreteDistribution(probs);
    }

    if(weight_uv_absorption_) {
        check_for_detailed_photon_tracking = true;
        PushFrame(frame);
    } else {
        calibration_frame_ = frame;
    }
}

double PMTResponse::NormalDistributionCDF(double mu, double sigma, double x) {
    return 0.5 * (1.0 + std::erf((x - mu) / (sigma * std::sqrt(2.0))));
}

double PMTResponse::NormalDistributionInverseCDF(double p, double mu, double sigma) {
    if(p <= 0.0 or p >= 1.0) {
        return std::numeric_limits<double>::infinity();
    } else {
        return mu + sigma * sqrt(2.0) * boost::math::erf_inv(2 * p - 1.0);
    }
}

double PMTResponse::ExponentialInverseCDF(double p, double tau) {
    return - tau * std::log(1.0 - p);
}

double PMTResponse::IntermediateInverseCDF(double p, double tau) {
    return tau * ((1.0 / (1.0 - p)) - 1.0);
}

void PMTResponse::Geometry(I3FramePtr frame) {
    // ok from our geometry file we need to grab the trigger copy map
    // this map is pmt key : board
    // and we need to add random offsets for each board
    CCMGeometry const & geo = frame->Get<CCMGeometry const>("CCMGeometry");
    trigger_copy_map_ = geo.trigger_copy_map;

    // grab our transit times
    pmt_transit_times_ = frame->Get<boost::shared_ptr<I3MapPMTKeyDouble const>>("PMTTransitTimes");

    PushFrame(frame);
}

void PMTResponse::DAQ(I3FramePtr frame) {
    // check to make sure we've seen a calibration frames
    if(not seen_cal_frame) {
        log_fatal("Did not see calibration frame -- please provide calibration frame first so I can reweight");
    }

    if(check_for_detailed_photon_tracking) {
        int has_detailed_photon_tracking = CheckForDetailedPhotonTracking(frame);
        if(has_detailed_photon_tracking < 0) {
            // No photons in this event
            // Cannot use it to check for detailed photon weighting
            // So we just skip it
        } else if(has_detailed_photon_tracking) {
            // We have detailed photon tracking
            // Don't need to check for it again
            check_for_detailed_photon_tracking = false;
        } else {
            log_fatal("Simulation does not have detailed photon tracking, but it is needed to reweight for UV absorption. Either resimulate with DetailedPhotonTracking=True, or disable UV absorption reweighting by setting WeightUVAbsorption=False.");
        }
    }

    // Update the calibration frame if necessary
    // This is done when UV absorption is handled by Geant4 instead of in this module
    // We ideally want the SimulationCalibration object to reflect the settings used by Geant4 if we are not overriding them
    if(calibration_frame_) {
        UpdateCalibrationFrame();
        PushFrame(calibration_frame_);
        calibration_frame_ = nullptr;
    }

    // before we start, for this event, let's take our trigger copy map and assign some time smearing based on boards

    // ok so this frame should contain a PMTMCHitsMap
    // we want to read that in, apply various efficiencies
    I3MapPMTKeyDouble this_event_board_time_offset;
    I3MapPMTKeyDouble this_event_board_time_error;
    I3Map<CCMTriggerKey, double> board_offsets;
    I3Map<CCMTriggerKey, double> board_errors;
    I3MapPMTKeyDoublePtr this_event_total_time_offsets = boost::make_shared<I3MapPMTKeyDouble>();

    for(std::pair<CCMPMTKey const, CCMTriggerKey> const & key : trigger_copy_map_) {
        if(board_offsets.find(key.second) != board_offsets.end()) {
            // ok we already have an offset time for this board, let's just save
            this_event_board_time_offset.insert(std::make_pair(key.first, board_offsets[key.second]));
            this_event_board_time_error.insert(std::make_pair(key.first, board_errors[key.second]));
        } else {
            // first time seeing this board! let's make an offset
            double b_offset = randomService_->Uniform(-1.0, 1.0);
            double b_error = randomService_->Uniform(-0.1, 0.1);
            board_offsets.insert(std::make_pair(key.second, b_offset));
            this_event_board_time_offset.insert(std::make_pair(key.first, board_offsets[key.second]));
            board_errors.insert(std::make_pair(key.second, b_error));
            this_event_board_time_error.insert(std::make_pair(key.first, board_errors[key.second]));
        }
    }

    std::map<CCMTriggerKey, std::tuple<std::vector<CCMPMTKey>, double, double>> _offsets;

    for(std::pair<CCMPMTKey const, CCMTriggerKey> const & key : trigger_copy_map_) {
        if(_offsets.find(key.second) != _offsets.end()) {
            std::get<0>(_offsets[key.second]).push_back(key.first);
        } else {
            std::vector<CCMPMTKey> pmts;
            pmts.push_back(key.first);
            _offsets.insert(std::make_pair(key.second, std::make_tuple(pmts, this_event_board_time_offset[key.first], this_event_board_time_error[key.first])));
        }
    }

    // set up object to hold re-weighed simulation
    CCMRecoPulseSeriesMapPtr mcpeseries_dest = boost::make_shared<CCMRecoPulseSeriesMap>();

    // grab the simulation output ccm mcpe series map
    boost::shared_ptr<CCMMCPESeriesMap const> mcpeseries_source = frame->Get<boost::shared_ptr<CCMMCPESeriesMap const>>(input_hits_map_name_);

    // grab a random number for the overall event time offset
    double event_time_offset = randomService_->Uniform(-10.0, 10.0);

    // Iterate over PMTs in source map
    for(CCMMCPESeriesMap::const_iterator it = mcpeseries_source->begin(); it != mcpeseries_source->end(); ++it) {
        std::stringstream ss;
        ss << it->first;

        // Find the corresponding PMT in the destination map
        CCMRecoPulseSeriesMap::iterator it_dest = mcpeseries_dest->find(it->first);

        // If the PMT is not in the destination, then insert an empty vector
        if(it_dest == mcpeseries_dest->end()) {
            mcpeseries_dest->insert(std::make_pair(it->first, CCMRecoPulseSeries()));
            // Update the iterator so it points to our new entry
            it_dest = mcpeseries_dest->find(it->first);
        }

        // Reference to the destination
        CCMRecoPulseSeries & dest_series = it_dest->second;

        CCMSimulationPMTCalibration const & pmt_cal = pmt_cal_.at(it->first);

        // While we are looping over tubes, grab the PMT efficiency
        double pmt_efficiency = pmt_cal.pmt_efficiency;
        if(flat_eff_) {
            pmt_efficiency = average_pmt_eff_;
        }

        double this_tube_board_time_offset = this_event_board_time_offset.at(it->first);
        double this_tube_board_time_error = this_event_board_time_error.at(it->first);

        double this_tube_transit_time = 0.0;
        I3MapPMTKeyDouble::const_iterator pmt_tt_it = pmt_transit_times_->find(it->first);
        if(pmt_tt_it != pmt_transit_times_->end()) {
            this_tube_transit_time = pmt_tt_it->second;
            this_event_total_time_offsets->insert(std::make_pair(it->first, this_tube_transit_time + this_tube_board_time_offset));
        }

        double this_tube_spe_threshold = 0.0;
        I3MapPMTKeyDouble::const_iterator spe_threshold_it = spe_lower_threshold_->find(it->first);
        if(spe_threshold_it != spe_lower_threshold_->end()) {
            this_tube_spe_threshold = spe_threshold_it->second;
        }

        double this_tube_spe_threshold_efficiency = 1.0 - NormalDistributionCDF(pmt_cal.pmt_spe_mu, pmt_cal.pmt_spe_sigma, this_tube_spe_threshold);

        DiscreteDistribution & pulse_selector = pulse_selector_.at(it->first);
        std::vector<CCMPulseTimeDistributionParameters> const & pulse_types = pulse_types_.at(it->first);

        // Iterate over the vector of CCMMCPE in the source map for this PMT
        std::vector<CCMRecoPulse> temp_series; temp_series.reserve(size_t(it->second.size() * 0.5));
        for(CCMMCPESeries::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            CCMMCPE const & pe = *it2;

            // If we want to remove cherenkov photons, do so now
            if(remove_cherenkov_ and pe.photon_source == CCMMCPE::PhotonSource::Cherenkov) {
                continue;
            }

            // first, grab random number 0 - 1 to determine if this photon survives
            double survival = randomService_->Uniform(0.0, 1.0);

            // Check if we surive wavelength response of the tube
            double wavelength = pe.wavelength / I3Units::nanometer;

            // find the first element in wavelength_qe_wavelength that is greater than or equal to wavelength
            double wavelength_qe_weighting;
            std::vector<double>::iterator lower_it = std::lower_bound(wavelength_qe_wavelength.begin(), wavelength_qe_wavelength.end(), wavelength);
            if(lower_it == wavelength_qe_wavelength.begin() or lower_it == wavelength_qe_wavelength.end()) {
                wavelength_qe_weighting = 0.0;
            } else {
                size_t lower_index = std::distance(wavelength_qe_wavelength.begin(), lower_it);
                double wl_below = *lower_it;
                double wlqe_below = wavelength_qe_efficiency.at(lower_index);
                if(wl_below > wavelength) {
                    lower_index -= 1;
                    wl_below = wavelength_qe_wavelength.at(lower_index);
                    wlqe_below = wavelength_qe_efficiency.at(lower_index);
                }
                double wl_above =  wavelength_qe_wavelength.at(lower_index + 1);
                double wlqe_above =  wavelength_qe_efficiency.at(lower_index + 1);
                wavelength_qe_weighting = wlqe_below + (wavelength - wl_below) * ((wlqe_above - wlqe_below) / (wl_above - wl_below));
            }

            // Check if survive uv absorption cuts
            double uv_abs_probability = 1.0;
            std::cout << "uv_absorption_a: " << uv_absorption_a_ << ", " << simulated_uv_absorption_a_ << std::endl;
            std::cout << "uv_absorption_b: " << uv_absorption_b_ << ", " << simulated_uv_absorption_b_ << std::endl;
            std::cout << "uv_absorption_d: " << uv_absorption_d_ << ", " << simulated_uv_absorption_d_ << std::endl;
            std::cout << "uv_absorption_scaling: " << uv_absorption_scaling_ << ", " << simulated_uv_absorption_scaling_ << std::endl;
            if(weight_uv_absorption_) {
                double original_wavelength = pe.original_wavelength / I3Units::nanometer;
                double wavelength = pe.wavelength / I3Units::nanometer;
                double distance_travelled_before_wls = pe.distance_uv;
                double distance_travelled_after_wls = pe.distance_visible;

                double scaling_before = ResponseUVAbsorptionScalingCM(original_wavelength, distance_travelled_before_wls);
                double scaling_after = ResponseUVAbsorptionScalingCM(wavelength, distance_travelled_after_wls);
                std::cout << "absorption probability: " << scaling_before * scaling_after << std::endl;
                uv_abs_probability *= scaling_before * scaling_after;

                if(simulated_enable_uv_absorption_) {
                    double simulated_scaling_before = SimulatedUVAbsorptionScalingCM(original_wavelength, distance_travelled_before_wls);
                    double simulated_scaling_after = SimulatedUVAbsorptionScalingCM(wavelength, distance_travelled_after_wls);
                    std::cout << "simulated absorption probability: " << simulated_scaling_before * simulated_scaling_after << std::endl;
                    uv_abs_probability /= (simulated_scaling_before * simulated_scaling_after);
                }
            }
            std::cout << "uv_absorption probability: " << uv_abs_probability << std::endl;

            std::cout << "wavelength: " << wavelength << std::endl;
            std::cout << "wavelength_qe_weighting: " << wavelength_qe_weighting << std::endl;
            std::cout << "pmt_efficiency: " << pmt_efficiency << std::endl;
            std::cout << "this_tube_spe_threshold_efficiency: " << this_tube_spe_threshold_efficiency << std::endl;
            std::cout << "this_tube_spe_threshold: " << this_tube_spe_threshold << std::endl;
            std::cout << "photon_sampling_factor_: " << photon_sampling_factor_ << std::endl;

            double survival_probability = this_tube_spe_threshold_efficiency * pmt_efficiency * wavelength_qe_weighting * uv_abs_probability / photon_sampling_factor_;
            std::cout << "survival probability: " << survival_probability << std::endl;

            double charge_scale_factor = 1.0;
            if(survival_probability > 1.0) {
                charge_scale_factor = survival_probability;
            } else if(survival > survival_probability) {
                continue;
            }

            // ok our photon has survived! yay! let's apply our various time offsets
            // 1 -- overall event time offset
            // 2 -- time jitter due to electron transit times
            // 3 -- physical time offset (either due to scintillation + tpb or to late pulse)

            double total_time_offset = event_time_offset + (pe.time / I3Units::nanosecond);

            CCMPulseTimeDistributionParameters const & pmt_pulse_type = pulse_types.at(pulse_selector(randomService_));

            double pulse_time_offset = SampleGEV(pmt_pulse_type);

            total_time_offset += pulse_time_offset;

            double scint_time_offset = 0.0;
            if(pe.photon_source == CCMMCPE::PhotonSource::Scintillation) {
                // now let's get a random number to see if we're in singlet, triplet, or intermediate times
                double time_distribution = randomService_->Uniform(0.0, 1.0);
                if(time_distribution < (ratio_singlet_ + ratio_triplet_)) {
                    // ok singlet or triplet! throw another random number to figure it out
                    double t = randomService_->Uniform(0.0, ratio_singlet_ + ratio_triplet_);
                    if(t < ratio_singlet_) {
                        // in singlet!
                        double singlet_rand = randomService_->Uniform(0.0, 1.0);
                        scint_time_offset = ExponentialInverseCDF(singlet_rand, singlet_time_constant_);
                        total_time_offset += scint_time_offset;
                    } else {
                        // in triplet
                        double triplet_rand = randomService_->Uniform(0.0, 1.0);
                        scint_time_offset += ExponentialInverseCDF(triplet_rand, triplet_time_constant_);
                        total_time_offset += scint_time_offset;
                    }
                } else {
                    // ok we are in the intermediate time component
                    double intermediate_rand = randomService_->Uniform(0.0, 1.0);
                    scint_time_offset = IntermediateInverseCDF(intermediate_rand, intermediate_time_constant_);
                    total_time_offset += scint_time_offset;
                }
            }

            // now, finally, sample our spe shape for this tube to assign the appropiate magnitude to our pulse
            double pulse_amplitude_rand = randomService_->Uniform(1.0 - this_tube_spe_threshold_efficiency, 1.0);
            double pulse_amplitude = NormalDistributionInverseCDF(pulse_amplitude_rand, pmt_cal.pmt_spe_mu, pmt_cal.pmt_spe_sigma);

            // check to make sure we dont have any infinities and our pulse amplitude is appropriate
            if(std::isinf(pulse_amplitude) or std::isinf(total_time_offset)) {
                continue;
            }

            // If we undersimulated the number of photons, then the next best thing we can do is scale up the charge
            pulse_amplitude *= charge_scale_factor;

            // note -- we are hacking the pmt width to be cherenkov spe
            double width = 0.0;
            if(pe.photon_source == CCMMCPE::PhotonSource::Cherenkov) {
                width = pulse_amplitude;
            }

            // ok time to save as a reco pulse!
            // add transit time for this tube
            total_time_offset += this_tube_transit_time + this_tube_board_time_offset + this_tube_board_time_error;

            // bin
            double binned_reco_time = static_cast<int>(total_time_offset / 2.0) * 2.0;

            // and now subtract off our transit time with error
            binned_reco_time -= this_tube_transit_time + this_tube_board_time_offset;

            temp_series.emplace_back();
            CCMRecoPulse & pulse = temp_series.back();
            pulse.SetCharge(pulse_amplitude);
            pulse.SetTime(binned_reco_time);
            pulse.SetWidth(width);
        }

        if(temp_series.empty()) {
            continue;
        }

        // Sort and merge the pulses
        std::sort(temp_series.begin(), temp_series.end(), [](const CCMRecoPulse& a, const CCMRecoPulse& b) { return a.GetTime() < b.GetTime(); });
        dest_series.push_back(temp_series.front());
        for(CCMRecoPulseSeries::const_iterator it = temp_series.begin() + 1; it != temp_series.end(); ++it) {
            CCMRecoPulse & dest_pulse = dest_series.back();
            if(it->GetTime() == dest_pulse.GetTime()) {
                dest_pulse.SetCharge(dest_pulse.GetCharge() + it->GetCharge());
                dest_pulse.SetWidth(dest_pulse.GetWidth() + it->GetWidth());
            } else {
                dest_series.push_back(*it);
            }
        }
    }

    frame->Put(output_true_event_time_name_, boost::make_shared<I3Double>(event_time_offset));
    frame->Put(output_time_offsets_name_, this_event_total_time_offsets);
    frame->Put(output_reco_pulse_name_, mcpeseries_dest);
    PushFrame(frame);
}

void PMTResponse::Finish() {}


