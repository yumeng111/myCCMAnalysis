#include <dataclasses/calibration/CCMSimulationCalibration.h>
#include <ostream>
#include <icetray/serialization.h>

CCMSimulationCalibration::CCMSimulationCalibration(){}


/*
    void SetPMTCalibration(CCMPMTKey key, const CCMSimulationPMTCalibration& pmt_cal);
    CCMSimulationPMTCalibration const & GetPMTCalibration(CCMPMTKey key) const;
    //CCMSimulationPMTCalibration & GetPMTCalibration(CCMPMTKey key);

    void SetPMTCalibration(CCMSimulationPMTCalibrationMap const & pmt_cal);
    CCMSimulationPMTCalibrationMap const & GetPMTCalibration() const;
    //CCMSimulationPMTCalibrationMap & GetPMTCalibration();
*/

void CCMSimulationCalibration::SetPMTCalibration(CCMPMTKey key, const CCMSimulationPMTCalibration& pmt_cal) {
    pmt_calibration[key] = pmt_cal;
}

CCMSimulationPMTCalibration const & CCMSimulationCalibration::GetPMTCalibration(CCMPMTKey key) const {
    return pmt_calibration.at(key);
}

//CCMSimulationPMTCalibration & CCMSimulationCalibration::GetPMTCalibration(CCMPMTKey key) {
//    return pmt_calibration[key];
//    }

void CCMSimulationCalibration::SetPMTCalibration(CCMSimulationPMTCalibrationMap const & pmt_cal) {
    pmt_calibration = pmt_cal;
}

CCMSimulationPMTCalibrationMap const & CCMSimulationCalibration::GetPMTCalibration() const {
    return pmt_calibration;
}

//CCMSimulationPMTCalibrationMap & CCMSimulationCalibration::GetPMTCalibration() {
//    return pmt_calibration;
//    }

template <class Archive>
void CCMSimulationPMTCalibration::serialize(Archive& ar, unsigned version) {
    if (version > ccmsimulationpmtcalibration_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMSimulationPMTCalibration class.", version, ccmsimulationpmtcalibration_version_);
    ar & make_nvp("pmt_efficiency", pmt_efficiency);
    ar & make_nvp("pmt_spe_mu", pmt_spe_mu);
    ar & make_nvp("pmt_spe_sigma", pmt_spe_sigma);
    ar & make_nvp("pmt_spe_threshold", pmt_spe_threshold);
    ar & make_nvp("main_pulse_mu", main_pulse_mu);
    ar & make_nvp("main_pulse_sigma", main_pulse_sigma);
    ar & make_nvp("late_pulses", late_pulses);
}

template <class Archive>
void CCMPulseTimeDistributionParameters::serialize(Archive& ar, unsigned version) {
    if (version > ccmlatepulseparameters_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMPulseTimeDistributionParameters class.", version, ccmlatepulseparameters_version_);
    ar & make_nvp("mu", mu);
    ar & make_nvp("sigma", sigma);
    ar & make_nvp("fraction", fraction);
}

template <class Archive>
void CCMSimulationCalibration::save(Archive& ar, unsigned version) const {
    if (version > ccmsimulationcalibration_version_)
        log_fatal("Attempting to save version %u from file but running version %u of CCMSimulationCalibration class.", version, ccmsimulationcalibration_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
    ar & make_nvp("PMTCalibration",pmt_calibration);
    ar & make_nvp("Rs",Rs);
    ar & make_nvp("Rt",Rt);
    ar & make_nvp("tau_s",tau_s);
    ar & make_nvp("tau_t",tau_t);
    ar & make_nvp("tau_other",tau_other);
    ar & make_nvp("uv_absorption_a",uv_absorption_a);
    ar & make_nvp("uv_absorption_b",uv_absorption_b);
    ar & make_nvp("uv_absorption_d",uv_absorption_d);
    ar & make_nvp("uv_absorption_scaling",uv_absorption_scaling);
}

template <class Archive>
void CCMSimulationCalibration::load(Archive& ar, unsigned version) {
    if (version > ccmsimulationcalibration_version_)
        log_fatal("Attempting to read version %u from file but running version %u of CCMSimulationCalibration class.", version, ccmsimulationcalibration_version_);

    ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));

    if(version <= 2) {
        I3MapPMTKeyDouble PMTEfficiencies;
        I3MapPMTKeyDouble LatePulseMu;
        I3MapPMTKeyDouble LatePulseSigma;
        I3MapPMTKeyDouble LatePulseScale;
        I3MapPMTKeyDouble PMTSPEMu;
        I3MapPMTKeyDouble PMTSPESigma;
        ar & make_nvp("PMTEfficiencies",PMTEfficiencies);
        ar & make_nvp("LatePulseMu",LatePulseMu);
        ar & make_nvp("LatePulseSigma",LatePulseSigma);
        ar & make_nvp("LatePulseScale",LatePulseScale);
        ar & make_nvp("PMTSPEMu",PMTSPEMu);
        ar & make_nvp("PMTSPESigma",PMTSPESigma);
        for(auto it = PMTEfficiencies.begin(); it != PMTEfficiencies.end(); ++it) {
            CCMPMTKey key = it->first;
            double eff = it->second;
            double late_mu = LatePulseMu[key];
            double late_sigma = LatePulseSigma[key];
            double late_scale = LatePulseScale[key];
            double mu = PMTSPEMu[key];
            double sigma = PMTSPESigma[key];
            CCMSimulationPMTCalibration cal;
            if(pmt_calibration.find(key) != pmt_calibration.end()) {
                cal = pmt_calibration[key];
            } else {
                cal = CCMSimulationPMTCalibration();
            }
            cal.pmt_efficiency = eff;
            cal.pmt_spe_mu = mu;
            cal.pmt_spe_sigma = sigma;
            cal.pmt_spe_threshold = 0.0;
            cal.main_pulse_mu = -0.45;
            cal.main_pulse_sigma = 0.9;
            cal.late_pulses.push_back(CCMPulseTimeDistributionParameters(late_mu, late_sigma, late_scale));
            pmt_calibration[key] = cal;
        }
    } else {
        ar & make_nvp("PMTCalibration",pmt_calibration);
    }
    ar & make_nvp("Rs",Rs);
    if(version == 0) {
        Rt = 1.0 - Rs;
    } else {
        ar & make_nvp("Rt",Rt);
    }
    ar & make_nvp("tau_s",tau_s);
    ar & make_nvp("tau_t",tau_t);
    ar & make_nvp("tau_other",tau_other);
    if(version >= 2) {
        ar & make_nvp("uv_absorption_a",uv_absorption_a);
        ar & make_nvp("uv_absorption_b",uv_absorption_b);
        ar & make_nvp("uv_absorption_d",uv_absorption_d);
        ar & make_nvp("uv_absorption_scaling",uv_absorption_scaling);
    } else {
        uv_absorption_a = 0.245319;
        uv_absorption_b = 111.817;
        uv_absorption_d = 5.8;
        uv_absorption_scaling = 0.0418985;
    }
}

std::ostream& operator<<(std::ostream& oss, const CCMPulseTimeDistributionParameters& c) {
    return(c.Print(oss));
}

std::ostream& CCMPulseTimeDistributionParameters::Print(std::ostream& os) const{
    os << "[ CCMPulseTimeDistributionParameters::"
        << "\n  mu :" << mu
        << "\n  sigma :" << sigma
        << "\n  fraction :" << fraction
        << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& oss, const CCMSimulationPMTCalibration& c) {
    return(c.Print(oss));
}

std::ostream& operator<<(std::ostream& os, const CCMSimulationPMTCalibrationMap& m) {
    os << "[ CCMSimulationPMTCalibrationMap::";
    for (const auto& it : m) {
        os << "\n    " << it.first << " : " << it.second;
    }
    os << " ]";
    return os;
}

std::ostream& CCMSimulationPMTCalibration::Print(std::ostream& os) const{
    os << "[ CCMSimulationPMTCalibration::"
        << "\n  pmt_efficiency :" << pmt_efficiency
        << "\n  pmt_spe_mu :" << pmt_spe_mu
        << "\n  pmt_spe_sigma :" << pmt_spe_sigma
        << "\n  pmt_spe_threshold :" << pmt_spe_threshold
        << "\n  main_pulse_mu :" << main_pulse_mu
        << "\n  main_pulse_sigma :" << main_pulse_sigma
        << "\n  late_pulses :" << late_pulses
        << " ]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const CCMSimulationCalibration& pe) {
    return(pe.Print(os));
}

std::ostream& CCMSimulationCalibration::Print(std::ostream& os) const{
    os << "[ CCMSimulationCalibration::"
        << "\n  Rs :" << Rs
        << "\n  Rt :" << Rt
        << "\n  Tau_s :" << tau_s
        << "\n  Tau_t :" << tau_t
        << "\n  Tau_other :" << tau_other
        << "\n  uv_absorption_a :" << uv_absorption_a
        << "\n  uv_absorption_b :" << uv_absorption_b
        << "\n  uv_absorption_d :" << uv_absorption_d
        << "\n  uv_absorption_scaling :" << uv_absorption_scaling
        << "\n  PMTCalibration :";
        for (const auto& it : pmt_calibration) {
            os << "\n    " << it.first << " : " << it.second;
        }
        os << " ]";
    return os;
}

I3_SPLIT_SERIALIZABLE(CCMSimulationCalibration);


