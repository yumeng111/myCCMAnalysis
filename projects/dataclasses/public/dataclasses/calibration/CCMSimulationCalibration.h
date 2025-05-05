/**
 * copyright  (C) 2013
 * the icecube collaboration
 * @version $Id: $
 */

#ifndef CCMSimulationCalibration_H_INCLUDED
#define CCMSimulationCalibration_H_INCLUDED

#include <dataclasses/I3Position.h>
#include <dataclasses/I3Direction.h>
#include <dataclasses/Utility.h>
#include <dataclasses/I3Map.h>
#include <dataclasses/I3Vector.h>
#include <icetray/I3FrameObject.h>

#include <string>
#include <iostream>
#include <sstream>

#include <icetray/CCMPMTKey.h>

static const unsigned ccmsimulationcalibration_version_ = 4;
static const unsigned ccmsimulationpmtcalibration_version_ = 0;
static const unsigned ccmlatepulseparameters_version_ = 0;

struct CCMPulseTimeDistributionParameters {
    double mu;
    double sigma;
    double fraction;
    CCMPulseTimeDistributionParameters() = default;
    CCMPulseTimeDistributionParameters(CCMPulseTimeDistributionParameters const & pulse) = default;
    CCMPulseTimeDistributionParameters(double mu, double sigma, double fraction) : mu(mu), sigma(sigma), fraction(fraction) {}
    CCMPulseTimeDistributionParameters(std::array<double, 3> const & pulse) : mu(pulse[0]), sigma(pulse[1]), fraction(pulse[2]) {}
    // Initializer list
    CCMPulseTimeDistributionParameters(std::initializer_list<double> const & pulse) {
        auto it = pulse.begin();
        mu = *it++;
        sigma = *it++;
        fraction = *it++;
    }

    bool operator==(const CCMPulseTimeDistributionParameters& rhs) const {
        return mu == rhs.mu
            && sigma == rhs.sigma
            && fraction == rhs.fraction;
    }

    std::ostream& Print(std::ostream&) const;

    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

I3_CLASS_VERSION(CCMPulseTimeDistributionParameters, ccmlatepulseparameters_version_);
typedef std::vector<CCMPulseTimeDistributionParameters> CCMPulseTimeDistributionParametersSeries;

std::ostream& operator<<(std::ostream& oss, const CCMPulseTimeDistributionParameters& c);

I3_POINTER_TYPEDEFS(CCMPulseTimeDistributionParametersSeries);

struct CCMSimulationPMTCalibration {
public:
    // Parameters of the SPE Charge distribution
    double pmt_efficiency;
    double pmt_spe_mu;
    double pmt_spe_sigma;
    double pmt_spe_threshold;

    double main_pulse_mu;
    double main_pulse_sigma;
    std::vector<CCMPulseTimeDistributionParameters> late_pulses;

    bool operator==(const CCMSimulationPMTCalibration& rhs) const {
        return pmt_efficiency == rhs.pmt_efficiency
            && pmt_spe_mu == rhs.pmt_spe_mu
            && pmt_spe_sigma == rhs.pmt_spe_sigma
            && pmt_spe_threshold == rhs.pmt_spe_threshold
            && main_pulse_mu == rhs.main_pulse_mu
            && main_pulse_sigma == rhs.main_pulse_sigma
            && late_pulses == rhs.late_pulses;
    }

    std::ostream& Print(std::ostream&) const;

    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

I3_CLASS_VERSION(CCMSimulationPMTCalibration, ccmsimulationpmtcalibration_version_);

std::ostream& operator<<(std::ostream& oss, const CCMSimulationPMTCalibration& c);

typedef std::map<CCMPMTKey, CCMSimulationPMTCalibration> CCMSimulationPMTCalibrationMap;
I3_POINTER_TYPEDEFS(CCMSimulationPMTCalibrationMap);

std::ostream& operator<<(std::ostream& oss, const CCMSimulationPMTCalibrationMap& m);

class CCMSimulationCalibration: public I3FrameObject {
public:
    CCMSimulationPMTCalibrationMap pmt_calibration;

    double Rs;
    double Rt;
    double tau_s;
    double tau_t;
    double tau_other;

    double uv_absorption_a;
    double uv_absorption_b;
    double uv_absorption_d;
    double uv_absorption_scaling;

    double normalization;

    SET_LOGGER("CCMSimulationCalibration");

    bool operator==(const CCMSimulationCalibration& rhs) const {
        return pmt_calibration == rhs.pmt_calibration
            && Rs == rhs.Rs
            && Rt == rhs.Rt
            && tau_s == rhs.tau_s
            && tau_t == rhs.tau_t
            && tau_other == rhs.tau_other
            && uv_absorption_a == rhs.uv_absorption_a
            && uv_absorption_b == rhs.uv_absorption_b
            && uv_absorption_d == rhs.uv_absorption_d
            && uv_absorption_scaling == rhs.uv_absorption_scaling
            && normalization == rhs.normalization;
    }

    CCMSimulationCalibration();

    std::ostream& Print(std::ostream&) const;

    void SetPMTCalibration(CCMPMTKey key, const CCMSimulationPMTCalibration& pmt_cal);
    CCMSimulationPMTCalibration const & GetPMTCalibration(CCMPMTKey key) const;
    //CCMSimulationPMTCalibration & GetPMTCalibration(CCMPMTKey key);

    void SetPMTCalibration(CCMSimulationPMTCalibrationMap const & pmt_cal);
    CCMSimulationPMTCalibrationMap const & GetPMTCalibration() const;
    //CCMSimulationPMTCalibrationMap & GetPMTCalibration();

friend class icecube::serialization::access;
    template<class Archive> void save(Archive& ar, unsigned version) const;
    template<class Archive> void load(Archive& ar, unsigned version);
    I3_SERIALIZATION_SPLIT_MEMBER();
};

I3_CLASS_VERSION(CCMSimulationCalibration,ccmsimulationcalibration_version_);

std::ostream& operator<<(std::ostream&, const CCMSimulationCalibration&);

I3_POINTER_TYPEDEFS(CCMSimulationCalibration);
#endif
