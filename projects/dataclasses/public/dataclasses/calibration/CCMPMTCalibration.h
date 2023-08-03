/**
 *
 * Definition of CCMPMTCalibration class
 *
 */

#ifndef CCMPMTCalibration_H_INCLUDED
#define CCMPMTCalibration_H_INCLUDED

#include <string>
#include <map>
#include <vector>
#include <stdint.h>
#include <sstream>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include <dataclasses/calibration/I3DOMCalibration.h>
#include <dataclasses/external/CompareFloatingPoint.h>

#include <icetray/I3Units.h>
#include <dataclasses/Utility.h>
#include <icetray/CCMPMTKey.h>



static const unsigned ccmsinglepulseparameters_version_ = 0;

struct CCMSinglePulseParameters {
    double peak_height;
    double relative_peak_time;
    double rise_time; // b1
    double duration; // b2

    bool operator==(CCMSinglePulseParameters const & other) const;

    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

I3_CLASS_VERSION(CCMSinglePulseParameters, ccmsinglepulseparameters_version_);

std::ostream& operator<<(std::ostream& oss, const CCMSinglePulseParameters& c);

static const unsigned ccmspetemplate_version_ = 0;

class CCMSPETemplate {
private:
    std::vector<CCMSinglePulseParameters> pulse_parameters_;
public:
    static double EvaluateSinglePulse(CCMSinglePulseParameters const & single_pulse_params, double const & time);
    std::vector<CCMSinglePulseParameters> const & GetPulseParameters() const;
    void SetPulseParameters(std::vector<CCMSinglePulseParameters> const & params);
    double Evaluate(double time) const;
    double GetTemplatePeak() const;

    bool operator==(CCMSPETemplate const & other) const;

    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

I3_CLASS_VERSION(CCMSPETemplate, ccmspetemplate_version_);

std::ostream& operator<<(std::ostream& oss, const CCMSPETemplate& c);

static const unsigned ccmpmtcalibration_version_ = 0;

class CCMPMTCalibration {
public:
    CCMPMTCalibration();
    bool operator==(CCMPMTCalibration const & other) const;

private:
    // parameters useful for calibration
    double droopTimeConstant_;
    double pmtGain_;
    double pmtDeltaT_;
    double relativePMTEff_;

    double pulse_start_time_;
    double pulse_end_time_;
    CCMSPETemplate speTemplate_;

    SPEChargeDistribution combinedSPEFit_;
    double meanPMTCharge_;
public:
    void SetPulseStartTime(double start_time);
    double GetPulseStartTime() const;

    void SetPulseEndTime(double end_time);
    double GetPulseEndTime() const;
    
    void SetDroopTimeConstant(double tau);
    double GetDroopTimeConstant() const;

    void SetPMTGain(double gain);
    double GetPMTGain() const;

    void SetPMTDeltaT(double delta_t);
    double GetPMTDeltaT() const;

    void SetPMTRelEff(double rel_eff);
    double GetPMTRelEff() const;

    void SetMeanPMTCharge(double mean_charge);
    double GetMeanPMTCharge() const;

    void SetSPEChargeDistribution(SPEChargeDistribution const & spe_fit);
    SPEChargeDistribution const & GetSPEChargeDistribution() const;

    void SetSPETemplate(CCMSPETemplate const & speTemplate_);
    CCMSPETemplate const & GetSPETemplate() const;

    void SetPMTMeanCharge(double mean_charge);
    double GetPMTMeanCharge() const;

    template <class Archive>
    void serialize(Archive& ar, unsigned version);
};

I3_CLASS_VERSION(CCMPMTCalibration, ccmpmtcalibration_version_);

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibration& c);

typedef std::map<CCMPMTKey, CCMPMTCalibration> CCMPMTCalibrationMap;
I3_POINTER_TYPEDEFS(CCMPMTCalibrationMap);

typedef std::vector<CCMSinglePulseParameters> CCMSinglePulseParametersSeries;
I3_POINTER_TYPEDEFS(CCMSinglePulseParametersSeries);

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibrationMap& m);

#endif // CCMPMTCalibration_H_INCLUDED
