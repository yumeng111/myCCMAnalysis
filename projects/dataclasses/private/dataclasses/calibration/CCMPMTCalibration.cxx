
//
//  $Id$
//
//
#include <map>
#include <tuple>
#include <limits>
#include <icetray/serialization.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <icetray/I3Units.h>

#include <gsl/gsl_integration.h>


bool CCMSinglePulseParameters::operator==(CCMSinglePulseParameters const & other) const {
    return std::tie(
        peak_height,
        relative_peak_time,
        rise_time, // b1
        duration // b2
            ) == std::tie(
        other.peak_height,
        other.relative_peak_time,
        other.rise_time, // b1
        other.duration // b2
    );
}

bool CCMSPETemplate::operator==(CCMSPETemplate const & other) const {
    return pulse_parameters_ == other.pulse_parameters_;
}

bool CCMPMTCalibration::operator==(CCMPMTCalibration const & other) const {
    return std::tie(
        droopTimeConstant_,
        pmtGain_,
        pmtDeltaT_,
        relativePMTEff_,
        pulse_start_time_,
        pulse_end_time_,
        speTemplate_,
        combinedSPEFit_,
        meanPMTCharge_
        ) == std::tie(
        other.droopTimeConstant_,
        other.pmtGain_,
        other.pmtDeltaT_,
        other.relativePMTEff_,
        other.pulse_start_time_,
        other.pulse_end_time_,
        other.speTemplate_,
        other.combinedSPEFit_,
        other.meanPMTCharge_
    );

}

CCMPMTCalibration::CCMPMTCalibration() :
    droopTimeConstant_(std::numeric_limits<double>::quiet_NaN()),
    pmtGain_(std::numeric_limits<double>::quiet_NaN()),
    pmtDeltaT_(std::numeric_limits<double>::quiet_NaN()),
    relativePMTEff_(std::numeric_limits<double>::quiet_NaN()),
    pulse_start_time_(std::numeric_limits<double>::quiet_NaN()),
    pulse_end_time_(std::numeric_limits<double>::quiet_NaN()),
    speTemplate_(),
    combinedSPEFit_(),
    meanPMTCharge_(std::numeric_limits<double>::quiet_NaN())
{}

void CCMPMTCalibration::SetMeanPMTCharge(double mean_charge){
    meanPMTCharge_ = mean_charge;
}

double CCMPMTCalibration::GetMeanPMTCharge() const {
    return meanPMTCharge_;
}

void CCMPMTCalibration::SetSPEChargeDistribution(SPEChargeDistribution const & combined_fit) {
    combinedSPEFit_ = combined_fit;
}

SPEChargeDistribution const & CCMPMTCalibration::GetSPEChargeDistribution() const {
    return combinedSPEFit_;
}

void CCMPMTCalibration::SetSPETemplate(CCMSPETemplate const & spe_params) {
    speTemplate_ = spe_params;
}

CCMSPETemplate const & CCMPMTCalibration::GetSPETemplate() const {
    return speTemplate_;
}

// set pulse start time
void CCMPMTCalibration::SetPulseStartTime(double start_time){
    pulse_start_time_ = start_time;
}

// get pulse start time
double CCMPMTCalibration::GetPulseStartTime() const {
    return pulse_start_time_;
}

// set pulse end time
void CCMPMTCalibration::SetPulseEndTime(double end_time){
    pulse_end_time_ = end_time;
}

// get pulse end time
double CCMPMTCalibration::GetPulseEndTime() const {
    return pulse_end_time_;
}
// set droop time constant
void CCMPMTCalibration::SetDroopTimeConstant(double tau) {
    droopTimeConstant_ = tau;
}

// get droop time constant
double CCMPMTCalibration::GetDroopTimeConstant() const {
    return droopTimeConstant_;
}

// set pmt gain
void CCMPMTCalibration::SetPMTGain(double gain) {
    pmtGain_ = gain;
}

// get pmt gain
double CCMPMTCalibration::GetPMTGain() const {
    return pmtGain_;
}

// set pmt delta t
void CCMPMTCalibration::SetPMTDeltaT(double delta_t){
    pmtDeltaT_ = delta_t;
}

// get pmt delta t
double CCMPMTCalibration::GetPMTDeltaT() const {
    return pmtDeltaT_;
}

// set pmt relative efficiency
void CCMPMTCalibration::SetPMTRelEff(double rel_eff){
    relativePMTEff_ = rel_eff;
}

// get pmt relative efficiency
double CCMPMTCalibration::GetPMTRelEff() const {
    return relativePMTEff_;
}

// set PMT mean charge
void CCMPMTCalibration::SetPMTMeanCharge(double mean_charge){
    meanPMTCharge_ = mean_charge;
}

// get PMT mean charge
double CCMPMTCalibration::GetPMTMeanCharge() const {
    return meanPMTCharge_;
}
/*
 * Pulse template functions for use in simulation and reconstruction.
 *
 */

double CCMSPETemplate::EvaluateSinglePulse(CCMSinglePulseParameters const & single_pulse_params, double const & time){
    // let's evaluate a single pulse
    double peak_height = single_pulse_params.peak_height; 
    double relative_peak_time = single_pulse_params.relative_peak_time;
    double rise_time = single_pulse_params.rise_time;
    double duration = single_pulse_params.duration;

    // now we need to convert between these parameters and the c, t0, b1, and b2 parameters use for pulse shape
    double b2 = duration;
    double b1 = rise_time;
    double t0 = double(time) + relative_peak_time + (b1 * b2) * (std::log(b1) - std::log(b2)) / (b1 + b2);
    double c = peak_height / ( std::pow(b1, 8*b1/(b1+b2) ) * std::pow(b2, 8*b2/(b1+b2) ) / std::pow((b1+b2), 8) );

    // now let's evaluate the pulse in our time bin
    double w = c / (std::pow(std::exp(-1 * (time - t0)/b1) + std::exp((time - t0)/b2), 8));
    return w;
}

std::vector<CCMSinglePulseParameters> const & CCMSPETemplate::GetPulseParameters() const {
    return pulse_parameters_; 
}

void CCMSPETemplate::SetPulseParameters(std::vector<CCMSinglePulseParameters> const & params) {
    pulse_parameters_ = params;
}

double CCMSPETemplate::Evaluate(double time) const {
    double total = 0;
    for(size_t i =0; i<pulse_parameters_.size(); ++i) {
        total += EvaluateSinglePulse(pulse_parameters_[i], time);
    }
    return total;
}


namespace GSL{

    ///\brief Compute a one-dimensional integral using GSL.
    ///\param f Function to integrate.
    ///\param a Lower integration limit.
    ///\param b Upper integration limit.
    ///\param acc Accuracy parameter.
    ///\param max_iter Maximum number of iterations to perform the integral.
    template<typename FunctionType>
        double integrate(FunctionType&& f, double a, double b, double rtol=1e-7, unsigned int max_iter=10000, size_t memory_alloc=10000){
            using IntegrateWorkspace=std::unique_ptr<gsl_integration_workspace, void(*)(gsl_integration_workspace*)>;
            IntegrateWorkspace ws(gsl_integration_workspace_alloc(memory_alloc), &gsl_integration_workspace_free);

            using FPtr=decltype(&f);
            double (*wrapper)(double,void*)=[](double x, void* params)->double{
                auto& f=*static_cast<FPtr>(params);
                return(f(x));
            };

            double result, error;
            gsl_function F;
            F.function = wrapper;
            F.params = const_cast<typename std::remove_const<typename std::remove_reference<decltype(f)>::type>::type*>(&f);

            gsl_integration_qag(&F, a, b, 0, rtol, max_iter, GSL_INTEG_GAUSS15, ws.get(), &result, &error);

            return(result);
        }

} //namespace GSL

template <class Archive>
void CCMSPETemplate::serialize(Archive& ar, unsigned version) {
    ar & make_nvp("pulseParameters", pulse_parameters_);
}

std::ostream& operator<<(std::ostream& oss, const CCMSPETemplate& c) {
    oss << "[ CCMSPETemplate  :: " << std::endl
        << "       PulseParameters : " << c.GetPulseParameters() << std::endl
        << "]" ;
    return oss;
}

I3_SERIALIZABLE(CCMSPETemplate);
    
template <class Archive>
void CCMSinglePulseParameters::serialize(Archive& ar, unsigned version) {
    ar & make_nvp("peakHeight", peak_height);
    ar & make_nvp("peakTime", relative_peak_time);
    ar & make_nvp("riseTime", rise_time);
    ar & make_nvp("duration", duration);
}

std::ostream& operator<<(std::ostream& oss, const CCMSinglePulseParameters& c) {
    oss << "[ CCMSinglePulseParameters  :: " << std::endl
        << "       PeakHeight : " << c.peak_height << std::endl
        << "            PeakTime : " << c.relative_peak_time << std::endl
        << "          RiseTime : " << c.rise_time << std::endl
        << "      Duration : " << c.duration << std::endl
        << "]" ;
    return oss;
}

/*
std::ostream& operator<<(std::ostream& oss, const CCMSinglePulseParametersMap& m) {
    oss << "[ CCMSinglePulseParametersMap :: " << std::endl;
    CCMSinglePulseParametersMap::const_iterator iter = m.begin();
    for (;iter != m.end();iter++)
    {
        oss << "  " << iter->first << " : " << iter->second << std::endl;
    }
    oss << "]" ;
    return oss;
}*/

I3_SERIALIZABLE(CCMSinglePulseParameters);

template <class Archive>
void CCMPMTCalibration::serialize(Archive& ar, unsigned version) {
    ar & make_nvp("droopTimeConstant", droopTimeConstant_);
    ar & make_nvp("pmtGain",pmtGain_);
    ar & make_nvp("pmtDeltaT", pmtDeltaT_);
    ar & make_nvp("relativePMTEff", relativePMTEff_);
    ar & make_nvp("pulseStartTime", pulse_start_time_);
    ar & make_nvp("pulseEndTime", pulse_end_time_);
    ar & make_nvp("speParameters", speTemplate_);
    ar & make_nvp("combinedSPEFit", combinedSPEFit_);
    ar & make_nvp("meanPMTCharge", meanPMTCharge_);
}

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibration& c) {
    oss << "[ CCMPMTCalibration  :: " << std::endl
        << "       DroopTimeConstant : " << c.GetDroopTimeConstant() << std::endl
        << "            PMTGain : " << c.GetPMTGain() << std::endl
        << "          PMTDeltaT : " << c.GetPMTDeltaT() << std::endl
        << "      RelativePMTEff : " << c.GetPMTRelEff() << std::endl
        << "     PulseStartTime : " << c.GetPulseStartTime() << std::endl
        << "     PulseEndTime : " << c.GetPulseEndTime() << std::endl
        << "    SPETemplate : " << c.GetSPETemplate() << std::endl
        << "    SPEChargeDistribution : " << c.GetSPEChargeDistribution() << std::endl
        << "      MeanPMTCharge : " << c.GetPMTMeanCharge() << std::endl
        << "]" ;
    return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibrationMap& m) {
    oss << "[ CCMPMTCalibrationMap :: " << std::endl;
    CCMPMTCalibrationMap::const_iterator iter = m.begin();
    for (;iter != m.end();iter++)
    {
        oss << "  " << iter->first << " : " << iter->second << std::endl;
    }
    oss << "]" ;
    return oss;
}

I3_SERIALIZABLE(CCMPMTCalibration);
