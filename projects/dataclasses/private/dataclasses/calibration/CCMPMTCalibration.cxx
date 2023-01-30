
//
//  $Id$
//
//
#include <map>
#include <icetray/serialization.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <icetray/I3Units.h>

#include <gsl/gsl_integration.h>

CCMPMTCalibration::~CCMPMTCalibration() { }

CCMPMTCalibration::CCMPMTCalibration()
  : temperature_(NAN),
    pmtGain_(NAN),
    pmtDeltaT_(NAN),
    relativePMTEff_(NAN),
    meanPMTCharge_(NAN)
{
  droopTimeConstants_[0] = NAN;
  droopTimeConstants_[1] = NAN;
}


void
CCMPMTCalibration::SetTemperature(double temperature)
{
  temperature_ = temperature;
  SetTauParameters(tauparameters_);
}

// May need to change this for CCM
void
CCMPMTCalibration::SetTauParameters(TauParam tauparameters)
{
  tauparameters_ = tauparameters;

  droopTimeConstants_[0] = tauparameters_.P0 + tauparameters_.P1 /
      (1.0 + exp(-((temperature_-273.0)/tauparameters_.P2)));
  droopTimeConstants_[1] = tauparameters_.P3 + tauparameters_.P4 /
      (1.0 + exp(-((temperature_-273.0)/tauparameters_.P5)));
}


/*
 * Pulse template functions for use in simulation and reconstruction.
 *
 * Per discussion with C. Wendt Jan. 13, 2011, his FADC fits are offset 50 ns
 * and ATWD fits offset 5 ns. Other parameters from web page at
 * http://www.icecube.wisc.edu/~chwendt/dom-spe-waveform-shape/
 * or, in the case of the denominator of c, the result of numerical integrals.
 *
 * Templates for ATWD channels 1 and 2 were bootstrapped from the channel 0
 * using the method described here:
 * http://www.icecube.wisc.edu/~jvansanten/docs/atwd_pulse_templates/
 */

// These numbers will have to be updated for CCM, new SPE template

const double causalityShift = -11.5;  /* Nanoseconds from peak center to photon */

const SPETemplate PMTTemplate(25.12 / 71.363940160184669,
			       61.27 - 50 - causalityShift,
			       30.,
			       186.);

const SPETemplate PMTDroopTemplate(-2.8837584956162883,
				    0.57888025049064207,
				    0.81965713180496758,
				    0.04299648444652391);


CCMPMTCalibration::DroopedSPETemplate
CCMPMTCalibration::PMTPulseTemplate(bool droopy) const
{
  if(!droopy)
    return(DroopedSPETemplate(PMTTemplate));
  return(DroopedSPETemplate(PMTTemplate,
                            PMTDroopTemplate,
                            tauparameters_.TauFrac,
                            droopTimeConstants_[0],
                            droopTimeConstants_[1]));
}

bool CCMPMTCalibration::DroopedSPETemplate::operator==(const DroopedSPETemplate& templ) const
{
  return(pulse.c==templ.pulse.c &&
         pulse.x0==templ.pulse.x0 &&
         pulse.b1==templ.pulse.b1 &&
         pulse.b2==templ.pulse.b2 &&
         droopy==templ.droopy &&
         droop.pulse.c==templ.droop.pulse.c &&
         droop.pulse.x0==templ.droop.pulse.x0 &&
         droop.pulse.b1==templ.droop.pulse.b1 &&
         droop.pulse.b2==templ.droop.pulse.b2 &&
         droop.tauFrac==templ.droop.tauFrac &&
         droop.time1==templ.droop.time1 &&
         droop.time2==templ.droop.time2);
}

bool CCMPMTCalibration::DroopedSPETemplate::operator!=(const DroopedSPETemplate& templ) const{
  return !operator==(templ);
}

bool CCMPMTCalibration::DroopedSPETemplate::operator<(const DroopedSPETemplate& templ) const
{
#define DSPET_compare_member(member) \
  if(member<templ.member) return(true); \
  if(member>templ.member) return(false);

  DSPET_compare_member(pulse.c);
  DSPET_compare_member(pulse.x0);
  DSPET_compare_member(pulse.b1);
  DSPET_compare_member(pulse.b2);
  //we'll say that non-droopy templates are smaller than droopy templates
  if(droopy!=templ.droopy)
    return(!droopy);
  if(!droopy)
    return(false); //neither is droopy, so the templates are equal
  //both are droopy, so we go on comparing the components of the droop
  DSPET_compare_member(droop.pulse.c);
  DSPET_compare_member(droop.pulse.x0);
  DSPET_compare_member(droop.pulse.b1);
  DSPET_compare_member(droop.pulse.b2);
  DSPET_compare_member(droop.tauFrac);
  DSPET_compare_member(droop.time1);
  DSPET_compare_member(droop.time2);
  return(false); //templates are equal

#undef DSPET_compare_member
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
void
CCMPMTCalibration::serialize(Archive& ar, unsigned version)
{
  ar & make_nvp("temperature",temperature_);
  ar & make_nvp("pmtGain",pmtGain_);
  ar & make_nvp("pmtBaseline",pmtBaselineFit_);
  ar & make_nvp("pmtDeltaT", pmtDeltaT_);
  ar & make_nvp("tauparameters", tauparameters_);
  SetTauParameters(tauparameters_);
  ar & make_nvp("relativePMTEff", relativePMTEff_);
  ar & make_nvp("combinedSPEFit", combinedSPEFit_);
  ar & make_nvp("meanPMTCharge", meanPMTCharge_);
}

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibration& c)
{
  oss << "[ CCMPMTCalibration  :: " << std::endl
      << "       TauParameters : " << c.GetTauParameters() << std::endl
      << "         Temperature : " << c.GetTemperature() << std::endl
      << "            PMTGain : " << c.GetPMTGain() << std::endl
      << "     PMTBaselineFit : " << c.GetPMTBaselineFit() << std::endl
      << "          PMTDeltaT : " << c.GetPMTDeltaT() << std::endl
      << "      RelativePMTEff : " << c.GetRelativePMTEff() << std::endl
      << "CombinedSPEChargeDistribution : " << c.GetCombinedSPEChargeDistribution() << std::endl
      << "      MeanPMTCharge : " << c.GetMeanPMTCharge() << std::endl
      << "]" ;
  return oss;
}

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibrationMap& m)
{
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
