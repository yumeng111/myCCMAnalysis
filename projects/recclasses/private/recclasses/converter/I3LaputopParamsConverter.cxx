/**
 * copyright  (C) 2010
 * The Icecube Collaboration
 *
 * $Id: I3LaputopParamsConverter.cxx 124033 2014-10-01 19:49:45Z jgonzalez $
 *
 * @version $Revision: 124033 $
 * @date $LastChangedDate: 2014-10-01 15:49:45 -0400 (Wed, 01 Oct 2014) $
 * @author Fabian Kislat <fabian.kislat@desy.de>, last changed by $LastChangedBy: jgonzalez $
 */

#include "I3LaputopParamsConverter.h"
#include <icetray/I3Units.h>

I3TableRowDescriptionPtr I3LaputopParamsConverter::CreateDescription(const I3LaputopParams& params)
{
    I3TableRowDescriptionPtr desc(new I3TableRowDescription());

    // For now, I'll keep these essentially the same as we do now, but with a few additions and subtractions

    // Fit parameters
    desc->AddField<double> ("s125", "VEM", "Pulse height in 125m distance from the core (Fit parameter: log10(S125))");
    desc->AddField<double> ("beta", "", "The slope parameter of the lateral fit. (DLP function)");
    desc->AddField<double> ("age", "", "The age parameter of the shower. (NKG function)");
    desc->AddField<double> ("A", "ns/m^2", "KR's curvature parabola-A free parameter (modified gausspar)");
    desc->AddField<double> ("D", "", "KR's curvature gaussian-D free parameter (modified gausspar)");
    desc->AddField<double> ("N", "ns", "KR's curvature gaussian-N free parameter (modified gausspar)");
    // Errors
    desc->AddField<double> ("xc_err", "m", "Uncertainty of the core X position");
    desc->AddField<double> ("yc_err", "m", "Uncertainty of the core Y position");
    desc->AddField<double> ("ny_err", "m", "Uncertainty of the direction unit vector X component");
    desc->AddField<double> ("nx_err", "m", "Uncertainty of the drection unit vector Y component");
    desc->AddField<double> ("tc_err", "m", "Uncertainty of the core time");
    desc->AddField<double> ("log10_s125_err", "1", "Fit uncertainty of log10(S125/VEM)");
    desc->AddField<double> ("beta_err", "", "Uncertainty of Beta");
    desc->AddField<double> ("A_err", "", "Uncertainty of curvature-A");
    desc->AddField<double> ("D_err", "", "Uncertainty of curvature-D");
    desc->AddField<double> ("N_err", "", "Uncertainty of curvature-N");

    // Convenient extras
    desc->AddField<double> ("s50", "VEM", "Pulse height in 50m distance from the core");
    desc->AddField<double> ("s70", "VEM", "Pulse height in 70m distance from the core");
    desc->AddField<double> ("s80", "VEM", "Pulse height in 80m distance from the core");
    desc->AddField<double> ("s100", "VEM", "Pulse height in 100m distance from the core");
    desc->AddField<double> ("s150", "VEM", "Pulse height in 150m distance from the core");
    desc->AddField<double> ("s180", "VEM", "Pulse height in 180m distance from the core");
    desc->AddField<double> ("s250", "VEM", "Pulse height in 250m distance from the core");
    desc->AddField<double> ("s500", "VEM", "Pulse height in 500m distance from the core");
    desc->AddField<double> ("llh", "", "The Log-Likelihood of the fit, divided by the degrees of freedom");
    desc->AddField<double> ("llh_silent", "", "The no-hit part of the likelihood");
    desc->AddField<double> ("chi2_ldf", "", "The Chi square, divided by the degrees of freedom, for the charge fit");
    desc->AddField<double> ("chi2_time", "", "The Chi square, divided by the degrees of freedom, for the timing fit");
    desc->AddField<int32_t>("ndf", "", "The number of degrees of freedom (after all cuts)");
    desc->AddField<double> ("rlogl", "", "The reduced log-likelihood");
    desc->AddField<int32_t>("nmini", "", "The number of minimization steps required to converge");
    return desc;
}
    
size_t I3LaputopParamsConverter::FillRows(const I3LaputopParams& params,
						I3TableRowPtr rows)
{
 
  // Nonstd fit parameters
  rows->Set<double> ("s125", (params.Has(Laputop::Parameter::Log10_S125) ? pow(10, params.GetValue(Laputop::Parameter::Log10_S125)) : NAN));
  rows->Set<double> ("beta", (params.Has(Laputop::Parameter::Beta) ? params.GetValue(Laputop::Parameter::Beta) : NAN));
  rows->Set<double> ("age", (params.Has(Laputop::Parameter::Age) ? params.GetValue(Laputop::Parameter::Age) : NAN));
  rows->Set<double> ("A", (params.Has(Laputop::Parameter::CurvParabA) ? params.GetValue(Laputop::Parameter::CurvParabA) : NAN));
  rows->Set<double> ("D", (params.Has(Laputop::Parameter::CurvGaussD) ? params.GetValue(Laputop::Parameter::CurvGaussD) : NAN));
  rows->Set<double> ("N", (params.Has(Laputop::Parameter::CurvGaussN) ? params.GetValue(Laputop::Parameter::CurvGaussN) : NAN));
  // Errors
  rows->Set<double> ("xc_err", (params.Has(Laputop::Parameter::Xc) ? params.GetError(Laputop::Parameter::Xc) : NAN));
  rows->Set<double> ("yc_err", (params.Has(Laputop::Parameter::Yc) ? params.GetError(Laputop::Parameter::Yc) : NAN));
  rows->Set<double> ("nx_err", (params.Has(Laputop::Parameter::Nx) ? params.GetError(Laputop::Parameter::Nx) : NAN));
  rows->Set<double> ("ny_err", (params.Has(Laputop::Parameter::Ny) ? params.GetError(Laputop::Parameter::Ny) : NAN));
  rows->Set<double> ("tc_err", (params.Has(Laputop::Parameter::Tc) ? params.GetError(Laputop::Parameter::Tc) : NAN));
  rows->Set<double> ("log10_s125_err", (params.Has(Laputop::Parameter::Log10_S125) ? params.GetError(Laputop::Parameter::Log10_S125) : NAN));
  rows->Set<double> ("beta_err", (params.Has(Laputop::Parameter::Beta) ? params.GetError(Laputop::Parameter::Beta) : NAN));
  rows->Set<double> ("A_err", (params.Has(Laputop::Parameter::CurvParabA) ? params.GetError(Laputop::Parameter::CurvParabA) : NAN));
  rows->Set<double> ("D_err", (params.Has(Laputop::Parameter::CurvGaussD) ? params.GetError(Laputop::Parameter::CurvGaussD) : NAN));
  rows->Set<double> ("N_err", (params.Has(Laputop::Parameter::CurvGaussN) ? params.GetError(Laputop::Parameter::CurvGaussN) : NAN));
  // Convenient extras
  rows->Set<double> ("s50", params.ExpectedSignal(50.*I3Units::m));
  rows->Set<double> ("s70", params.ExpectedSignal(70.*I3Units::m));
  rows->Set<double> ("s80", params.ExpectedSignal(80.*I3Units::m));
  rows->Set<double> ("s100", params.ExpectedSignal(100.*I3Units::m));
  rows->Set<double> ("s150", params.ExpectedSignal(150.*I3Units::m));
  rows->Set<double> ("s180", params.ExpectedSignal(180.*I3Units::m));
  rows->Set<double> ("s250", params.ExpectedSignal(250.*I3Units::m));
  rows->Set<double> ("s500", params.ExpectedSignal(500.*I3Units::m));
  rows->Set<double> ("llh", params.logLikelihood_);
  rows->Set<double> ("llh_silent", params.logLikelihood_Silent_);
  rows->Set<double> ("chi2_ldf", params.chi2_LDF_);
  rows->Set<double> ("chi2_time", params.chi2_Time_);
  rows->Set<int32_t>("ndf", params.ndf_); // Not sure about this one... how can ndf_Time be different?
  rows->Set<double> ("rlogl", params.logLikelihood_ / params.ndf_);
  rows->Set<int32_t>("nmini", params.nMini_);
  
  return 1;
}
