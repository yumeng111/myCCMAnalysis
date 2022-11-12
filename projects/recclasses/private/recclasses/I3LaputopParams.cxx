/**
 * Copyright (C) 2015
 * The IceCube collaboration
 * ID: $Id$
 *
 * @file I3LaputopParams.cxx
 * @version $Rev$
 * @date $Date$
 * @author $Author$
 */

#include "recclasses/I3LaputopParams.h"
#include "recclasses/Utility.h"
#include <icetray/serialization.h>
#include <icetray/I3Units.h>
#include <dataclasses/physics/I3Particle.h>
#include <dataclasses/I3Position.h>
#include <dataclasses/ostream_overloads.hpp>

#include <iomanip>
#include <limits>
#include <cmath>


using namespace Laputop;

namespace {
  static const double qNaN = std::numeric_limits<double>::quiet_NaN();
  static const double deps = std::numeric_limits<double>::epsilon();

  // ----------- A copy of the DLP function from toprec ------------------
  // We could just call this directly using somethine like:
  // LateralFitFunctions::top_ldf_dlp(100.*I3Units::m, par)
  // ... but that would require dependence on toprec, and we don't want that, right?
  double top_ldf_dlp(double r, double log10_s125, double beta) {

    // These stolen from LateralFitFunctions.h
    double R0_PARAM = 125.0 * I3Units::m;  // parameter of the fit
    double KAPPA = 0.30264;                // constant of the DLP function (set to zero for powerlaw)

    // Now the actual function... adapted from LateralFitFunctions.cxx
    double local_x = log10(r/R0_PARAM);
    return log10_s125 - beta * local_x - KAPPA*local_x*local_x;
  }


  // ------------- The "gausspar" family of CURVATURE FUNCTIONS ---------------
  // also copied from from toprec/LateralFitFunctions,
  // but they should live here now. --KR
  // The public functions "ExpectedShowerFrontDelay" and "ExpectedShowerFrontDelayError" will call these
  
  // This is the general form of all "gausspar" curvatures
  // It can describe the whole family, as well as be used when 
  // these parameters are fit as free parameters
  double top_curv_gausspar_generic(double r, double N, double S, double A) {
    return N * (exp(-r*r/(S*S)) - 1.) - A*r*r; 
  }
  
  // The "IT-26 default" function:
  static const double CURV_GAUSSPAR_N = 19.41;
  static const double CURV_GAUSSPAR_D = 118.1;
  static const double CURV_GAUSSPAR_A = 4.823e-4;
  double top_curv_gausspar(double r){
    // this must be a function that is never NAN!
    // 118.1*118.1 = 2 * sigma*sigma with sigma = 83.5 (formalism in NIM)
    //return 19.41 * (exp(-r*r/(118.1*118.1)) - 1.) - 4.823e-4*r*r; // from TW
    return top_curv_gausspar_generic(r, CURV_GAUSSPAR_N, CURV_GAUSSPAR_D, CURV_GAUSSPAR_A);
  }
  
  // Emily's version of the function:
  static const double CURV_EMILY_N = 9.7233;
  static const double CURV_EMILY_D = 45.2176;
  static const double CURV_EMILY_A = 4.897e-4;
  double top_curv_gausspar_emily(double r){
    return top_curv_gausspar_generic(r, CURV_EMILY_N, CURV_EMILY_D, CURV_EMILY_A);    
  }

  // Temporary: a version with only "freeA"
  // (this should be made generic in the future)
  //double top_curv_gausspar_freeA(double r, double A){
  //  return top_curv_gausspar_generic(r, CURV_GAUSSPAR_N, CURV_GAUSSPAR_S, A);    
  //}

  // The sigma functions 
  double top_curv_sigma(double r) {
    return 2.92 + 3.77e-4*r*r; // from TW
  }

  double top_curv_sigma_emily(double r) {
    return 0.777125 + 1.8623e-4*r*r; // from TW
  }


  //--------------------------------------------------------------------------


  template <typename T>
  T sqr(T t) { return t*t; }


  double calibration(unsigned n,
                     const double* table_cosz,
                     const double* table_par,
                     double cosz, double log10_s125) {
    unsigned i = 0;
    for (; i < n && table_cosz[i] <= cosz; ++i);
    --i;
    if (i >= n || std::isnan(cosz))
      return qNaN;
    return std::pow(10.0, table_par[2 * i] + table_par[2 * i + 1] * log10_s125) * I3Units::GeV;
  }
}


I3LaputopParams::~I3LaputopParams() {}

template <class Archive>
void I3LaputopParams::serialize(Archive& ar, unsigned version)
{
  if (version > i3laputopparams_version_)
    log_fatal("Attempting to read version %u from file but running version %u of I3LaputopParams.",
              version, i3laputopparams_version_);

  // no fix for reading version 1, was not used in production

  ar & make_nvp("I3FrameObject", base_object<I3FrameObject>(*this));
  ar & make_nvp("ParameterStorage", base_object<ParameterStorageType>(*this));
  ar & make_nvp("curLDFType", typeLDF_);
  ar & make_nvp("curFrontDelayType", typeFrontDelay_);
  ar & make_nvp("curEnergyType", typeEnergy_);
  ar & make_nvp("zPos", zPos_);
  ar & make_nvp("logLikelihood", logLikelihood_);
  ar & make_nvp("logLikelihood_Silent", logLikelihood_Silent_) ;
  ar & make_nvp("chi2_LDF", chi2_LDF_);
  ar & make_nvp("chi2_Time", chi2_Time_);

  if (version == 2) {
    double dummy = 0;
    ar & make_nvp("ndf_LDF", ndf_);
    ar & make_nvp("ndf_Time", dummy);      
    chi2_LDF_ = qNaN; // cannot be recovered
    chi2_Time_ = qNaN; // cannot be recovered
  } else {
    ar & make_nvp("ndf", ndf_);
  }

  ar & make_nvp("nMini", nMini_);

  ar & make_nvp("userData", userData_);
}

I3_SERIALIZABLE(I3LaputopParams);

void I3LaputopParams::FillFromI3Particle(const I3Particle& p)
{
  SetValue(Parameter::Xc, p.GetPos().GetX());
  SetValue(Parameter::Yc, p.GetPos().GetY());
  zPos_ = p.GetPos().GetZ();

  SetValue(Parameter::Tc, p.GetTime());
  SetValue(Parameter::Nx, p.GetDir().GetX());
  SetValue(Parameter::Ny, p.GetDir().GetY());
}

double I3LaputopParams::ExpectedSignal(double r, double xi,
                                       LDF::Enum type) const
{
  if (type == LDF::None)
    type = typeLDF_;

  switch (type) {
    case LDF::DLP:
      // from toprec/private/laputop/LateralFitFunctions.cxx
      return std::pow(10.,
        top_ldf_dlp(r,
                    GetValue(Laputop::Parameter::Log10_S125),
                    GetValue(Laputop::Parameter::Beta)));
    default: break; // do nuthin'
  }

  return qNaN;
}

double I3LaputopParams::ExpectedSignalError(double r, double xi,
                                            LDF::Enum type) const
{
  if (type == LDF::None)
    type = typeLDF_;

  switch (type) {
    case LDF::DLP: {
      const double s = ExpectedSignal(r, xi, type);
      const double lgs = std::min(std::log10(s), 2.077);
      const double lgs_err =  
        std::pow(10, lgs < 0.340 ?
                 -0.552 - 0.078 * lgs :
                 -0.373 - 0.658 * lgs + 0.158 * lgs*lgs);
      // y = 10^x
      // sy = |dy/dx| sx = 10^x log(10) * sx
      return s * std::log(10.0) * lgs_err;
    }
    default: break; // do nuthin'
  }

  return qNaN;
}

double I3LaputopParams::ExpectedShowerFrontDelay(double r, double xi,
                                                 FrontDelay::Enum type) const
{
  if (type == FrontDelay::None)
    type = typeFrontDelay_;

  switch (type) {
    case FrontDelay::GaussParabola:
      return top_curv_gausspar(r);
    case FrontDelay::GaussParabolaEmily:
      return top_curv_gausspar_emily(r);
    case FrontDelay::GaussParabolaFromParams:
      { double a = Has(Parameter::CurvParabA) ? GetValue(Parameter::CurvParabA) : CURV_GAUSSPAR_A;
	double d = Has(Parameter::CurvGaussD) ? GetValue(Parameter::CurvGaussD) : CURV_GAUSSPAR_D;
	double n = Has(Parameter::CurvGaussN) ? GetValue(Parameter::CurvGaussN) : CURV_GAUSSPAR_N;
	return top_curv_gausspar_generic(r, n, d, a);
      }
    default: break; // do nuthin'
  }

  return qNaN;
}

double I3LaputopParams::ExpectedShowerFrontDelayError(double r, double xi,
                                                      FrontDelay::Enum type) const
{
  if (type == FrontDelay::None)
    type = typeFrontDelay_;

  switch (type) {
    case FrontDelay::GaussParabola:
      return top_curv_sigma(r);
    case FrontDelay::GaussParabolaEmily:
      return top_curv_sigma_emily(r);
    case FrontDelay::GaussParabolaFromParams:
      return top_curv_sigma(r);   // for now, just use the IT-26 one 
    default: break; // do nuthin'
  }

  return qNaN;
}

double I3LaputopParams::Energy(Energy::Enum type) const
{
  if (type == Energy::None)
    type = typeEnergy_;

  const double cosz = std::cos(GetValue(Par::Theta));

  const double lgs125 = GetValue(Parameter::Log10_S125);

  const double tcosz[] = {0.80, 0.85, 0.90, 0.95};
  switch (type) {
    case Energy::IC73SpectrumPaper: {
      // H4a
      const double tpar[] = { 6.182, 0.914,
                              6.117, 0.921,
                              6.062, 0.929,
                              6.018, 0.938 };
      return calibration(4, tcosz, tpar, cosz, lgs125);
    }
    case Energy::IC73SpectrumPaperProton: {
      const double tpar[] = { 6.139, 0.923,
                              6.081, 0.936,
                              6.034, 0.948,
                              5.998, 0.962 };
      return calibration(4, tcosz, tpar, cosz, lgs125);
    }
    case Energy::IC73SpectrumPaperIron: {
      const double tpar[] = { 6.288, 0.878,
                              6.202, 0.888,
                              6.130, 0.900,
                              6.069, 0.913 };
      return calibration(4, tcosz, tpar, cosz, lgs125);
    }

     case Energy::ICRC2015_H4a_E27: {
      const double tpar[] = { 6.177271, 0.907456,
                              6.109777, 0.914971,
                              6.054677, 0.923860,
                              6.010569, 0.933316 };
      return calibration(4, tcosz, tpar, cosz, lgs125);
    }
     
    default: break; // do nuthin'
  }
  return qNaN;
}

bool I3LaputopParams::Has(Laputop::Parameter::Enum par) const
{
  switch (par) {
    case Par::Theta:
      return Base::Has(Par::Theta) || 
        (Base::Has(Par::Nx) && Base::Has(Par::Ny));
    case Par::Phi:
      return Base::Has(Par::Phi) ||
        (Base::Has(Par::Nx) && Base::Has(Par::Ny));
    default:
      return Base::Has(par);
  }
  assert(!"never arrive here");
  return false;
}

double I3LaputopParams::GetValue(Laputop::Parameter::Enum par) const
{
  switch (par) {
    case Par::Theta: {
      if (Base::Has(Par::Theta))
        return Base::GetValue(Par::Theta);
      // compute on the fly
      const double nx = Base::GetValue(Par::Nx);
      const double ny = Base::GetValue(Par::Ny);
      const double nz = std::sqrt(std::max(0.0, 1.0 - sqr(nx) - sqr(ny)));
      return std::acos(nz);
    }
    case Par::Phi: {
      if (Base::Has(Par::Phi))
        return Base::GetValue(Par::Phi);
      // compute on the fly
      const double nx = Base::GetValue(Par::Nx);
      const double ny = Base::GetValue(Par::Ny);
      return std::atan2(ny, nx);
    }
    default:
      return Base::GetValue(par);
  }
  assert(!"never arrive here");
  return qNaN;
}

double I3LaputopParams::GetError(Laputop::Parameter::Enum par) const
{
  return std::sqrt(GetCovariance(par, par));
}

double I3LaputopParams::GetCovariance(Laputop::Parameter::Enum par1,
                                      Laputop::Parameter::Enum par2) const
{
  if (Base::Has(par1) && Base::Has(par2))
    return Base::GetCovariance(par1, par2);

  double c12 = 0.0;
  for (int k = 0; k < int(Par::size); ++k) {    
    for (int l = 0; l < int(Par::size); ++l) {
      if (Base::Has(Par::Enum(k)) && Base::Has(Par::Enum(l))) {
        const double d1k = Jacobi(par1, Par::Enum(k));
        const double d2l = Jacobi(par2, Par::Enum(l));
        c12 += d1k * d2l * Base::GetCovariance(Par::Enum(k), Par::Enum(l));        
      }
    }
  }
  return c12;
}

double I3LaputopParams::Jacobi(Laputop::Parameter::Enum par1,
                               Laputop::Parameter::Enum par2) const
{
  // Jacobi(par1, par2) = dpar1/dpar2
  switch (par1) {
    case Par::Theta: {
      // Based on: theta = arccos(sqrt(1 - nx^2 - ny^2))
      const double nx = Base::GetValue(Par::Nx);
      const double ny = Base::GetValue(Par::Ny);
      const double nz2 = std::max(0.0, 1.0 - sqr(nx) - sqr(ny));
      const double c = std::sqrt(nz2 * (sqr(nx) + sqr(ny)));
      switch (par2) {
        case Par::Nx: return nx / (c + deps); // don't divide by zero
        case Par::Ny: return ny / (c + deps); // don't divide by zero
        default: return par1 == par2;
      }
    } break;
    case Par::Phi: {
      // Based on: phi = arctan(ny/nx)
      const double nx = Base::GetValue(Par::Nx);
      const double ny = Base::GetValue(Par::Ny);
      const double c = sqr(nx) + sqr(ny);
      switch (par2) {
        case Par::Nx: return -ny / (c + deps); // don't divide by zero
        case Par::Ny: return  nx / (c + deps); // don't divide by zero
        default: return par1 == par2;
      }
    } break;
    case Par::Nx: {
      switch (par2) {
        case Par::Theta: assert(!"not implemented"); return qNaN;
        case Par::Phi: assert(!"not implemented"); return qNaN;
        default: return par1 == par2;        
      }
    } break;
    case Par::Ny: {
      switch (par2) {
        case Par::Theta: assert(!"not implemented"); return qNaN;
        case Par::Phi: assert(!"not implemented"); return qNaN;
        default: return par1 == par2;        
      }
    } break;
    default: return par1 == par2;
  }
  assert(!"never arrive here");
  return qNaN;
}

double I3LaputopParams::GetAngularResolution() const {
  // We compute this from the covariance matrix of 
  // from Theta and Phi, except if the direction is almost
  // vertical, where the covariance matrix of Theta and Phi
  // becomes degenerate because of the pole. For those cases, 
  // where we use an approximation based directly on Nx and Ny.

  // To compute the angular resolution, we add the major and
  // minor axis of covariance matrix in theta and phi.
  // For a 2x2-matrix we can compute this sum without a library.
  // We need to solve the eigenvalue problem det(A - lambda I) = 0.

  // For cov = ((a c), (c b)) we obtain the condition
  // lambda^2 - (a + b) lambda + (a b - c^2) = 0, which we
  // solve with the pq-formula with p = a + b. We are only
  // interested in the sum of the solutions var_r1 and var_r2,
  // so the answer turns out to be simply var_r1 + var_r2 = p.
  const double x = Base::GetValue(Par::Nx);
  const double y = Base::GetValue(Par::Ny);  
  const double cxx = Base::GetCovariance(Par::Nx, Par::Nx);
  const double cyy = Base::GetCovariance(Par::Ny, Par::Ny);
  double p = 0.0;
  if (sqr(x) < cxx && sqr(y) < cyy) {
    // near-vertical case, an approximation that's better than the code below
    p = cxx + cyy;
  } else {
    // standard case
    const double sin_theta_sqr = 1.0 - sqr(x) - sqr(y);
    const double ctt = GetCovariance(Par::Theta, Par::Theta);
    const double cpp = GetCovariance(Par::Phi, Par::Phi);
    p = ctt + sin_theta_sqr * cpp;
  }

  // We compute a circular radius now, which covers the same confidence
  // region as the ellipse given by var_r1, var_r2. For near vertical
  // showers, var_r1 is close to var_r2 and this simplification is good.
  const double var_r = 0.5 * p; // since p = var_r1 + var_r2

  // Taken individually, var_r1 and var_r2 each already cover a
  // 68 % confidence interval along their respective axis in 1d.
  // In order to obtain a 2d confidence region for the combination
  // of two random variables with 68 % coverage, we need to apply a
  // factor given by the chi2-distribution for 2 degrees of freedom.
  const double f = 2.27886856638; // scipy.stats.chi2(2).ppf(0.68)
  const double sigma_arc = std::sqrt(f * var_r);
  // Finally, we need to convert the arc length into an angle via
  // angle = arc length / radius, but we computed all this on
  // the unit-sphere, whose radius is 1 and so we are done.
  return sigma_arc;
}

std::ostream& operator<<(std::ostream& os, const I3LaputopParams& p) {
  return(p.Print(os));
}

std::ostream& I3LaputopParams::Print(std::ostream& os) const{
  os << "[ I3LaputopParams::";

  #define PAROUT(Enum)             \
    if (Has(Parameter::Enum)) {                  \
      os << "\n" << std::setw(20) << #Enum << ": " \
         << GetValue(Parameter::Enum) << " +/- " \
         << GetError(Parameter::Enum);           \
    }
  PAROUT(Log10_S125);
  PAROUT(Beta);
  PAROUT(Age);
  PAROUT(Nmu);
  PAROUT(Nx);
  PAROUT(Ny);
  PAROUT(Tc);
  PAROUT(Xc);
  PAROUT(Yc);
  PAROUT(Theta);
  PAROUT(Phi);
  PAROUT(CurvParabA);   // KR adds: curvature free parameters 6/2017
  PAROUT(CurvGaussD);   
  PAROUT(CurvGaussN);
  #undef PAROUT

  #define OUT(name) \
    os << "\n" << std::setw(20) << #name << ": " << name##_
  OUT(typeLDF);
  OUT(typeFrontDelay);
  OUT(typeEnergy);
  OUT(zPos);
  OUT(logLikelihood);
  OUT(logLikelihood_Silent);
  OUT(chi2_LDF);
  OUT(chi2_Time);
  OUT(ndf);
  OUT(nMini);
  OUT(userData);
  #undef OUT

  os << "\n" << std::setw(20) << "energy" << ": " << Energy();

  os << " ]";

  return os;
}

bool I3LaputopParams::operator==(const I3LaputopParams& other) const
{
  return
    ParameterStorageType::operator==(other) &&
    typeLDF_ == other.typeLDF_ &&
    typeFrontDelay_ == other.typeFrontDelay_ &&
    typeEnergy_ == other.typeEnergy_ &&
    nan_aware_equality(zPos_, other.zPos_) &&
    nan_aware_equality(logLikelihood_, other.logLikelihood_) &&
    nan_aware_equality(logLikelihood_Silent_, other.logLikelihood_Silent_) &&
    nan_aware_equality(chi2_LDF_, other.chi2_LDF_) &&
    nan_aware_equality(chi2_Time_, other.chi2_Time_) &&
    ndf_ == other.ndf_ &&
    nMini_ == other.nMini_ &&
    userData_ == other.userData_;
}
