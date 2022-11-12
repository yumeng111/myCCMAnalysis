/**
 * copyright (C) 2006
 * $Id: I3CscdLlhFitParams.cxx 41418 2008-01-25 11:20:12Z schlenst $
 *
 * @file I3CscdLlhFitParams.cxx
 * @version:
 * @date 17Apr2006
 * @author Doug Rutledge
 */

#include "icetray/serialization.h"
#include "recclasses/I3CscdLlhFitParams.h"

I3CscdLlhFitParams::~I3CscdLlhFitParams() {}

std::ostream& I3CscdLlhFitParams::Print(std::ostream& os) const{
  os << "[I3CscdLlhFitParams HitCount : " << hitCount_ << '\n'
     << "                  HitOMCount : " << hitOmCount_ << '\n'
     << "                UnhitOMCount : " << unhitOmCount_ << '\n'
     << "                   FitStatus : " << status_ << '\n'
     << "                   TimeError : " << errT_ << '\n'
     << "                      XError : " << errX_ << '\n'
     << "                      YError : " << errY_ << '\n'
     << "                      ZError : " << errZ_ << '\n'
     << "                 ZenithError : " << errTheta_ << '\n'
     << "                AzimuthError : " << errPhi_ << '\n'
     << "                 EnergyError : " << errEnergy_ << '\n'
     << "                 NegativeLlh : " << negLlh_ << '\n'
     << "                  ReducedLlh : " << reducedLlh_ << ']';
  return(os);
}

std::ostream& operator<<(std::ostream& os, const I3CscdLlhFitParams& fp){
  return(fp.Print(os));
}

template <class Archive>
void I3CscdLlhFitParams::serialize(Archive& archive, unsigned version)
{
  archive & make_nvp("I3FrameObject",base_object<I3FrameObject>(*this));
  archive & make_nvp("HitCount",hitCount_);
  archive & make_nvp("HitOMCount",hitOmCount_);
  archive & make_nvp("UnhitOMCount",unhitOmCount_);
  archive & make_nvp("FitStatus",status_);
  archive & make_nvp("TimeError",errT_);
  archive & make_nvp("XError",errX_);
  archive & make_nvp("YError",errY_);
  archive & make_nvp("ZError",errZ_);
  archive & make_nvp("ZenithError",errTheta_);
  archive & make_nvp("AzimuthError",errPhi_);
  archive & make_nvp("EnergyError",errEnergy_);
  archive & make_nvp("NegativeLlh",negLlh_);
  archive & make_nvp("ReducedLlh",reducedLlh_);
  //archive & make_nvp("",);
}

I3_SERIALIZABLE(I3CscdLlhFitParams);
