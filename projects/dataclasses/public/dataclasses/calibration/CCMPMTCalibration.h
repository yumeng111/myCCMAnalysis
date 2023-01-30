/**
 *
 * Definition of CCMPMTCalibration class
 *
 * copyright  (C) 2004
 * the IceCube collaboration
 * @version $Id$
 * @file CCMPMTCalibration.h
 * @date $Date$
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


static const unsigned ccmpmtcalibration_version_ = 0;

/**
 * @brief Class that stores the calibration information for a DOM
 * 
 * This class stores the information from the Calibration stream.
 *
 * ATTENTION:
 * Calibration information is set assuming the bin number as it is in
 * the calibration database (reversed in time for bins 0-127).
 * Calibration information is fetched by the time-ordered bin numbers.
 *
 * Nothing here really has units.  Most are fits/offsets that
 * convert values for later storage (at that SHOULD use I3Units)
 *
 * @author Tom McCauley <tpmccauley@lbl.gov>
 * @author Erik Blaufuss \<blaufuss at icecube umd edu\>
 *
 * 
 */

class CCMPMTCalibration {
  
 public:
  CCMPMTCalibration();
  ~CCMPMTCalibration();
  
  /**
   * Get MB Temperature at time of calibration
   */
  double GetTemperature() const { return temperature_; }
  
  /**
   * Set MB Temperature at time of calibration
   */
  void SetTemperature(double temperature);
    

  /**
   * Get PMT Gain- relation between measured counts and mV.
   */
  double GetPMTGain() const { return pmtGain_; }

  /**
   * Get PMT Pedestal- baseline point from which waveforms start.
   */
  LinearFit GetPMTBaselineFit() const { return pmtBaselineFit_ ; }

  /**
   *  Get the PMT intrinsic time offset (in units of time)
   */
  double GetPMTDeltaT() const { return pmtDeltaT_; }
  
  /**
   * Get parameters for droop correction on the baseline
   */
  TauParam GetTauParameters() const { return tauparameters_; }
  
  /**
   * Set parameters for droop correction on the baseline
   */
  void SetTauParameters(TauParam tauparameters);



  /**
   * Set PMT calibration parameters. Currently the FADC
   * calibration is a work in progress and a moving target
   * so this is only a tentative interface -tpm
   */
  void SetPMTGain(double gain)
  {
    pmtGain_     = gain;
  };

  void SetPMTBaselineFit(LinearFit basefit)
    {
      pmtBaselineFit_ = basefit;
    }
    
  void SetPMTDeltaT(double deltaT)
    {
      pmtDeltaT_ = deltaT;
    }
    


  /**
   *  Get/set for relativePMTEff
   */
  double GetRelativePMTEff() const { return relativePMTEff_ ; }
  void SetRelativePMTEff(double relaeff)
  {
    relativePMTEff_ = relaeff;
  }

  //On the assumption that this will be evaulated many times, we copy all data into it. 
  //This makes the object larger but hopefully avoids extra pointer derefences
  class DroopedSPETemplate{
  public:
    SPETemplate pulse;
    struct droopParams{
      SPETemplate pulse;
      double tauFrac, time1, time2;
      
      droopParams(){}
      droopParams(const SPETemplate& templ, 
            double tauFrac, double time1, double time2):
      pulse(templ),tauFrac(tauFrac),time1(time1),time2(time2){}
    } droop;
    bool droopy;
    
    DroopedSPETemplate(const SPETemplate& templ):
    pulse(templ),droopy(false){}
    
    DroopedSPETemplate(const SPETemplate& templ,
               const SPETemplate& droopTempl, 
               double tauFrac, double time1, double time2):
    pulse(templ),droop(droopTempl,tauFrac,time1,time2),droopy(true){}
    
    double operator()(double t){
      if (!droopy)
        return SPEPulseShape(t);
      
      double norm = (1.0 - droop.tauFrac)*droop.time1 + droop.tauFrac*droop.time2;
      double c1 = (1.0 - droop.tauFrac)*droop.time1/norm;
      double c2 = droop.tauFrac*droop.time2/norm;
      
      return SPEPulseShape(t) +
      c1*DroopReactionShape(t, droop.time1) +
      c2*DroopReactionShape(t, droop.time2);
    }
    
    bool operator==(const DroopedSPETemplate& templ) const;
    bool operator!=(const DroopedSPETemplate& templ) const;
    bool operator<(const DroopedSPETemplate& templ) const;
    
  private:
    double SPEPulseShape(double t) const {
      return pulse.c/std::pow(exp(-(t - pulse.x0)/pulse.b1) + 
                         exp((t - pulse.x0)/pulse.b2),8);
    }
    
    double DroopReactionShape(double t, double tau) const{
      return (pulse.c*droop.pulse.c/tau)/
        std::pow(exp(-(t - pulse.x0*droop.pulse.x0)/(pulse.b1*droop.pulse.b1)) + 
            exp((t - pulse.x0*droop.pulse.x0)/(pulse.b2*droop.pulse.b2*tau)),8);
    }
  };

  DroopedSPETemplate PMTPulseTemplate(bool droopy = false) const;
 
  template <class Archive>
    void serialize(Archive& ar, unsigned version);
    
  
  /**
   *  PMT-specific corrections to the SPE charge distribution
   */
  double GetMeanPMTCharge() const {return meanPMTCharge_;}

  /**
   * In dataclasses we use NaN to denote "invalid" however in the JSON
   * file invalid entries are set to 0.  To cover both cases we check 
   * that it's finite and greater than 0.
   */ 
  bool IsMeanPMTChargeValid() const {
    return ((std::isfinite(meanPMTCharge_)) && (meanPMTCharge_ > 0.));
  }

  void SetMeanPMTCharge(double charge) {
    meanPMTCharge_ = charge;
  }

  const SPEChargeDistribution& GetCombinedSPEChargeDistribution() const {
    return combinedSPEFit_;
  }

  void SetCombinedSPEChargeDistribution(const SPEChargeDistribution& fit) {
    combinedSPEFit_ = fit;
  }
  
  bool operator==(const CCMPMTCalibration& rhs) const
  {    
    return (CompareFloatingPoint::Compare_NanEqual(droopTimeConstants_[0],rhs.droopTimeConstants_[0]) &&
        CompareFloatingPoint::Compare_NanEqual(droopTimeConstants_[1],rhs.droopTimeConstants_[1]) &&
        CompareFloatingPoint::Compare_NanEqual(temperature_,rhs.temperature_) &&
        CompareFloatingPoint::Compare_NanEqual(pmtGain_,rhs.pmtGain_) &&
        pmtBaselineFit_ == rhs.pmtBaselineFit_ &&
        CompareFloatingPoint::Compare_NanEqual(pmtDeltaT_,rhs.pmtDeltaT_) &&
        tauparameters_ == rhs.tauparameters_ &&
        CompareFloatingPoint::Compare_NanEqual(relativePMTEff_,rhs.relativePMTEff_) &&
        combinedSPEFit_ == rhs.combinedSPEFit_ &&
        CompareFloatingPoint::Compare_NanEqual(meanPMTCharge_,rhs.meanPMTCharge_));
  }
  bool operator!=(const CCMPMTCalibration& rhs) const
  {
    return !operator==(rhs);
  }

 private:

  double SPEPulseTemplate(double t, const SPETemplate& templ,
    const SPETemplate& droop, bool droopy) const;

  double droopTimeConstants_[2];

  double  temperature_;
 
  /**
   * Gain and pedestal values
   */
  double pmtGain_;
  LinearFit pmtBaselineFit_;
  
  /**
   *  Inherent time offset (ns)
   */
  double pmtDeltaT_;
  
  /**
   *   Parameters for droop correction   
   */
  TauParam tauparameters_;
  
  /**
   *  Relative PMT efficiency, normalized to 1.0 for the average pmt.
   */
  double relativePMTEff_;

  /**
   *  Combined-fit SPE distribution function (exponential1 + exponential2 + Gaussian)
   */
  SPEChargeDistribution combinedSPEFit_;

  /**
   *  PMT-specific corrections to the SPE charge distribution
   */
  double meanPMTCharge_;

  /**
   *  Allow the Diff compression class to directly use private data
   */
  friend class CCMPMTCalibrationDiff;
};

typedef std::map<CCMPMTKey, CCMPMTCalibration> CCMPMTCalibrationMap;
I3_POINTER_TYPEDEFS(CCMPMTCalibrationMap);

I3_CLASS_VERSION(CCMPMTCalibration, ccmpmtcalibration_version_);
I3_POINTER_TYPEDEFS(CCMPMTCalibration);

std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibration& c);
std::ostream& operator<<(std::ostream& oss, const CCMPMTCalibrationMap& m);

#endif //CCMPMTCalibration_H_INCLUDED
