/**
 * copyright  (C) 2004
 * the icecube collaboration
 * @version $Id$
 * @file CCMRecoPulse.h
 * @date $Date$
 */

#ifndef CCMRECOPULSE_H_INCLUDED
#define CCMRECOPULSE_H_INCLUDED

#include "dataclasses/Utility.h"
#include "dataclasses/I3Vector.h"
#include "icetray/CCMPMTKey.h"
#include "dataclasses/I3Map.h"
#include "icetray/I3Frame.h"

/**
 * @brief A storage class for extracted pulses from a feature extractor
 * A readout independent representation of a waveform feature or Analog
 *  readout.
 */
static const unsigned ccmrecopulse_version_ = 1;

class CCMRecoPulse 
{
  double time_;
  float charge_;
  float width_;

  public:

  /**
   * Construct a pulse with all properties initialized to Not-a-Number.
   */
  CCMRecoPulse() : time_(NAN), charge_(NAN), width_(NAN) {}

  /**
   * Print a string representation of this pulse
   */
  std::ostream& Print(std::ostream&) const;
  
  /**
   * @brief Get the start time of the pulse.
   */
  double GetTime() const {return time_;}

  /**
   * @brief Set the start time of the pulse.
   */
  void SetTime(double time) {time_ = time;}

  /**
   * @brief Get the amplitude of the pulse in units of ideally amplified
   *        photoelectrons.
   */
  float GetCharge() const {return charge_;}

  /**
   * @brief Set the amplitude of the pulse in units of ideally amplified
   *        photoelectrons.
   */
  void SetCharge(float charge) {charge_ = charge;}

  /**
   * @brief The time between this pulse and the subsequent basis function
   *        used in the pulse unfolding.
   *
   * This quantity can be used to approximate the time interval within which the
   * pulse charged was observed as [ GetTime() , GetTime() + GetWidth() ]. This
   * approximation becomes dubious, however, when the charge density is high
   * (such that several adjacent basis functions in the unfolding were assigned
   * nonzero amplitudes).
   */
  float GetWidth() const {return width_;}

  /**
   * @brief Store the time between this pulse and the subsequent basis function
   *        used in the pulse unfolding.
   */
  void SetWidth(float width) {width_ = width;}

  ~CCMRecoPulse();

  bool operator==(const CCMRecoPulse& rhs) const;
  bool operator!=(const CCMRecoPulse& rhs) const { return !operator==(rhs); }

  private:
  friend class icecube::serialization::access;
  template <class Archive> void serialize(Archive & ar, unsigned version);
};

I3_POINTER_TYPEDEFS(CCMRecoPulse);
I3_CLASS_VERSION(CCMRecoPulse, ccmrecopulse_version_);

typedef std::vector<CCMRecoPulse> CCMRecoPulseSeries;
typedef I3Map<CCMPMTKey, CCMRecoPulseSeries> CCMRecoPulseSeriesMap;
typedef I3Map<CCMPMTKey, CCMRecoPulse> CCMRecoPulseMap;

std::ostream& operator<<(std::ostream& oss, const CCMRecoPulse& p);

I3_POINTER_TYPEDEFS(CCMRecoPulseSeries);
I3_POINTER_TYPEDEFS(CCMRecoPulseSeriesMap);
I3_POINTER_TYPEDEFS(CCMRecoPulseMap);

/*
 * Specialize I3Frame::Get() to turn convert various objects
 * in the frame into CCMRecoPulseSeriesMaps.
 */

// need to hide this from ROOT
#ifndef __CINT__
#include "icetray/I3Frame.h"

template <>
CCMRecoPulseSeriesMapConstPtr
I3Frame::Get(const std::string& name, void*, void*) const;
#endif //__CINT__

#endif //CCMRECOPULSE_H_INCLUDED


