/**
 *
 * Definition of CCMCalibration class
 *
 * copyright  (C) 2004
 * the IceCube collaboration
 * @version $Id$
 * @file CCMCalibration.h
 * @date $Date$
 */

#ifndef CCMCalibration_H_INCLUDED
#define CCMCalibration_H_INCLUDED

#include <map>

#include "dataclasses/Utility.h"
#include "dataclasses/calibration/CCMPMTCalibration.h"
#include "dataclasses/I3Time.h"
#include "icetray/OMKey.h"
#include <icetray/I3FrameObject.h>
#include <icetray/I3DefaultName.h>

static const unsigned ccmcalibration_version_ = 1;

class CCMCalibration : public I3FrameObject {
public:
  I3Time startTime;
  I3Time endTime;

  CCMCalibration();
    
  ~CCMCalibration();
    
  CCMPMTCalibrationMap pmtCal;
  
  bool operator==(const CCMCalibration& rhs)
  {
    return (startTime == rhs.startTime &&
            endTime == rhs.endTime &&
            pmtCal == rhs.pmtCal);
  }
  bool operator!=(const CCMCalibration& rhs)
  {
    return !operator==(rhs);
  }
  
private:
  friend class icecube::serialization::access;
  template <class Archive> void load(Archive & ar, unsigned version);
  template <class Archive> void save(Archive & ar, unsigned version) const;
  I3_SERIALIZATION_SPLIT_MEMBER();
};

I3_CLASS_VERSION(CCMCalibration, ccmcalibration_version_);
I3_DEFAULT_NAME(CCMCalibration);
I3_POINTER_TYPEDEFS(CCMCalibration);

#endif // CCMCalibration_H_INCLUDED
    

