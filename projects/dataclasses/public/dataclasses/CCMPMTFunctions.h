/**
 * copyright  (C) 2006
 * the icecube collaboration
 * @version $Id$
 * @file I3DOMFunctions.h
 * @date $Date$
 */
#ifndef CCMPMTFunctions_H_INCLUDED
#define CCMPMTFunctions_H_INCLUDED

//includes
#include "dataclasses/status/CCMPMTStatus.h"
#include "dataclasses/calibration/CCMPMTCalibration.h"
#include "dataclasses/physics/I3DOMLaunch.h"

#include <vector>

/**
 * Return the calculated SPEMean (the ADC "charge", in GV ns, corresponding to
 * the calibrated value of 1 PE). This converts ADC amplitudes into units of
 * photoelectrons.
 */
double SPEMean (const CCMPMTStatus& , const CCMPMTCalibration&);

#endif //CCMPMTFunctions_H_INCLUDED
