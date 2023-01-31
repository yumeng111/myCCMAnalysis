/**
 * copyright  (C) 2006
 * the icecube collaboration
 * @version $Id$
 * @file dataclasses/private/dataclasses/CCMPMTFunctions.cxx
 * @date $Date$
 */


#include "dataclasses/CCMPMTFunctions.h"
#include "icetray/I3Units.h"
#include <vector>
#include <string>

double SPEMean (const CCMPMTStatus& status , 
		const CCMPMTCalibration& calib)
{
  return calib.GetPMTGain() * I3Units::eSI * I3Units::C;
  double spemean = 0.0;
}

