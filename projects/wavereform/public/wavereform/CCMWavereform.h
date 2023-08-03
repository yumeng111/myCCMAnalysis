/**
 * @file 
 * @brief 
 *
 * (c) 2011 the IceCube Collaboration
 *
 * $Id:$
 * @version $Revision$
 * @date $Date$
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 */

#ifndef WAVEREFORM_CCMWAVEREFORM_H_INCLUDED
#define WAVEREFORM_CCMWAVEREFORM_H_INCLUDED

#include <cctype>
#include <string>
#include <vector>
#include <float.h>
#include <utility>
#include <cholmod.h>

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include <icetray/I3Frame.h>
#include <icetray/I3Units.h>
#include <icetray/I3Module.h>
#include <icetray/I3Logging.h>
#include <icetray/CCMPMTKey.h>
#include <icetray/CCMTriggerKey.h>
#include <icetray/I3ConditionalModule.h>

#include <dataclasses/physics/CCMWaveform.h>
#include <dataclasses/physics/NIMLogicPulse.h>
#include <dataclasses/geometry/CCMGeometry.h>
#include <dataclasses/calibration/BaselineEstimate.h>

#include <dataclasses/CCMPMTFunctions.h>
#include <dataclasses/I3TimeWindow.h>
#include <dataclasses/calibration/CCMCalibration.h>
#include <dataclasses/calibration/CCMPMTCalibration.h>
#include <dataclasses/physics/CCMRecoPulse.h>
#include <dataclasses/status/CCMDetectorStatus.h>
#include <dataclasses/status/CCMPMTStatus.h>

class CCMCalibration;
class CCMPMTCalibration;

class CCMWavereform : public I3ConditionalModule {
public:
	CCMWavereform(const I3Context&);
	void Configure();
	void Calibration(I3FramePtr);
	void DAQ(I3FramePtr);
	
    void Geometry(I3FramePtr frame);
    void GetRefoldedWf(std::vector<double> & refolded_wf, std::vector<CCMRecoPulse> const & pulses, CCMPMTCalibration const & calib);
    void GetChi(double & chi, std::vector<double> const & refolded_wf, std::vector<double> const & samples);
	
private:
	std::string waveform_name_, pulse_name_, chi_name_, flag_name_;
    I3Map<CCMPMTKey, uint32_t> pmt_channel_map_;
    std::string geometry_name_;
	double chi_threshold_;
	
	SET_LOGGER("CCMWavereform");
};

#endif /* WAVEREFORM_I3WAVEREFORM_H_INCLUDED */

