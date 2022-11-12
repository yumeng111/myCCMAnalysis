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

#ifndef WAVEREFORM_I3WAVEREFORM_H_INCLUDED
#define WAVEREFORM_I3WAVEREFORM_H_INCLUDED

#include "icetray/I3ConditionalModule.h"
#include "dataclasses/physics/I3Waveform.h"
#include "dataclasses/physics/I3RecoPulse.h"

class I3Calibration;
class I3DetectorStatus;
class I3DOMCalibration;
class I3DOMStatus;

class I3Wavereform : public I3ConditionalModule {
public:
	I3Wavereform(const I3Context&);
	void Configure();
	void Calibration(I3FramePtr);
	void DetectorStatus(I3FramePtr);
	void DAQ(I3FramePtr);
	
	static void RefoldPulses(const std::vector<I3RecoPulse> &pulses,
	    I3Waveform::Source source, unsigned channel, const I3DOMCalibration &calib,
	    const I3DOMStatus &status, double *times, double *bins, size_t n_bins,
	    bool is_simulation = false);
	
	static void RefoldPulses(I3RecoPulseSeries::const_iterator begin,
	    I3RecoPulseSeries::const_iterator end,
	    I3Waveform::Source source, unsigned channel, const I3DOMCalibration &calib,
	    const I3DOMStatus &status, double *times, double *bins, size_t n_bins,
	    bool is_simulation = false);
	
	static unsigned GetChannel(const I3Waveform &wf);
	
	static std::vector<double> GetDigitizerSteps(const I3Waveform &wf,
	    const I3DOMCalibration &calib);
	
	struct Refolded {
		double waveform, refolded, step;
		Refolded(double w, double r, double s) : waveform(w), refolded(r), step(s) {};
	};
	
	static std::vector<Refolded> GetRefolded(const std::vector<I3RecoPulse> &pulses,
	    const I3Waveform &wf, const I3DOMCalibration &calib,
	    const I3DOMStatus &status, bool is_simulation = false);

	// Get the maximum charge of pulses in an integration window of twindow
	static double GetMaxCharge(const std::vector<I3RecoPulse> &pulses,
	    double twindow);
	
private:
	std::string waveform_name_, pulse_name_, chi_name_, flag_name_;
	double chi_threshold_;
	bool use_domsimulator_hacks_;
	boost::shared_ptr<const I3Calibration> calibration_;
	boost::shared_ptr<const I3DetectorStatus> status_;
	
	SET_LOGGER("I3Wavereform");
};

#endif /* WAVEREFORM_I3WAVEREFORM_H_INCLUDED */

