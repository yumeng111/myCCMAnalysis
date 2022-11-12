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

#include "wavereform/I3Wavereform.h"
#include "dataclasses/I3DOMFunctions.h"
#include "dataclasses/physics/I3RecoPulse.h"
#include "dataclasses/physics/I3Waveform.h"
#include "dataclasses/status/I3DOMStatus.h"
#include "dataclasses/calibration/I3DOMCalibration.h"

#include "icetray/I3Units.h"

#include <boost/foreach.hpp>

/*
 * A little template abuse: pulse template functions that know about their
 * bounds of support. Even though the template functions technically extend
 * over the entire real line, we truncate them at reasonably small values
 * for efficiency's sake.
 */

struct ATWDPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		return calib.ATWDPulseTemplate()(t);
	}
	static double t_min, t_max;
};

double ATWDPulse::t_min = -10*I3Units::ns; /* 1e-8 */
double ATWDPulse::t_max = 100*I3Units::ns; /* 1e-8 */

struct ATWDChannel0Pulse : public ATWDPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		return calib.ATWDPulseTemplate(0)(t);
	}
};

struct ATWDChannel1Pulse : public ATWDPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		return calib.ATWDPulseTemplate(1)(t);
	}
};

struct ATWDChannel2Pulse : public ATWDPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		return calib.ATWDPulseTemplate(2)(t);
	}
};

struct FADCPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		return calib.FADCPulseTemplate()(t);
	}
	static double t_min, t_max;
};

double FADCPulse::t_min = -50*I3Units::ns; /* 8e-10 */
double FADCPulse::t_max = 500*I3Units::ns; /* 4e-10 */

/* 
 * XXX: These are also present in wavedeform. If they are changed here,
 * they also need to be changed in wavedeform (or vice versa)
 */

namespace {

struct SPETemplate { double c, x0, b1, b2; };

const SPETemplate DOMsimulatorATWDNewToroid = {
	1.07138375,
	2.47418736,
	5.06631611,
	42.18168102
};

const SPETemplate DOMsimulatorATWDOldToroid = {
	2.11429879,
	4.2001904,
	5.4693169,
	32.60217968
};

const SPETemplate DOMsimulatorFADC = {
	1.75632114e-01,
	16.3087924,
	2.17527146e+01,
	2.13038550e+02
};

inline double
SPEPulseShape(double t, const SPETemplate &p)
{
	return p.c*pow(exp(-(t - p.x0)/p.b1) + exp((t - p.x0)/p.b2),-8);
}

struct DOMsimulatorATWDPulse : ATWDPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		if (calib.GetToroidType()==I3DOMCalibration::NEW_TOROID)
			return SPEPulseShape(t, DOMsimulatorATWDNewToroid);
		else
			return SPEPulseShape(t, DOMsimulatorATWDOldToroid);
	}
};

struct DOMsimulatorFADCPulse : FADCPulse {
	inline double operator()(const I3DOMCalibration &calib, double t)
	{
		return SPEPulseShape(t, DOMsimulatorFADC);
	}
};

} // namespace

template <typename PulseTemplate>
static void
RefoldPulsesImpl(I3RecoPulseSeries::const_iterator pulse, I3RecoPulseSeries::const_iterator pend,
    const I3DOMCalibration &calib, const I3DOMStatus &status,
    double *times, double *bins, size_t n_bins)
{
	const double spe_charge =
	    SPEMean(status, calib)*calib.GetFrontEndImpedance();
	double *t_start = times;
	double *t_end = times + 1;
	PulseTemplate templ;
	
	for ( ; pulse != pend; pulse++) {
		/* Find the bounds of support for this pulse template */
		t_start = std::lower_bound(t_start, times + n_bins,
		    pulse->GetTime() + PulseTemplate::t_min);
		t_end = std::upper_bound(t_end, times + n_bins,
		    pulse->GetTime() + PulseTemplate::t_max);
		
		/* Compute the corresponding indices */
		const unsigned begin = t_start-times;
		const unsigned end = t_end-times;
		
		/* Convolve! */
		const double norm = spe_charge*pulse->GetCharge();
		for (unsigned i = begin; i < end; i++)
			bins[i] += norm*templ(calib, times[i]-pulse->GetTime());
	}
}

void
I3Wavereform::RefoldPulses(const I3RecoPulseSeries &pulses,
    I3Waveform::Source source, unsigned channel, const I3DOMCalibration &calib,
    const I3DOMStatus &status, double *times, double *bins, size_t n_bins, bool is_simulation)
{
	RefoldPulses(pulses.begin(), pulses.end(), source, channel, calib, status,
	    times, bins, n_bins, is_simulation);
}

void
I3Wavereform::RefoldPulses(I3RecoPulseSeries::const_iterator begin, I3RecoPulseSeries::const_iterator end,
    I3Waveform::Source source, unsigned channel, const I3DOMCalibration &calib,
    const I3DOMStatus &status, double *times, double *bins, size_t n_bins, bool is_simulation)
{	
	memset(bins, 0, n_bins*sizeof(double));
	
	switch (source) {
		case I3Waveform::ATWD:
			if (is_simulation)
				RefoldPulsesImpl<DOMsimulatorATWDPulse>(begin, end, calib, status, times, bins, n_bins);
			else
				switch (channel) {
					case 0:
						RefoldPulsesImpl<ATWDChannel0Pulse>(begin, end, calib, status, times, bins, n_bins);
						break;
					case 1:
						RefoldPulsesImpl<ATWDChannel1Pulse>(begin, end, calib, status, times, bins, n_bins);
						break;
					case 2:
						RefoldPulsesImpl<ATWDChannel2Pulse>(begin, end, calib, status, times, bins, n_bins);
						break;
					default:
						log_fatal("Unknown ATWD channel %u!", channel);
				}
			break;
		case I3Waveform::FADC:
			if (is_simulation)
				RefoldPulsesImpl<DOMsimulatorFADCPulse>(begin, end, calib, status, times, bins, n_bins);
			else
				RefoldPulsesImpl<FADCPulse>(begin, end, calib, status, times, bins, n_bins);
			break;
		default:
			log_fatal("Unkown waveform source!");
	}
}

std::vector<double>
I3Wavereform::GetDigitizerSteps(const I3Waveform &waveform,
    const I3DOMCalibration &calib)
{
	const size_t nbins = waveform.GetWaveform().size();
	std::vector<double> steps(nbins, 0.5*I3Units::mV);
	const std::vector<I3Waveform::StatusCompound> &wf_info
	    = waveform.GetWaveformInformation();
	
	if (waveform.GetDigitizer() == I3Waveform::ATWD) {
		unsigned atwd_id = waveform.GetSourceIndex();		
		std::vector<I3Waveform::StatusCompound>::const_iterator
		    status_it = wf_info.begin();
		
		unsigned i = 0;
		int channel = 0;
		
		while (i < nbins) {
			unsigned end = 0;
			
			while (end == 0 && status_it != wf_info.end()) {
				if (i < status_it->GetInterval().first) {
					end = status_it->GetInterval().first;
					channel = 0;
				} else if (i < status_it->GetInterval().second) {
					end = status_it->GetInterval().second;
					channel = status_it->GetChannel();
					if (channel < 0) {
						log_warn("Negative channel number found in I3Waveform status! "
						    "The step size in this range is meaningless.");
						channel = 0;
					}
				} else
					status_it++;
			}
			
			if (status_it == wf_info.end()) {
				end = nbins;
				channel = 0;
			}
			
			for ( ; i < end; i++)
				steps[i] = calib.GetATWDBinCalibSlope(atwd_id, channel, i)/calib.GetATWDGain(channel);
		}
	} else if (waveform.GetDigitizer() == I3Waveform::FADC) {
		steps = std::vector<double>(nbins, calib.GetFADCGain());
	} else {
		log_error("I can only handle ATWD and FADC waveforms.");
	}
	
	return steps;
}

inline bool
IsBorked(const I3Waveform &wf, const I3Waveform::StatusCompound &status)
{
	const I3Waveform::Status stat = status.GetStatus();
	switch (wf.GetSource()) {
		case I3Waveform::ATWD:
			return ((stat == I3Waveform::SATURATED) || (stat == I3Waveform::UNDERSHOT));
			break;
		case I3Waveform::FADC:
			return (stat != I3Waveform::VIRGINAL);
			break;
		default:
			return false;
	}
}

unsigned
I3Wavereform::GetChannel(const I3Waveform &wf)
{
	unsigned channel = 0;
	const std::vector<I3Waveform::StatusCompound> &winfo = wf.GetWaveformInformation();
	std::vector<I3Waveform::StatusCompound>::const_iterator it = winfo.begin();
	if (it != winfo.end())
		channel = it->GetChannel();
	else
		return channel;
		
	it++;
	for ( ; it != winfo.end(); it++)
		if (channel != unsigned(it->GetChannel()))
			log_error("Mixed channels %u and %d found while refolding ATWD waveform! "
			    "Defaulting to highest-gain channel.", channel, it->GetChannel());
			
	return channel;	
}

std::vector<I3Wavereform::Refolded>
I3Wavereform::GetRefolded(const std::vector<I3RecoPulse> &pulses,
    const I3Waveform &waveform, const I3DOMCalibration &calib,
    const I3DOMStatus &status, bool is_simulation)
{
	unsigned i;
	const double t0 = waveform.GetStartTime();
	const double dt = waveform.GetBinWidth();
	const size_t n_bins = waveform.GetWaveform().size();
	const std::vector<double> &wf = waveform.GetWaveform();
	std::vector<double> edge_times(n_bins, 0.0);
	std::vector<double> refolded(n_bins, 0.0);
	std::vector<double> steps = GetDigitizerSteps(waveform, calib);
	std::vector<Refolded> result;
	
	for (i = 0; i < n_bins; i++)
		edge_times[i] = t0 + (i+1)*dt;
		
	RefoldPulses(pulses, waveform.GetDigitizer(), waveform.GetChannel(), calib, status,
	    &edge_times.front(), &refolded.front(), n_bins, is_simulation);
	
	std::vector<I3Waveform::StatusCompound>::const_iterator
		status_it = waveform.GetWaveformInformation().begin();
	
	for (i = 0; i < n_bins; i++) {
		/* Ignore saturated or undershot regions. */
		if (status_it != waveform.GetWaveformInformation().end()) {
			if (IsBorked(waveform, *status_it) &&
			    i >= status_it->GetInterval().first) {
				if (status_it->GetInterval().second == n_bins)
					break;
				i = status_it->GetInterval().second;
				status_it++;
			} else if (i >= status_it->GetInterval().second)
				status_it++;
		}
		
		result.push_back(Refolded(wf[i], refolded[i], steps[i]));
	}
		
	return result;
}

typedef I3RecoPulseSeries::const_iterator pulse_it_t;
typedef std::pair<pulse_it_t, pulse_it_t> pulse_it_pair_t;

static
std::vector<pulse_it_pair_t>
GenerateTimeWindows(const I3RecoPulseSeries &pulses, const double twindow)
{
	std::vector<pulse_it_pair_t> ranges;

	pulse_it_t begin = pulses.begin(), end = begin+1;
	if (twindow > 0)
		for ( ; begin < pulses.end(); begin++, end = begin+1) {
			while (end < pulses.end() &&
			    end->GetTime() - begin->GetTime() < twindow)
				end++;
			ranges.push_back(std::make_pair(begin, end));
		}
	else
		ranges.push_back(std::make_pair(pulses.begin(), pulses.end()));

	return ranges;
}

double
I3Wavereform::GetMaxCharge(const I3RecoPulseSeries &pulses, double twindow)
{
	double qtot_max = 0;

	BOOST_FOREACH(const pulse_it_pair_t &range,
	    GenerateTimeWindows(pulses, twindow)) {
		double qtot = 0;
		for (pulse_it_t p = range.first; p != range.second; p++)
			qtot += p->GetCharge();

		if (qtot > qtot_max)
			qtot_max = qtot;
	}

	return qtot_max;
}

