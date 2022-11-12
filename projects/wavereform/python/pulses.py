"""
Python implementations of the pulse templates, for testing and plotting.
"""

import numpy
from icecube.icetray import I3Units

params = {
	'fadc' : [25.12 / 71.363940160184669, 61.27 - 50, 30., 186.],
	'atwd' : {
		'old' : [
			[15.47 / 13.292860653948139, -3.929 - 5, 4.7, 39.],
			[2.07399312, -10.95781298, 4.86019733, 30.74826947],
			[1.35835821, -9.68624195,  3.5016398,  30.96897853]
		],
		'new' : [
			[17.899 / 14.970753076313095, -4.24 - 5, 5.5, 42.],
			[1.6581978,  -11.70227755, 5.4664884,  36.22319705],
			[0.70944364, -10.58782492, 3.48330553, 42.10873959]
		]
	},
	
	'domsimulator' : {
		'fadc' : [1.75632114e-01, 16.3087924, 2.17527146e+01, 2.13038550e+02],		
		'atwd' : {
			'old' : [[2.11429879, 4.2001904, 5.4693169, 32.60217968]]*3,
			'new' : [[1.07138375, 2.47418736, 5.06631611, 42.18168102]]*3,
		}
	}
}

def pulse_template(time, args):
	c, x0, b1, b2 = args
	t = time-11.5
	func = c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8
	return func

def spe_fadc(time):
	t = time - 11.5 # causality
	c = 25.12 / 71.363940160184669; x0 = 61.27 - 50; b1 = 30.;  b2 = 186.
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8

def spe_atwd_new(time):
	t = time - 11.5 # causality
	c = 17.899 / 14.970753076313095; x0 = -4.24 - 5; b1 = 5.5; b2 = 42.;
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8
	
def spe_atwd1_new(time):
	t = time - 11.5 # causality
	c = 1.6581978;  x0 = -11.70227755; b1 = 5.4664884;  b2 = 36.22319705;
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8
	
def spe_atwd2_new(time):
	t = time - 11.5 # causality
	c = 0.70944364;  x0 = -10.58782492; b1 = 3.48330553;  b2 = 42.10873959;
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8

def spe_atwd_old(time):
	t = time - 11.5 # causality
	c = 15.47 / 13.292860653948139; x0 = -3.929 - 5; b1 = 4.7; b2 = 39.;
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8
	
def spe_atwd1_old(time):
	t = time - 11.5 # causality
	c = 2.07399312;  x0 = -10.95781298; b1 = 4.86019733;  b2 = 30.74826947;
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8
	
def spe_atwd2_old(time):
	t = time - 11.5 # causality
	c = 1.35835821;  x0 = -9.68624195;  b1 = 3.5016398;   b2 = 30.96897853;
	return c*(numpy.exp(-(t - x0)/b1) + numpy.exp((t - x0)/b2))**-8

def wfdata(wf, title=None):
	"""Return waveform in mV"""
	if len(wf.waveform) == 0:
		return numpy.zeros(1), numpy.zeros(1)
	#x = wf.start_time + numpy.linspace(1, wf.bin_width*(len(wf.waveform)), len(wf.waveform))
	x = wf.time + wf.bin_width*numpy.arange(1, len(wf.waveform)+1, dtype=float)
	y = numpy.asarray(wf.waveform) # no unit conversions.
	return x, y

def refold_pulses(pulses, template, times, norm):
	"""Convolve pulses with SPE shape."""
	wf = numpy.zeros(len(times))
	for pulse in pulses:
		if pulse.time > (times[-1]+100*I3Units.ns) or pulse.time < (times[0] - 1*I3Units.ms):
			continue
		wf += pulse.charge* \
		     template(times - pulse.time)
	wf *= norm
	return wf
	
def pulseseries_to_wf(pulses, atwd_width, times, norm, title=None):
	"""Convolve pulses with SPE shape."""
	wf = numpy.zeros(len(times))
	for pulse in pulses:
		if pulse.time > (times[-1]+100*I3Units.ns) or pulse.time < (times[0] - 1*I3Units.ms):
			continue
		if atwd_width > 40:
			wf += pulse.charge* \
			     spe_atwd_new(times - pulse.time)
		elif atwd_width > 0:
			wf += pulse.charge* \
			     spe_atwd_old(times - pulse.time)
		else:
			wf += pulse.charge* \
			     spe_fadc(times - pulse.time)

	wf *= norm
	return wf

def droop(pulse, binSize, timeConstant1, timeConstant2, factor):
	"""Undo droop correction; braindead copy-paste from DOMsimulator"""
	
	Tau1OverSampleWidth = timeConstant1/binSize
	Tau2OverSampleWidth = timeConstant2/binSize

	Exp1 = numpy.exp(-1./Tau1OverSampleWidth);
	Exp2 = numpy.exp(-1./Tau2OverSampleWidth);

	A1 = Tau1OverSampleWidth*(1.- Exp1);
	A2 = Tau2OverSampleWidth*(1.- Exp2);

	C1 = (1. - factor) * Tau1OverSampleWidth;
	C2 =        factor * Tau2OverSampleWidth;

	w1 = C1/(C1+C2);
	w2 = C2/(C1+C2);

	coef0 = 1./(w1*A1 + w2*A2);
	coef1 = w1*A1*A1/Tau1OverSampleWidth;
	coef2 = w2*A2*A2/Tau2OverSampleWidth;
	
	S1 = 0.;
	S2 = 0.;

	X = 0.

	drooped = pulse.copy()

	for i in range(1, pulse.size):
		S1 = X + Exp1*S1;
		S2 = X + Exp2*S2;

		drooped[i] = pulse[i]/coef0 - coef1*S1 - coef2*S2
		X = pulse[i]

	return drooped
