#!/usr/bin/env python3

from I3Tray import *
import os, math, random
from icecube import icetray, dataio, dataclasses
from icecube import phys_services, wavedeform
from os.path import expandvars
from wavedeform_random_waveform_generator import RandomWaveforms

random.seed(0)

tray = I3Tray()

i3_testdata = expandvars("$I3_TESTDATA")

tray.AddModule('I3InfiniteSource', 'reader', Prefix = i3_testdata + '/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz')				
tray.AddModule(RandomWaveforms, 'random', Streams=[icetray.I3Frame.DAQ])
tray.AddModule('I3Wavedeform', 'deform')

import unittest
class TimeCheck(unittest.TestCase):
	fakePulseKey = "RandomPulses"
	pulseKey = "WavedeformPulses"
	
	
	def testCorrespondingPulses(self):
		
		def timeCorrespondence(fakePulse, pulse, margin):
			frac = math.fabs((pulse.time - fakePulse.time)/pulse.width)
			return frac < margin
		
		for om in self.frame[self.fakePulseKey].keys():
			if len(self.frame[self.fakePulseKey][om]) == 0:
				continue
			
			margin = 10
			
			for fakePulse in self.frame[self.fakePulseKey][om]:
				cond = any(x for x in self.frame[self.pulseKey][om] if
						            timeCorrespondence(fakePulse, x, margin))
				self.assert_(cond, ("Reco pulse time within "
								    "%d sigma of fake pulse time" % margin))

	def testLeadingEdgeTime(self):
		fracdiffs = []
		for om in self.frame[self.fakePulseKey].keys():
			if len(self.frame[self.fakePulseKey][om]) == 0:
				continue

			simfirst = self.frame[self.fakePulseKey][om][0].time
			for pulse in self.frame[self.pulseKey][om]:
				recofirst = pulse.time
				recowidth = pulse.width
				# Get a for-sure real pulse
				if pulse.charge > 0.4: break

			margin = 10 # sigma -- this is a crappy test, so it's
			            # big (better test below)
			frac = math.fabs((recofirst - simfirst)/recowidth)
			fracdiffs.append(frac)
			self.assert_(frac < margin, \
			    ("Sim time (%f) and reco time (%f) match to " +
			     "within %d sigma (%f ns) in OM %s") % (simfirst, \
			    recofirst, margin, margin*recowidth, om))

		self.assert_(sum(fracdiffs)/len(fracdiffs) < 1, \
		    "Simulation and reco times within 1 sigma on average")

tray.AddModule(icetray.I3TestModuleFactory(TimeCheck), 'testy', Streams=[icetray.I3Frame.DAQ])

tray.Execute(6)

