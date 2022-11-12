#!/usr/bin/env python3

from I3Tray import *
import os, math
from icecube import icetray, dataio, dataclasses
from icecube import phys_services, wavedeform
from os.path import expandvars
from wavedeform_random_waveform_generator import RandomWaveforms

tray = I3Tray()

i3_testdata = expandvars("$I3_TESTDATA") 

tray.AddModule('I3InfiniteSource', 'reader', Prefix = i3_testdata + '/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz')
tray.AddModule(RandomWaveforms, 'random', Streams=[icetray.I3Frame.DAQ])
tray.AddModule('I3Wavedeform', 'deform')

import unittest
class NPulsesCheck(unittest.TestCase):
	fakePulseKey = "RandomPulses"
	pulseKey = "WavedeformPulses"

	# Test that pulse merging/splitting are not over-enthusiastic and
	# the number of reco pulses is at least similar to the number of sim
	# pulses. Note that this is not important, nor guaranteed, but tends
	# to break if the algorithm is broken for other reasons.
	def testNPulsesSane(self):
		fracdiffs = []
		for om in self.frame[self.fakePulseKey].keys():
			nsim = len(self.frame[self.fakePulseKey][om])
			nreco = len(self.frame[self.pulseKey][om])

			denom = nreco
			if denom == 0:
				denom = nsim
			if denom == 0:
				frac = 0
			else:
				frac = math.fabs((float(nreco) - nsim)/denom)
			margin = 0.60
			fracdiffs.append(frac)
			condition = (frac < margin or
                                     math.fabs(nreco - nsim) < 1.1)
			self.assert_(condition, \
			    ("Sim npulses (%d) and reco (%d) match to " +
			     "within %d%% in OM %s") % (nsim, \
			    nreco, margin*100, om))

		print(sum(fracdiffs)/len(fracdiffs))
		self.assert_(sum(fracdiffs)/len(fracdiffs) < 0.15, \
		    "Simulation and reco N pulses within 15% on average")

tray.AddModule(icetray.I3TestModuleFactory(NPulsesCheck), 'testy', Streams=[icetray.I3Frame.DAQ])

tray.Execute(6)

