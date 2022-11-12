#!/usr/bin/env python3

from I3Tray import *
import os, math
from icecube import icetray, dataio, dataclasses
from icecube import phys_services, wavedeform
from os.path import expandvars
from wavedeform_random_waveform_generator import RandomWaveforms

i3_testdata = expandvars("$I3_TESTDATA")

tray = I3Tray()
tray.AddModule('I3InfiniteSource', 'reader', Prefix = i3_testdata + '/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz')				
tray.AddModule(RandomWaveforms, 'random', Streams=[icetray.I3Frame.DAQ])
tray.AddModule('I3Wavedeform', 'deform')

import unittest
class SanityCheck(unittest.TestCase):
	fakePulseKey = "RandomPulses"
	pulseKey = "WavedeformPulses"

	def testKeys(self):
		self.assert_(self.pulseKey in self.frame, "The pulses actually show up in the frame.")
	def testPositiveWidth(self):
		psm = self.frame[self.pulseKey]
		for om, vec in psm:
			for p in vec:
				self.assert_(p.width >= 0, "Pulse width (%g) is non-negative" % p.width)
	def testPulseSorting(self):
		psm = self.frame[self.pulseKey]
		for om, pulses in psm:
			for prev, next in zip(pulses[:-1], pulses[1:]):
				self.assert_(next.time - prev.time > 0, "Pulses at t=%g and %g are strictly ordered in time" % (prev.time, next.time))
	def testTimeRange(self):
		range = self.frame[self.pulseKey + 'TimeRange']
		sourcerange = self.frame['CalibratedWaveformRange']
		self.assert_(range.stop == sourcerange.stop, "Time ranges stop at the same time (%f, %f)" % (range.stop, sourcerange.stop))
		self.assert_(range.start == sourcerange.start - 25*I3Units.ns, "Time range starts 25 ns before waveforms (%f, %f)" % (range.start, sourcerange.start))

tray.AddModule(icetray.I3TestModuleFactory(SanityCheck), 'testy', Streams=[icetray.I3Frame.DAQ])

tray.Execute(6)

