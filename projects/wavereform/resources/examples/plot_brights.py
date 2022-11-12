#!/usr/bin/env python3

from optparse import OptionParser, Option, OptionValueError
from icecube import icetray

def check_omkey(option, opt, value):
	import re
	match = re.match('\(?(\d+),(\d+)\)?', value)
	if match:
		return icetray.OMKey(int(match.groups()[0]), int(match.groups()[1]))
	else:
		raise OptionValueError("Can't parse a string/om pair out of '%s'. Use e.g. --key=1,2 or --key=(1,2) to specify String 1, OM 2" % value)
	
Option.TYPES = Option.TYPES + ("omkey",)
Option.TYPE_CHECKER["omkey"] = check_omkey

parser = OptionParser()
parser.add_option('--waveforms', dest='waveforms', type=str, default=None,
    help='Name of calibrated waveforms in the frame. If None, calibrate from launches.')
parser.add_option('--pulses', dest='pulses', type=str, default=[],
    action='append',
    help='Add a pulse series to plot. If empty, extract from launches.')
parser.add_option('--sim', dest='simhacks', default=False, action='store_true',
    help='Enable simulation hacks in calibration and feature extraction')
parser.add_option('--top', dest='ndoms', type=int, default=None,
    help='Plot the NDOMS brightest DOMs for each event.')
parser.add_option('--key', dest='keys', type="omkey", default=[], action='append',
    help='Plot only waveforms from this DOM. This option may be specified multiple times')
parser.add_option('--decode', dest='decode', default=False,
    help='Decode raw data', action='store_true')
parser.add_option('--show', dest='show', default=True,
	help='Show plots on the screen')
parser.add_option('--save', dest='save', default=False,
	help='Save plots to a png file')

opts, infiles = parser.parse_args()

if len(opts.keys) > 0 and opts.ndoms is not None:
	parser.error("You may specify either --key or --top, but not both!")
if len(infiles) == 0:
	parser.error("You must specify at least one input file!")

import os, sys
from icecube import icetray, dataio, dataclasses, wavereform
from icecube.wavereform.plotter import WaveformPlotter, ShowPlots
import I3Tray

tray = I3Tray.I3Tray()

tray.AddModule("I3Reader", "reader", filenameList=infiles)

class SingleGCDFilter(icetray.I3Module):
	"""
	Sometimes files include GCD frames, and sometimes they're just wrong.
	Prepending a correct GCD file doesn't work, since the GCD frames in the
	second file overwrite the frame cache. Fix it by swallowing all repeated
	GCD frames that don't have at least 1 Q frame in between them.
	"""
	def __init__(self, ctx):
		super(SingleGCDFilter, self).__init__(ctx)
		self.AddOutBox("OutBox")
	def Configure(self):
		self.current = [False, False, False]
	def filter(self, frame, i):
		if not self.current[i]:
			self.PushFrame(frame)
			self.current[i] = True
	def Geometry(self, frame):
		self.filter(frame, 0)
	def Calibration(self, frame):
		self.filter(frame, 1)
	def DetectorStatus(self, frame):
		self.filter(frame, 2)
	def DAQ(self, frame):
		self.current = [False, False, False]
		self.PushFrame(frame)
tray.AddModule(SingleGCDFilter, 'kill_repeat_GCD')

if opts.decode:
	from icecube.payload_parsing import I3DOMLaunchExtractor
	tray.AddModule('QConverter', 'qify', WritePFrame=False)
	tray.AddModule(lambda fr: 'I3DAQData' in fr, 'rawdata_cut', Streams=[icetray.I3Frame.DAQ])
	tray.AddSegment(I3DOMLaunchExtractor, 'extract')

if len(opts.keys) > 0:
	def has_hlc(frame):
		dlsm = frame['InIceRawData']
		if not any([key in dlsm for key in opts.keys]):
			return False
		return any([any([d.lc_bit for d in dlsm[key]]) for key in opts.keys])
	tray.AddModule(has_hlc, 'key_select', Streams=[icetray.I3Frame.DAQ])

if opts.waveforms is None:
	
	def filter_borked_slc(frame):
		dlsm = frame['InIceRawData']
		ndlsm = dlsm.__class__()
		for om, dls in dlsm:
			ndls = []
			for dl in dls:
				if not ((not dl.lc_bit) and dl.charge_stamp_highest_sample == 0):
					ndls.append(dl)
			ndlsm[om] = dls.__class__(ndls)
		del frame['InIceRawData']
		frame['InIceRawData'] = ndlsm
	tray.AddModule(filter_borked_slc, 'daqbug?', Streams=[icetray.I3Frame.DAQ])
	
	opts.waveforms = "CalibratedWaveforms"
	
	tray.Add("Delete", keys=["CalibratedWaveformRange"])
	
	if opts.simhacks:
		from icecube.WaveCalibrator import DOMSimulatorCalibrator
		tray.AddSegment(DOMSimulatorCalibrator, "calibrate",
		    Waveforms=opts.waveforms,
		    Errata="DroopCorrectionFailure",
		)
		
	else:
		from icecube import WaveCalibrator
		tray.AddModule("I3WaveCalibrator", "calibrate",
		    Waveforms=opts.waveforms,
		    Errata="DroopCorrectionFailure",
		)
	
	tray.AddModule("I3PMTSaturationFlagger", "find_saturation",
	    Waveforms=opts.waveforms,
	    Output="PMTSaturation",
	)

if len(opts.pulses) == 0:
	icetray.load('wavedeform', False)
	opts.pulses.append("WavedeformPulses")
	if opts.simhacks:
		tray.AddModule("I3Wavedeform", "deform",
		    Waveforms=opts.waveforms,
		    Output=opts.pulses[-1],
		    UseDOMsimulatorTemplates=True,
		)
	else:
		tray.AddModule("I3Wavedeform", "deform",
		    Waveforms=opts.waveforms,
		    Output=opts.pulses[-1],
		)
	
tray.AddModule('Dump', 'dumpy')

def WaveformSelector(frame, pulses, waveforms, number=10):
	pulsemap = None
	if type(frame[pulses]) == dataclasses.I3RecoPulseSeriesMap:
		pulsemap = frame.Get(pulses)
	elif (type(frame[pulses]) == dataclasses.I3RecoPulseSeriesMapMask or 
		type(frame[pulses]) == dataclasses.I3RecoPulseSeriesMapUnion):
		pulsemap = frame[pulses].apply(frame)
	else:
		icetray.logging.log_fatal("Unable to handle pulses {} of type {}".format(pulses, type(pulses)))
	chargemap = dict([(k, sum([p.charge for p in v])) for k, v in pulsemap])
	keys = list(chargemap.keys())
	keys.sort(key=lambda k: chargemap[k])
	
	wfsm = frame[waveforms]
	subwfsm = wfsm.__class__()

	qtot = sum(chargemap.values())
	nchan = len(chargemap)

	if len(opts.keys) > 0:
		for key in opts.keys:
			if key in wfsm:
				subwfsm[key] = wfsm[key]
	else:
		for om in keys[::-1][:number]:
			subwfsm[om] = wfsm[om]
		
	if len(subwfsm) == 0:
		return False
	
	frame['Selected%s' % waveforms] = subwfsm
	return True

def PulseSummary(frame, pulsenames):
	print('-'*79)
	for name in pulsenames:
		pulsemap = None
		if type(frame[name]) == dataclasses.I3RecoPulseSeriesMap:
			pulsemap = frame.Get(name)
		elif (type(frame[name]) == dataclasses.I3RecoPulseSeriesMapMask or 
			type(frame[name]) == dataclasses.I3RecoPulseSeriesMapUnion):
			pulsemap = frame[name].apply(frame)
		else:
			icetray.logging.log_fatal("Unable to handle pulses {} of type {}".format(name, type(name)))
		hlccharge = dict()
		slccharge = dict()
		nhlc = 0
		nslc = 0
		flags = dataclasses.I3RecoPulse.PulseFlags
		for om, series in pulsemap:
			hlc = [p.charge for p in series if p.flags & flags.LC]
			slc = [p.charge for p in series if not p.flags & flags.LC]
			hlcq = sum(hlc)
			slcq = sum(slc)
			if hlcq > 0:
				hlccharge[om] = hlcq
			if slcq > 0:
				slccharge[om] = slcq
			nhlc += len(hlc)
			nslc += len(slc)

		print('%s: HLC: %d chan/%d pulses/%e PE SLC: %d chan/%d pulses/%e PE' % (name, len(hlccharge), nhlc, sum(hlccharge.values()), len(slccharge), nslc, sum(slccharge.values())))
	print('-'*79)
	if len(opts.keys):
		print([p.time for p in pulsemap[opts.keys[0]]])

tray.AddModule(WaveformSelector, "selecty_%s" % opts.pulses[0],
	waveforms=opts.waveforms,
	pulses=opts.pulses[0],
	number=opts.ndoms,
	Streams=[icetray.I3Frame.DAQ],
)

tray.AddModule(PulseSummary, 'summary',
	pulsenames=opts.pulses,
	Streams=[icetray.I3Frame.DAQ],
	)

for pulsename in opts.pulses:

	tray.AddModule(WaveformPlotter, "plotesque_%s" % pulsename,
		Waveforms="Selected%s" % opts.waveforms,
		Errata=["DroopCorrectionFailure", "PMTSaturation"],
		Pulses=pulsename,
		# ShowChiSquared=False,
		Block=False,
	)

tray.AddModule(ShowPlots, "show",
	DisplayPlots=opts.show,
	SavePlots=opts.save,
	OutFilename='{}_{}_{}_{}.png'.format(opts.keys[0].string, opts.keys[0].om, 
		opts.waveforms, opts.pulses[0])
	)



tray.Execute()


