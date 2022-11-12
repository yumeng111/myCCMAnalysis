#!/usr/bin/env python3

"""
An example of decoding, calibrating, feature-extracting, and sanity-checking raw DAQ data.
Raw DOMLaunches are included in the output whenever the I3Wavereform chi^2 sanity check fails.
"""

import os, sys
from icecube import icetray, dataio, WaveCalibrator, wavedeform, wavereform
from icecube.payload_parsing import I3DOMLaunchExtractor
import I3Tray

infiles = sys.argv[1:-1]
outfile = sys.argv[-1]

tray = I3Tray.I3Tray()

tray.AddModule("I3Reader", "reader", filenameList=infiles)

# Decode raw data
tray.AddModule('QConverter', 'qify', WritePFrame=False)
tray.AddModule(lambda fr: 'I3DAQData' in fr, 'rawdata_cut', Streams=[icetray.I3Frame.DAQ])
tray.AddSegment(I3DOMLaunchExtractor, 'extract')
tray.AddModule(lambda frame: "InIceRawData" in frame, "picky")

# Calibrate (all together now)
tray.AddModule('I3WaveCalibrator', 'domcal')

# Extract pulses.
tray.AddModule('I3Wavedeform', 'deform')

# Perform some sanity checks.
tray.AddModule("I3Wavereform", "wavereform_atwd",
	Waveforms="CalibratedWaveforms",
	Pulses="WavedeformPulses",
	Chi="Chi_WavedeformPulses_CalibratedWaveforms",
	ChiThreshold=1e4, # flag ~ 0.2% of waveforms, decrease if you want more flagged waveforms
	Flag="Borked_Pulses",
	)

# Save raw data for those OMs where calibration or feature extraction
# failed in an obvious way. Two keys are added to the frame:
# InIceErrata => raw data for screwy DOMs (keep only this if discarding the DAQ payload)
# InIceErrataKeys => compact list of screwy keys (keep only this if transmitting the DAQ payload) 
tray.AddModule("I3LaunchSelector", "seatbelt",
	Launches="InIceRawData",
	Flags=["CalibrationErrata", "Borked_Pulses"],
	Output="InIceErrata",
	)

# Optional: Write the chi^2 to an HDF5 table.
# from icecube import tableio, hdfwriter, dataclasses
# tabler=hdfwriter.I3HDFTableService(outfile + '.hdf5')
# tray.AddModule(tableio.I3TableWriter, "scribe",
# 	tableservice=tabler,
# 	types=[dataclasses.I3MapKeyVectorDouble],
# 	)

# Clean out the frame
keep_keys = ['I3EventHeader', 'InIceErrata']
tray.AddModule('Keep', 'trapper-keeper', keys=keep_keys)

tray.AddModule("I3Writer", "writer", streams=[icetray.I3Frame.DAQ], filename=outfile)


tray.Execute()


