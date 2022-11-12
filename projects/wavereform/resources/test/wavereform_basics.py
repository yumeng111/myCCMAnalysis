#!/usr/bin/env python3

from I3Tray import *
import os, random, math, sys
from icecube import icetray, dataio, dataclasses
from icecube import phys_services, wavereform
from os.path import expandvars
from copy import deepcopy

i3_testdata = expandvars("$I3_TESTDATA") 

tray = I3Tray()
tray.AddModule('I3InfiniteSource', 'reader', Prefix = i3_testdata + '/GCD/GeoCalibDetectorStatus_2012.56063_V0.i3.gz')

if sys.version_info[0] >= 3:
    random.seed(0)
else:
    random.seed(5)

def RandomWaveforms(fr):
    calib = fr['I3Calibration']
    status = fr['I3DetectorStatus']

    pulsemap = dataclasses.I3RecoPulseSeriesMap()
    borkedPulseMap = dataclasses.I3RecoPulseSeriesMap()
    wfmap = dataclasses.I3WaveformSeriesMap()

    for om in calib.dom_cal.keys():
        pulses = dataclasses.I3RecoPulseSeries()
        waveforms = dataclasses.I3WaveformSeries()
        borked = ( random.randint(0, 10) == 0 )

        fadc_templ = calib.dom_cal[om].fadc_pulse_template(False)
        atwd0_templ = calib.dom_cal[om].atwd_pulse_template(0, False)

        # Skip DOMs without status info
        if om not in status.dom_status:
            continue

        # Skip DOMs that are off
        if (status.dom_status[om]).pmt_hv == 0:
            continue 
        
        # Filter out DOMs with ATWD problems (e.g. Coxae)
        if ((calib.dom_cal[om]).atwd_bin_calib_slope[(0, 0, 0)] == 0.):
            continue

        spe_charge = dataclasses.spe_mean(status.dom_status[om],
            calib.dom_cal[om])
        if spe_charge < I3Units.pC:
            continue
        spe_charge *= calib.dom_cal[om].front_end_impedance

        for launch in range(0, random.randint(1, 4)):
            npulses = random.randint(0, 5)
            launchtime = launch*10000
            # Make 30% of SPE launches SLC
            slc = (npulses == 1 and random.uniform(0,1) < 0.3)

            # ATWD Waveform
            atwd0_wf = dataclasses.I3Waveform()
            atwd0_wf.waveform = [ \
                random.normalvariate(0, 0.3)*I3Units.mV for \
                i in range(0, 128)]
            atwd0_wf.digitizer = dataclasses.I3Waveform.ATWD
            atwd0_wf.bin_width = 3.3
            atwd0_wf.hlc = not slc
            atwd0_wf.time = launchtime

            # FADC Waveform
            if slc:
                fadc_nbins = 3
            else:
                fadc_nbins = 256
            fadc_wf = dataclasses.I3Waveform()
            fadc_wf.waveform = [ \
                random.normalvariate(0, 0.1)*I3Units.mV for \
                i in range(0, fadc_nbins)]
            fadc_wf.digitizer = dataclasses.I3Waveform.FADC
            fadc_wf.bin_width = 25
            fadc_wf.hlc = not slc
            fadc_wf.time = launchtime

            for p in range(0, npulses):
                pulse = dataclasses.I3RecoPulse()
                pulse.charge = random.randint(1, 3)
                if not slc:
                    pulse.time = launchtime + \
                        random.gammavariate(2.5, 80)
                    pulse.flags = pulse.PulseFlags.LC
                else:
                    pulse.time = launchtime + \
                        random.uniform(-25, 25)
                # Simulate wavedeform failure
                finalPulse = deepcopy(pulse)
                if borked:
                    finalPulse.time += 25
                    finalPulse.charge += 100
                pulses.append(finalPulse)
                

                norm = spe_charge * pulse.charge
                for i in range(0, len(fadc_wf.waveform)):
                    fadc_wf.waveform[i] += norm * \
                        fadc_templ((i+1)*fadc_wf.bin_width - \
                        (pulse.time - launchtime))
                for i in range(0, len(atwd0_wf.waveform)):
                    atwd0_wf.waveform[i] += norm * \
                        atwd0_templ((i+1)*atwd0_wf.bin_width - \
                        (pulse.time - launchtime))

            waveforms.append(fadc_wf)
            if not slc:
                waveforms.append(atwd0_wf)
        wfmap[om] = waveforms
        pulsemap[om] = dataclasses.I3RecoPulseSeries(sorted(pulses,
            key=lambda pulse: pulse.time))
        if borked:
            borkedPulseMap[om] =  dataclasses.I3RecoPulseSeries(sorted(pulses,
                key=lambda pulse: pulse.time))
    fr['RandomPulses'] = pulsemap
    fr['BorkedPulses'] = borkedPulseMap
    fr['CalibratedWaveforms'] = wfmap
    fr['CalibratedWaveformRange'] = dataclasses.I3TimeWindow(0,10000)
                
    
tray.AddModule(RandomWaveforms, 'random', Streams=[icetray.I3Frame.DAQ])

# Perform some sanity checks.
tray.AddModule("I3Wavereform", "wavereform",
    Waveforms="CalibratedWaveforms",
    Pulses="RandomPulses",
    Chi="Chi",
    ChiThreshold=1e4, # flag ~ 0.2% of waveforms
    Flag="Borked",
)

import unittest
class SanityCheck(unittest.TestCase):
    fakePulseKey = "RandomPulses"
    borkedPulseKey = "BorkedPulses"
    flaggedKey = "Borked"
    chiKey = "Chi"

    def testKeys(self):
        self.assert_(self.chiKey in self.frame, "Chi is in the frame.")
        self.assert_(self.flaggedKey in self.frame, "Flagged OM list is in the frame.")
    def testChi(self):
        psm = self.frame[self.fakePulseKey]
        chi = self.frame[self.chiKey]
        borkedPulses = self.frame[self.borkedPulseKey]
        flaggedOMs = self.frame[self.flaggedKey]
        for om, vec in psm:
            if om in borkedPulses and len(psm[om]) > 0:
                self.assert_(any(x > 1e4 for x in chi[om]), "Chi for borked pulses is above threshold")
                self.assert_(om in flaggedOMs, "OMs with borked pulses are flagged")
            else:
                self.assert_(all(x < 1e4 for x in chi[om]), "Chi for good pulses is below threshold")

tray.AddModule(icetray.I3TestModuleFactory(SanityCheck), 'testy', Streams=[icetray.I3Frame.DAQ])
tray.AddModule('I3Writer', 'writer', Filename='wavereform_test.i3')

tray.Execute(6)


os.unlink('wavereform_test.i3')

