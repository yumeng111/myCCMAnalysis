
from I3Tray import *
import random, sys
from icecube import icetray, dataclasses

random.seed(0)


def RandomWaveforms(fr):
    calib = fr['I3Calibration']
    status = fr['I3DetectorStatus']

    pulsemap = dataclasses.I3RecoPulseSeriesMap()
    wfmap = dataclasses.I3WaveformSeriesMap()

    for om in calib.dom_cal.keys():
        # Only use string-21 DOMs for historical reasons.  The current tests
        # are designed too sensitive to random fluctuations that would be
        # considered as normal data.  See ticket #1679.
        if om.string != 21:
            continue
        pulses = dataclasses.I3RecoPulseSeries()
        waveforms = dataclasses.I3WaveformSeries()

        fadc_templ = calib.dom_cal[om].fadc_pulse_template(False)
        atwd0_templ = calib.dom_cal[om].atwd_pulse_template(0, False)

        spe_charge = dataclasses.spe_mean(status.dom_status[om],
            calib.dom_cal[om])
        if spe_charge < I3Units.pC: continue
        spe_charge *= calib.dom_cal[om].front_end_impedance

        for launch in range(0, random.randint(0, 4)):
            npulses = random.randint(0, 5)
            launchtime = launch*10000
            # Make 30% of SPE launches SLC
            slc = (npulses == 1 and random.uniform(0,1) < 0.3)

            # ATWD Waveform
            atwd0_wf = dataclasses.I3Waveform()
            atwd0_wf.waveform = [ \
                random.normalvariate(0, 0.2)*I3Units.mV for \
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
                    if (pulse.time > (launchtime + 125 * atwd0_wf.bin_width)):
                        # We missed it
                        continue
                    pulse.flags = pulse.PulseFlags.LC
                else:
                    pulse.time = launchtime + \
                        random.uniform(-25, 25)
                pulses.append(pulse)

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
    fr['RandomPulses'] = pulsemap
    fr['CalibratedWaveforms'] = wfmap
    fr['CalibratedWaveformRange'] = dataclasses.I3TimeWindow(0,10000)
