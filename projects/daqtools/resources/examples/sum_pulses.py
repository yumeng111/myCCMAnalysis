#!/usr/bin/env python3
import numpy as np
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.dataclasses import CCMOMGeo
from icecube.icetray import CCMTriggerKey, CCMPMTKey
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)

class FilterSumPulses(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.n_daq_frames = 0
        self.n_frames = 0
        self.n_pulses = {}

    def Configure(self):
        self.output_files = {}
        pass

    def Geometry(self, frame):
        geo = frame['CCMGeometry']
        self.pmt_channel_map = geo.pmt_channel_map
        self.trigger_copy_map = geo.trigger_copy_map

        #self.output_files = {}
        #for pmt, channel in self.pmt_channel_map.items():
        #    file = open(f"pulse_samples_{pmt[0]}_{pmt[1]}.txt", "a")
        #    self.output_files[pmt] = file
        #    self.n_pulses[pmt] = 0

        self.PushFrame(frame)

    def DAQ(self, frame):
        #if len(frame["NIMPulses"][CCMTriggerKey(CCMTriggerKey.TriggerType.StrobeTrigger, 1)]) == 0:
        #    self.n_frames += 1
        #    self.PushFrame(frame)
        #    return
        if self.n_daq_frames % 10 == 0:
            print("N Frames:", np.median(list(self.n_pulses.values())), "/", self.n_daq_frames, end="\r")
        self.n_daq_frames += 1
        self.n_frames += 1
        if not frame.Has("SummedPulses"):
            return
        if not frame.Has("SummedPulsesPeakPositions"):
            return
        if not frame.Has("SummedPulsesCounts"):
            return
        print("\n")
        self.PushFrame(frame)

    def Finish(self):
        #for pmt, file in self.output_files.items():
        #    file.close()
        pass

if __name__ == "__main__":

    import argparse
    # Input arguments
    parser = argparse.ArgumentParser(
            prog = "Sum waveforms",
            description = "Sum the waveforms from certain PMTs together, accounting for timing differences, and normalizing to a single channel",
            )
    parser.add_argument("--files",
            type=str,
            nargs="+",
            action="append",
            default=[],
            help="Files to Process")
    parser.add_argument("--num-events",
            type=int,
            default=0,
            help="Number of events to process")
    args = parser.parse_args()

    if len(args.files) == 0:
        args.files = [["/lustre/scratch4/turquoise/aschneider/data/2022/merged/run012569/*.i3.zst"]]

    import glob
    print("Glob:", args.files)
    fnames = [s for l in args.files for g in l for s in sorted(glob.glob(g))]
    print("Found these files:", fnames)

    tray = I3Tray()
    tray.Add("I3Reader", "reader", FilenameList=fnames)
    tray.Add("PulseCollector")
    tray.Add("SumPulses")
    tray.Add(FilterSumPulses)
    tray.Add("Keep", Keys=["SummedPulses", "SummedPulsesPeakPositions", "SummedPulsesCounts"])
    tray.Add("I3Writer", "writer", FileName="test.i3.zst", Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

    #tray.Add("Dump") # Prints out the names of the objects in every frame
    if args.num_events < 1:
        tray.Execute() # Process all frames
    else:
        tray.Execute(args.num_events + 1) # Number of frames to process is num_events DAQ frames plus one Geometry frame

