#!/usr/bin/env python3
import numpy as np
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.dataclasses import CCMOMGeo
from icecube.icetray import CCMTriggerKey, CCMPMTKey
import time
import uuid
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)

def get_filename(prefix=None):
    timestr = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    uid = str(uuid.uuid4())
    fname = timestr + '_' + uid
    if prefix is not None:
        fname = prefix + '_' + fname
    return fname


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
    parser.add_argument("--output-dir",
            type=str,
            help="Output directory")
    parser.add_argument("--daq-directories",
            type=str,
            nargs="+",
            action="append",
            default=[],
            help="Sub directories for runs from each DAQ machine")
    parser.add_argument("--input-directory",
            type=str,
            help="The top-level input directory")
    parser.add_argument("--run",
            type=int,
            help="The run number")
    parser.add_argument("--output-file",
            type=str,
            help="Output file name and path")
    parser.add_argument("--num-events",
            type=int,
            default=0,
            help="Number of events to process")
    args = parser.parse_args()

    if len(args.daq_directories) == 0:
        args.daq_directories = ["mills", "wills"]
    else:
        args.daq_directories = [y for x in args.daq_directories for y in x]

    if args.input_directory is None:
        args.input_directory = "/lustre/scratch4/turquoise/aschneider/data/2022/separate_daqs"

    if args.run is None:
        args.run = 12569

    if args.output_file is None and args.run is not None:
        run_prefix = "run%06d" % args.run
        args.output_file = get_filename(run_prefix) + ".i3.zst"
    else:
        args.output_file = get_filename() + ".i3.zst"

    if len(args.files) == 0:
        run_prefix = "run%06d" % args.run
        args.files = [[args.input_directory + "/" + s + "/" + run_prefix + "/*.i3.zst"] for s in args.daq_directories]

    import glob
    print("Glob:", args.files)
    fnames = [[s for g in l for s in sorted(glob.glob(g))] for l in args.files]
    print("Found these files:", fnames)

    tray = I3Tray()
    tray.Add("MergedSource", "reader", FileLists=fnames, MaxTimeDiff=32)
    tray.Add("Delete", Keys=["CCMGeometry"]) # Remove keys that we are computing in this example
    # Generate the CCMGeometry object based on the CCMDAQConfig and the (PMT name --> position) mapping
    # SquashDuplicateConfigs=True prevents us from writing duplicate Geometry frames
    tray.Add("CCMGeometryGenerator", "geometry_generator", SquashDuplicateConfigs=True)
    tray.Add("CCMTriggerMerger", "trigger_merger")
    tray.Add("PulseCollector")
    tray.Add("SumPulses")
    tray.Add(FilterSumPulses)
    tray.Add("Keep", Keys=["SummedPulses", "SummedPulsesPeakPositions", "SummedPulsesCounts"])
    tray.Add("I3Writer", "writer", FileName=args.output_file, Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

    #tray.Add("Dump") # Prints out the names of the objects in every frame
    if args.num_events < 1:
        tray.Execute() # Process all frames
    else:
        tray.Execute(args.num_events + 1) # Number of frames to process is num_events DAQ frames plus one Geometry frame

