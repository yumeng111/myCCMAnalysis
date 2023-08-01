### This is an example script to run the module CosmicTriggerSumWaveforms.cxx
### on the unmerged data streams


#!/usr/bin/env python3
import numpy as np
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses, phys_services
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.dataclasses import CCMOMGeo
from icecube.icetray import CCMTriggerKey, CCMPMTKey
import time
import uuid
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)

def get_filename(output_dir, prefix=None):
    timestr = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    uid = str(uuid.uuid4())
    fname = timestr + '_' + uid
    if prefix is not None:
        fname = args.output_dir + prefix + '_' + fname
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
            default="",
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
    parser.add_argument("--merged",
            default=False,
            action="store_true",
            )
    args = parser.parse_args()

    if len(args.daq_directories) == 0:
        args.daq_directories = ["mills", "wills"]
    else:
        args.daq_directories = [y for x in args.daq_directories for y in x]

    if args.input_directory is None:
        if args.merged:
#            args.input_directory = "/lustre/scratch4/turquoise/aschneider/data/2022/merged"
            args.input_directory = "/lustre/scratch4/turquoise/marisolc/data/2022/merged_data"
        else:
#            args.input_directory = "/lustre/scratch4/turquoise/aschneider/data/2022/separate_daqs"
            args.input_directory = "/lustre/scratch4/turquoise/marisolc/data/2022/separate_daqs"

    if args.run is None:
        #args.run = 12569
        args.run = 12486

    if args.output_file is None and args.run is not None:
        run_prefix = "run%06d" % args.run
        args.output_file = get_filename(args.output_dir, run_prefix) + ".i3.zst"
    else:
        args.output_file = get_filename(args.output_dir) + ".i3.zst"

    if len(args.files) == 0:
        run_prefix = "run%06d" % args.run
        if args.merged:
            args.files = [[args.input_directory + "/" + run_prefix + "/*.i3.zst"]]
        else:
            args.files = [[args.input_directory + "/" + s + "/" + run_prefix + "/*.i3.zst"] for s in args.daq_directories]

    import glob
    print("Glob:", args.files)
    if args.merged:
        fnames = [s for l in args.files for g in l for s in sorted(glob.glob(g))]
    else:
        fnames = [[s for g in l for s in sorted(glob.glob(g))] for l in args.files]
    print("Found these files:", fnames)

    tray = I3Tray()
    tray.context["I3RandomService"] = phys_services.I3GSLRandomService(42)
    if args.merged:
        tray.Add("I3Reader", "reader", FilenameList=fnames)
        tray.Add(daqtools.GeometryReplacer.GeometryReplacer, GeometryFile="/usr/projects/w20_ccm_lanl/geometry/2022/replacements/2022_geometry_replacement_run012298-run012756_2023-07-22.i3.zst")
    else:
        tray.Add("MergedSource", "reader", FileLists=fnames, MaxTimeDiff=32)
        tray.Add("Delete", Keys=["CCMGeometry"]) # Remove keys that we are computing in this example
        # Generate the CCMGeometry object based on the CCMDAQConfig and the (PMT name --> position) mapping
        # SquashDuplicateConfigs=True prevents us from writing duplicate Geometry frames
        tray.Add(daqtools.GeometryReplacer.GeometryReplacer, GeometryFile="/usr/projects/w20_ccm_lanl/geometry/2022/replacements/2022_geometry_replacement_run012298-run012756_2023-07-22.i3.zst")
        tray.Add("CCMTriggerMerger", "trigger_merger")
        tray.Add("NIMLogicPulseFinder", "nim_pulses")
        tray.Add("BeamCurrentMonitorSummary", "bcm_summary")
    N_samples_before_baseline = 400
    N_samples_for_baseline = 50
    max_pulse_start_sample = N_samples_before_baseline + N_samples_before_baseline + 50
    tray.Add("FirstPulseFinder", "pulse_finder", NumSamplesBeforePulse=N_samples_before_baseline, MaxPulseStartSample=max_pulse_start_sample)
    tray.Add("BaselineEstimator", "baseline_estimator", NumSamples=N_samples_for_baseline, NumFramesForEstimate=25)
    tray.Add("Dump")
    tray.Add("Delete", Keys=["PulsePositions", "WaveformSmoothers"])
    tray.Add("SumWaveforms", AllowedTriggerKeys=[CCMTriggerKey(CCMTriggerKey.StrobeTrigger, 1)])
    tray.Add("I3Writer", "writer", FileName=args.output_file, Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
    #tray.Add("Dump") # Prints out the names of the objects in every frame
    if args.num_events < 1:
        tray.Execute() # Process all frames
    else:
        tray.Execute(args.num_events + 1) # Number of frames to process is num_events DAQ frames plus one Geometry frame

