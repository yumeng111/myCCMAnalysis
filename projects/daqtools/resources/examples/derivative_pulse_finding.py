#!/usr/bin/env python3
import numpy as np
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses, phys_services, daqtools
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.dataclasses import CCMOMGeo
from icecube.icetray import CCMTriggerKey, CCMPMTKey
from icecube import wavedeform, wavereform
import time
import uuid
import os

icetray.set_log_level('INFO')
icetray.logging.set_level("INFO")

load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)


class FindDesiredFrames(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter(
            "DesiredTriggerType",
            "Trigger type you want to select for and continue processing data",
            "",
        )
        self.AddParameter(
            "FrameTypesToFilter",
            "The frame types to filter on",
            [I3Frame.DAQ, I3Frame.Physics],
        )
        self.n_frames = 0

    def Configure(self):
        self.desired_trigger_type = self.GetParameter("DesiredTriggerType")
        self.frame_types = self.GetParameter("FrameTypesToFilter")

    def IsCorrectTriggerType(self, frame):
        ### let's get the nim pulse
        nim_pulses = frame["NIMPulses"]
        trigger_key = CCMTriggerKey(self.desired_trigger_type, 1)
        if trigger_key in nim_pulses and len(nim_pulses[trigger_key]) > 0:
            return True
        return False

    def Process(self):
        frame = self.PopFrame()
        if frame.Stop in self.frame_types:
            if self.IsCorrectTriggerType(frame):
                self.PushFrame(frame)
            else:
                return
        else:
            self.PushFrame(frame)

    def Finish(self):
        pass


class InjectCalibrationFrame(I3Module):
    def __init__(self, context):
        I3Module.__init__(self, context)
        self.AddParameter(
            "CalibrationFile",
            "The path to the i3 file containing the calibration frame to inject",
            "",
        )
        self.AddParameter(
            "FrameTypesAfterInjection",
            "The frame types that come after the Calibration frame",
            [I3Frame.DetectorStatus, I3Frame.DAQ, I3Frame.Physics],
        )
        self.AddParameter(
            "FrameTypesBeforeInjection",
            "The frame types that come before the Calibration frame",
            [I3Frame.Geometry],
        )

    def Configure(self):
        self.calibration_filename = self.GetParameter("CalibrationFile")
        self.frame_types_before = self.GetParameter("FrameTypesBeforeInjection")
        self.frame_types_after = self.GetParameter("FrameTypesAfterInjection")
        self.cached_frames = []

        i3file = dataio.I3File(self.calibration_filename)
        self.calibration_frame = None
        while i3file.more():
            calibration_frame = i3file.pop_frame()
            if calibration_frame.Stop == icetray.I3Frame.Calibration:
                self.calibration_frame = calibration_frame
                break
        if self.calibration_frame is None:
            raise RuntimeError("Could not find Calibration frame in supplied file")

    def Process(self):
        frame = self.PopFrame()
        if frame.Stop in self.frame_types_before:
            self.cached_frames.append(frame)
            return
        elif frame.Stop in self.frame_types_after and len(self.cached_frames) > 0:
            for cached_frame in self.cached_frames:
                self.PushFrame(cached_frame)
            self.PushFrame(self.calibration_frame)
            self.PushFrame(frame)
            self.cached_frames = []
        else:
            self.PushFrame(frame)

    def Finish(self):
        if len(self.cached_frames) > 0:
            for cached_frame in self.cached_frames:
                self.PushFrame(cached_frame)
            self.PushFrame(self.calibration_frame)
            self.PushFrame(frame)
            self.cached_frames = []


if __name__ == "__main__":
    import argparse

    # Input arguments
    parser = argparse.ArgumentParser(
        description="Read files from separate DAQs, run basic processing, and find pulses with the derivative method",
    )
    parser.add_argument(
        "--files",
        type=str,
        nargs="+",
        action="append",
        default=[],
        help="Files to Process",
    )
    parser.add_argument("--output-dir", type=str, default="", help="Output directory")
    parser.add_argument(
        "--daq-directories",
        type=str,
        nargs="+",
        action="append",
        default=[],
        help="Sub directories for runs from each DAQ machine",
    )
    parser.add_argument(
        "--input-directory",
        type=str,
        default=None,
        help="The top-level input directory",
    )
    parser.add_argument(
        "--in-run-folders",
        type=bool,
        default=False,
    )
    parser.add_argument("--run", type=int, help="The run number")
    parser.add_argument(
        "--output-prefix", type=str, default="%04d", help="Output file name"
    )
    parser.add_argument(
        "--num-events", type=int, default=0, help="Number of events to process"
    )
    parser.add_argument("--calibration-file", type=str, required=True)
    parser.add_argument("--geometry-file", type=str, required=True)
    args = parser.parse_args()

    if len(args.daq_directories) == 0:
        args.daq_directories = ["millstester", "willstester"]
    else:
        args.daq_directories = [y for x in args.daq_directories for y in x]

    if args.input_directory is None:
        args.input_directory = (
            "/lustre/scratch4/turquoise/aschneider/data/2022/separate_daqs"
        )

    output_prefix = args.output_prefix
    if len(args.files) > 0:
        output_dir = args.output_dir
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{run_prefix}_{output_prefix}.i3.zst")
        geo_output_file = args.output_dir + f"GCD_{run_prefix}_{output_prefix}.i3.zst"
    elif args.run is not None:
        run_prefix = "run%06d" % args.run
        if args.in_run_folders:
            args.files = [
                [args.input_directory + "/" + s + "/" + run_prefix + "/*.i3.zst"]
                for s in args.daq_directories
            ]
        else:
            args.files = [
                [args.input_directory + "/" + s + "/*.i3.zst"]
                for s in args.daq_directories
            ]
        output_dir = os.path.join(args.output_dir, run_prefix)
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{output_prefix}.i3.zst")
        geo_output_file = args.output_dir + f"GCD_{output_prefix}.i3.zst"

    import glob

    print("Glob:", args.files)
    fnames = [[s for g in l for s in sorted(glob.glob(g))] for l in args.files]
    print("Found these files:")
    for s in fnames:
        print(s)

    tray = I3Tray()
    # Read in files from two data streams and align them
    # This also produces a geometry object from the CCMDAQConfig contained in the files
    tray.Add("MergedSource", "reader", FileLists=fnames, MaxTimeDiff=32)
    # Replace the geometry and CCMDAQConfig with the patched version
    #tray.Add(
    #    daqtools.GeometryReplacer.GeometryReplacer,
    #    GeometryFile=args.geometry_file,
    #)
    # Merge triggers that are overlapping
    tray.Add("CCMGeometryGenerator")
    tray.Add("CCMTriggerMerger", "trigger_merger")
    # Add the calibration frame to the data stream since the MergedSource reader prevents us from reading in other files
    tray.Add(InjectCalibrationFrame, CalibrationFile=args.calibration_file)
    # Find the NIM pulses
    tray.Add("NIMLogicPulseFinder", "nim_pulses")
    # Examine the beam pulse
    tray.Add("BeamCurrentMonitorSummary", "bcm_summary")
    # Filter out DAQ frames for certain triggers
    tray.Add(
        FindDesiredFrames, DesiredTriggerType=CCMTriggerKey.TriggerType.CosmicTrigger
    )
    tray.Add("BaselineEstimator")
    tray.Add("Rename", Keys=["CCMPMTCalibration", "CCMCalibration"])
    tray.Add("ElectronicsCorrection")
    tray.Add("CCMDerivativePulseFinder")
    # Write the GCD frames to a separate file
    tray.AddModule(
        "I3Writer",
        "geo writer",
        FileName=geo_output_file,
        Streams=[
            icetray.I3Frame.Geometry,
            icetray.I3Frame.Calibration,
            icetray.I3Frame.DetectorStatus,
        ],
    )
    # Write all the DAQ frames to their own set of files
    if "%" in output_file:
        # Split the output into files of a desired size
        tray.AddModule(
            "I3MultiWriter",
            "multi writer",
            SizeLimit=(3 * 1024**3),
            FileName=output_file,
            MetadataStreams=[
                icetray.I3Frame.Geometry,
                icetray.I3Frame.Calibration,
                icetray.I3Frame.DetectorStatus,
            ],
            CounterStream=icetray.I3Frame.DAQ,
            SizeCheckInterval=1000,
            NWorkers=8,
            Streams=[
                icetray.I3Frame.DetectorStatus,
                icetray.I3Frame.DAQ,
                icetray.I3Frame.Physics,
            ],
        )
    else:
        # Write all the output to one file
        tray.AddModule(
            "I3Writer",
            "writer",
            FileName=output_file,
            Streams=[
                icetray.I3Frame.Calibration,
                icetray.I3Frame.DetectorStatus,
                icetray.I3Frame.DAQ,
                icetray.I3Frame.Physics,
            ],
        )
    tray.Add("Dump")  # Print out the names of the objects in every frame
    if args.num_events < 1:
        tray.Execute()  # Process all frames
    else:
        tray.Execute(
            args.num_events + 1
        )  # Number of frames to process is num_events DAQ frames plus one Geometry frame
