#!/usr/bin/env python3
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)

if __name__ == "__main__":

    import argparse
    # Input arguments
    parser = argparse.ArgumentParser(
            prog = "Merge DAQ streams",
            description = "Merges CCM i3 files from multiple DAQs",
            )
    parser.add_argument("--files",
            type=str,
            nargs="+",
            action="append",
            default=[],
            help="Files to Process")
    parser.add_argument("--output-file",
            type=str,
            required=True,
            help="Output file name")
    parser.add_argument("--geometry-output-file",
            type=str,
            required=True,
            help="Output file name")
    parser.add_argument("--num-events",
            type=int,
            default=0,
            help="Number of events to process")
    args = parser.parse_args()

    if len(args.files) == 0:
        args.files = [
            ["/lustre/scratch4/turquoise/aschneider/data/2022/separate_daqs/mills/run012569/*.i3.zst"],
            ["/lustre/scratch4/turquoise/aschneider/data/2022/separate_daqs/wills/run012569/*.i3.zst"],
        ]

    import glob
    print("Glob:", args.files)
    fnames = [[s for g in l for s in sorted(glob.glob(g))] for l in args.files]
    print("Found these files:")
    for s in fnames:
        print(s)

    tray = I3Tray()
    # Read events from multiple DAQ streams and merge them into a steam of frames
    tray.Add("MergedSource", "reader", FileLists=fnames, MaxTimeDiff=32)
    tray.Add("Delete", Keys=["CCMGeometry"]) # Remove keys that we are computing in this example
    # Generate the CCMGeometry object based on the CCMDAQConfig and the (PMT name --> position) mapping
    # SquashDuplicateConfigs=True prevents us from writing duplicate Geometry frames
    tray.Add("CCMGeometryGenerator", "geometry_generator", SquashDuplicateConfigs=True)
    # Write the geometry to a separate file
    tray.Add("I3Writer", "geometry_writer", FileName=args.geometry_output_file, Streams=[icetray.I3Frame.Geometry])
    tray.Add("CCMTriggerMerger", "trigger_merger")

    tray.Add("Delete", Keys=["NIMPulses", "BCMSummary"]) # Remove keys that we are computing in this example

    # Add the module that finds NIM pulses in the trigger channels
    # All parameters are set to there defaults here, so the lines below are quivalent to:
    #   tray.Add("NIMLogicPulseFinder", "nim_pulses")
    tray.Add("NIMLogicPulseFinder", "nim_pulses",
        CCMGeometryName="CCMGeometry", # Frame key for CCMGeometry
        CCMWaveformsName="CCMWaveforms", # Frame key to output vector of CCMWaveforms
        SampleStep=50, # The number of steps between samples used for the initial baseline estimate
        NFramesForBaseline=10, # The number of frames to use for a baseline estimate
        ConstantFraction=0.05, # The fraction of the pulse height to use for its start time
        MinimumPulseHeight=1000, # The minimum pulse height to consider a NIM pulse to be present
        NIMLogicPulseSeriesMapName="NIMPulses", # Name for the output nim pulses map
    )

    # The NIM pulses are required if we want to restrict computing the BCM summary to only beam frames

    # Add the module that analyzes the beam current monitor waveform and produces summary information
    # All parameters are set to there defaults here, so the lines below are quivalent to:
    #   tray.AddModule("BeamCurrentMonitorSummary", "bcm_summary")
    tray.AddModule("BeamCurrentMonitorSummary", "bcm_summary",
        CCMGeometryName="CCMGeometry", # Frame key for CCMGeometry
        CCMWaveformsName="CCMWaveforms", # Frame key to output vector of CCMWaveforms
        TimeBeforePeak=2000.0, # Time in ns before the BCM peak to consider when computing the baseline and looking for the BCM start time
        ExpSmoothingTau=10.0, # Time constant in ns for exponential smoothing
        DerivativeThreshold=0.3, # Theshold below which derivativ is considered to be zero in ADC/ns
        OnlyBeamFrames=True, # Only run on frames with a beam NIM pulse?
        NIMPulsesName="NIMPulses", # Key for NIMLogicPulseSeriesMap
        CCMBCMSummaryName="BCMSummary", # Name for the output CCMBCMSummary
    )

    if "%" in args.output_file:
        tray.Add("I3MultiWriter", "writer", SizeLimit=(3*1024**3), FileName=args.output_file, MetadataStreams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus], CounterStream=icetray.I3Frame.DAQ, SizeCheckInterval=1000, NWorkers=8, Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
    else:
        tray.Add("I3Writer", "writer", FileName=args.output_file, Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])

    if args.num_events < 1:
        tray.Execute() # Process all frames
    else:
        tray.Execute(args.num_events + 1) # Number of frames to process is num_events DAQ frames plus one Geometry frame

