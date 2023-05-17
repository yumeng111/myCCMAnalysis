#!/usr/bin/env python3
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)

class PrintNIM(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)

    def Configure(self):
        pass

    def Geometry(self, frame):
        pass

    def DAQ(self, frame):
        print("Printing frame[\"NIMPulses\"]")
        print(frame["NIMPulses"])
        print()

    def Finish(self):
        pass

if __name__ == "__main__":

    import argparse
    # Input arguments
    parser = argparse.ArgumentParser(
            prog = "NIM pulse printer",
            description = "Analyze trigger channels and print NIM pulses",
            )
    parser.add_argument("--files",
            type=str,
            nargs="+",
            action="append",
            default=[],
            help="Files to Process")
    parser.add_argument("--num-events",
            type=int,
            default=100,
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

    tray.Add("Dump") # Prints out the names of the objects in every frame
    tray.Add(PrintNIM, "print_nim") # Prints out the NIM pulses in each DAQ frame
    tray.Execute(args.num_events + 1) # Number of frames to process is num_events DAQ frames plus one Geometry frame

