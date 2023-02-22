
def run(files, output_file, sample_length):
    import glob
    print("Glob", files)
    fnames = [s for l in files for g in l for s in glob.glob(g)]
    print("Found these files:")
    for s in fnames:
        print(s)
    # exit(1)

    import I3Tray
    from icecube import icetray
    from icecube.icetray import load
    from icecube import dataclasses
    load("CCMBinary", False)
    load("daqtools", False)

    tray = I3Tray.I3Tray()
    tray.AddModule("I3Reader", "reader", FilenameList=fnames)
    tray.AddModule("CCMBaselineAnalyzer", MinimumSampleLength=sample_length)
    tray.AddModule("Keep", Keys=[
        "BaselineEstimateVariances",
        "BaselineEstimates",
        "BaselineSampleTimes",
        "BaselineSampleVariances",
        "BaselineSamples",
        "FirstTriggerTime",
        "TriggerTimes",
    ])
    if "%" in output_file:
        tray.AddModule("I3MultiWriter", "writer", SizeLimit=int(1.0*1024**3), FileName=output_file, MetadataStreams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus], CounterStream=icetray.I3Frame.DAQ, SizeCheckInterval=1000, NWorkers=8, Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
    else:
        tray.AddModule("I3Writer", "writer", FileName=output_file, Streams=[icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus, icetray.I3Frame.DAQ, icetray.I3Frame.Physics])
    tray.Execute()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            prog = "Analyze baseline",
            description = "Analyze PMT baselines for a run",
    )

    parser.add_argument("--files",
            type=str,
            nargs="+",
            action="append",
            default=[],
            required=True,
            help="Files to process")
    parser.add_argument("--output-file",
            type=str,
            required=True,
            help="Output file name")
    parser.add_argument("--sample-length",
            type=int,
            required=True,
            help="Minimum number of samples")

    args = parser.parse_args()

    run(args.files, args.output_file, args.sample_length)

