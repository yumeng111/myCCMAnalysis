
def merge(files, output_file):
    import glob
    print("Glob", files)
    fnames = [[s for g in l for s in glob.glob(g)] for l in files]
    print("Found these files:")
    for s in fnames:
        print(s)
    # exit(1)

    import I3Tray
    from icecube import icetray
    from icecube.icetray import load
    load("CCMBinary", False)

    tray = I3Tray.I3Tray()
    tray.AddModule("MergedSource", "reader", FileLists=fnames)
    if "%" in output_file:
        tray.AddModule("I3MultiWriter", "writer", SizeLimit=1, FileName=output_file, MetadataStreams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus], CounterStream=icetray.I3Frame.DAQ, SizeCheckInterval=1000, NWorkers=6)
    else:
        tray.AddModule("I3Writer", "writer", FileName=output_file)
    tray.Execute()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            prog = "Merge DAQ streams",
            description = "Merges CCM i3 files from multiple DAQs",
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

    args = parser.parse_args()

    merge(args.files, args.output_file)

