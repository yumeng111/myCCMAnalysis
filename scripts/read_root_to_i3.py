
def merge(files, output_file):
    import glob
    print("Glob", files)
    fnames = [s for g in files for s in sorted(glob.glob(g))]
    print("Found these files:")
    for s in fnames:
        print(s)
    # exit(1)

    import I3Tray
    from icecube import icetray
    from icecube.icetray import load
    load("CCMBinary", False)
    load("CCMDataclasses", False);
    load("CCMIO", False)

    tray = I3Tray.I3Tray()
    tray.AddModule("CCMI3RootReader", "reader", FilenameList=fnames, KeysToLoad=["pulses", "events", "accumWaveform", "mcTruth"])
    if "%" in output_file:
        tray.AddModule("I3MultiWriter", "writer", SizeLimit=(3*1024**3), FileName=output_file, MetadataStreams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus], CounterStream=icetray.I3Frame.DAQ, SizeCheckInterval=1000, NWorkers=8)
    else:
        tray.AddModule("I3Writer", "writer", FileName=output_file)
    tray.Execute()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            prog = "Read CCM root file",
            description = "Reads data from a CCM root file and saves it to the i3 format",
    )

    parser.add_argument("--files",
            type=str,
            nargs="+",
            default=[],
            required=True,
            help="Files to process")

    parser.add_argument("--output-file",
            type=str,
            required=True,
            help="Output file name")

    args = parser.parse_args()

    merge(args.files, args.output_file)

