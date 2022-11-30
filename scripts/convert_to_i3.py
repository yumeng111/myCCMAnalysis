
def convert(files, output_file):
    import glob
    print("Glob", files)
    fnames = [s for l in files for g in l for s in glob.glob(g)]
    print("Found these files:")
    for s in fnames:
        print(s)
    # exit(1)

    import I3Tray
    from icecube.icetray import load
    load("CCMBinary", False)

    tray = I3Tray.I3Tray()
    tray.AddModule("CCMBinaryReader", "reader", FilenameList=fnames)
    tray.AddModule("I3Writer", "writer", Filename=output_file)
    tray.Execute()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
            prog = "ConvertToI3",
            description = "Converts CCM versioned binary files to i3 files",
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

    convert(args.files, args.output_file)

