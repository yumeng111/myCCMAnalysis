
def run(files):
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
    tray.AddModule("CCMBaselineAnalyzer")
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

    args = parser.parse_args()

    run(args.files)

