
def run(files, output_file, N, idx):
    import glob
    import numpy as np
    print("Glob", files)
    fnames = [s for l in files for g in l for s in glob.glob(g)]
    print("Found these files:")
    for s in fnames:
        print(s)
    # exit(1)

    import I3Tray
    from icecube import icetray
    from icecube.icetray import load, I3Module, I3Frame
    from icecube import dataclasses
    load("CCMBinary", False)
    load("daqtools", False)
    load("dataclasses", False)

    sample_values = []
    wf_values = []
    first_samples = True
    first_wf = True

    def save(frame):
        nonlocal first_wf
        nonlocal first_samples
        if(frame.Has("BaselineEstimates")):
            n = len(frame.Get("BaselineEstimates"))
            if first_wf:
                for idx in range(n):
                    wf_values.append([])
            for idx in range(n):
                t = [
                    float(frame.Get("BaselineEstimates")[idx]),
                    float(frame.Get("BaselineEstimateVariances")[idx]),
                    float(frame.Get("FirstTriggerTime").value),
                ]
                wf_values[idx].append(t)
            first_wf = False

        if(frame.Has("BaselineSamples")):
            n = len(frame.Get("BaselineSamples"))
            if first_samples:
                for idx in range(n):
                    sample_values.append([])
            for idx in range(n):
                t = [
                    [float(x) for x in frame.Get("BaselineSamples")[idx]],
                    [float(x) for x in frame.Get("BaselineSampleVariances")[idx]],
                    [float(x) for x in frame.Get("BaselineSampleTimes")[idx]],
                ]
                sample_values[idx].extend(np.array(t).T)
            first_samples = False

    class selector(I3Module):
        def __init__(self, context):
            I3Module.__init__(self, context)
            self.AddParameter("N", "Number of frames to process", 100)
            self.AddParameter("Streams", "Streams to count", [I3Frame.DAQ])
            self.n = 0
            self.N = 0
            self.streams = [I3Frame.DAQ]

        def Configure(self):
            self.N = self.GetParameter("N")
            self.streams = self.GetParameter("Streams")
            if self.N <= 0:
                self.N = np.inf

        def Process(self):
            frame = self.PopFrame()
            if(frame.Stop in self.streams):
                self.n += 1
                if self.n > self.N:
                    self.RequestSuspension()
                    return
                else:
                    self.PushFrame(frame)
            else:
                self.PushFrame(frame)

    tray = I3Tray.I3Tray()
    tray.AddModule("I3Reader", "reader", FilenameList=fnames)
    tray.AddModule(selector, "selector", Streams=[icetray.I3Frame.DAQ], N=N)
    tray.AddModule(save, "save", Streams=[icetray.I3Frame.DAQ])
    tray.Execute()

    for i in range(len(sample_values)):
        np.savetxt(output_file + "/samples/test_samples_%04d.txt" % i, sample_values[i], delimiter=" ")
        np.savetxt(output_file + "/wf/test_wf_%04d.txt" % i, wf_values[i], delimiter=" ")


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

    parser.add_argument("--num-events",
            type=int,
            required=True,
            help="Number of events to process")

    parser.add_argument("--waveform-number",
            type=int,
            required=True,
            help="Index for waveform in CCMWaveformSeries")

    args = parser.parse_args()

    run(args.files, args.output_file, args.num_events, args.waveform_number)

