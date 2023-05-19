#!/usr/bin/env python3
import numpy as np
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.dataclasses import CCMOMGeo
from icecube.icetray import CCMTriggerKey, CCMPMTKey
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)

class SumWaveforms(I3ConditionalModule):
    # This constructor method is required, and is called once to instantiate the object
    # Registering which parameters the module has is required at this stage since parameters will be check immediately afterwards
    # The context argument is also required
    def __init__(self, context):
        # This line calls the constructor for the parent class and is also required
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("PMTTypes", "PMT types to sum",
                [
                    dataclasses.CCMOMGeo.CCM8inUncoated,
                    dataclasses.CCMOMGeo.CCM8inCoated
                ])
        self.AddParameter("OutputKey", "Key to save summed waveform to", "SummedWaveform")

    # This method is required and is where the module can get parameters passed to the I3Tray.Add method
    # It is called once by the I3Tray before processing begins
    def Configure(self):
        self.pmt_types = self.GetParameter("PMTTypes")
        self.output_name = self.GetParameter("OutputKey")

    # This method is called once for each "Geometry" frame in the data stream
    # The Geometry frame should be the first frame that the module receives
    # There is always one Geometry frame for a run, which contains the DAQ configuration and detector setup
    def Geometry(self, frame):
        # Get the geometry from the frame
        geo = frame['CCMGeometry']

        # We need the information about which PMT corresponds to which trigger
        self.pmt_channel_map = geo.pmt_channel_map
        # We need t know which trigger copy is relevant for each PMT
        self.trigger_copy_map = geo.trigger_copy_map

        # Gather information about which PMTs we will process,
        #  which trigger copy is relevant for that PMT,
        #  and which channel that PMT is on
        self.channels_to_sum = []
        for pmt_key, channel in self.pmt_channel_map:
            trigger_key = self.trigger_copy_map[pmt_key]
            pmt = geo.pmt_geo[pmt_key]
            if pmt.omtype in self.pmt_types:
                # Store PMT key, trigger key for corresponding board, and channel index
                self.channels_to_sum.append((pmt_key, trigger_key, channel))

        # Store all the trigger keys we will need to look at to determine waveform sizing and offsets
        self.unique_trigger_keys = list(sorted(set([trigger_key for _, trigger_key, _ in self.channels_to_sum])))
        # Here we push the frame to the next module so that any subsequent modules can process it
        self.PushFrame(frame)

    # Here we compute the sizes of each waveform that we intend to look at
    def get_sizes(self, waveforms):
        sizes = []
        for pmt_key, trigger_key, channel in self.channels_to_sum:
            size = len(waveforms[channel].waveform)
            sizes.append(size)
        sizes = np.unique(sizes)
        return sizes

    # Here we compute the extra elements required for the summed waveform
    # and the offsets of each waveform within the summed waveform,
    # or rather the offset for each board since waveforms on the same board share an offset
    def get_size_and_offsets(self, nim_pulses):
        offsets = {}
        for trigger_key in self.unique_trigger_keys:
            pulses = nim_pulses[trigger_key]
            if len(pulses) == 0:
                offsets[trigger_key] = None
                print("Warning! Could not fine NIM pulse for " + str(trigger_key))
                continue
            offset = int(pulses[0].time / 2.0)
            offsets[trigger_key] = offset
        offset_list = [offset for offset in offsets.values() if offset is not None]
        min_offset = np.amin(offset_list)
        max_offset = np.amax(offset_list)
        offsets = {trigger_key: max_offset - offset for trigger_key, offset in offsets.items()}

        return max_offset - min_offset, offsets

    # This method is called once for each "DAQ" frame in the data stream
    # Each DAQ frame corresponds to one continuous detector readout window
    def DAQ(self, frame):
        # Get the waveforms from the frame
        # Convert the C++ std::vector<CCMWaveformUInt16> object to a python list of CCMWaveformUInt16 objects
        waveforms = list(frame['CCMWaveforms'])

        # Check that the length of each waveform is the same
        # We want to skip events with missing data, and events with more than one trigger in the DAQ readout window
        sizes = self.get_sizes(waveforms)
        if len(sizes) > 1 or sizes[0] == 0:
            print("Skipping event because waveform sizes do not match or empty waveform")
            return
        size = sizes[0]

        # Grab the NIM pulses from the frame
        # This is of type std::map<CCMTriggerKey, std::vector<NIMLogicPulse>>
        nim_pulses = frame["NIMPulses"]

        # Compute the extra size we need in the summed waveform due to differing offsets on each board
        # Also compute the offsets of each waveform relative to the larger summed waveform
        extra_size, offsets = self.get_size_and_offsets(nim_pulses)

        summed_wf = np.zeros(size + extra_size)
        channels_summed = np.zeros(size + extra_size)
        for pmt_key, trigger_key, channel in self.channels_to_sum:
            offset = offsets[trigger_key]
            if offset is None:
                print("Skipping event because of missing board trigger copy NIM pulse")
                return

            # Convert the waveform to a numpy array by copying it
            # Invert the waveform so pulses are positive deviations from the baseline
            wf = -np.array(waveforms[channel].waveform).astype(np.int64)

            # Perform a crude baseline estimate using the median of the waveform and subtract this off
            wf -= np.int64(np.median(wf))

            # Find the position of this waveform in the (larger) summed waveform using the NIM pulses
            offset = offsets[trigger_key]

            # Add the waveform to the summed waveform
            #s = slice(offset, (offset - extra_size if offset < extra_size else None))
            summed_wf[offset:offset+size] += wf
            # Track how many waveforms have been added to each position in the summed waveform
            channels_summed[offset:offset+size] += 1

        # Normalize the waveform to a single channel
        summed_wf /= channels_summed

        # Store the summed waveform in the frame
        # Objects in the frame must be subclasses of I3FrameObject and cannot be raw python or numpy objects
        # I3VectorDouble is simply a version of std::vector<double> that derives from I3FrameObject
        frame[self.output_name] = dataclasses.I3VectorDouble(summed_wf)

        # Here we push the frame to the next module so that any subsequent modules can process it
        self.PushFrame(frame)

    # This optional method is called once by the I3Tray after all frames have been processed
    def Finish(self):
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
    parser.add_argument("--num-events",
            type=int,
            default=0,
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
    tray.Add(SumWaveforms, PMTTypes=[dataclasses.CCMOMGeo.CCM8inUncoated, dataclasses.CCMOMGeo.CCM8inCoated], OutputKey="InnerVolumeSummedWaveform")
    #tray.Add(SumWaveforms, PMTTypes=[dataclasses.CCMOMGeo.CCM8inCoated], OutputKey="Coated8inPMTSummedWaveform")
    #tray.Add(SumWaveforms, PMTTypes=[dataclasses.CCMOMGeo.CCM8inUncoated], OutputKey="Uncoated8inPMTSummedWaveform")
    #tray.Add(SumWaveforms, PMTTypes=[dataclasses.CCMOMGeo.CCM1in], OutputKey="VetoSummedWaveform")

    tray.Add("Dump") # Prints out the names of the objects in every frame
    if args.num_events < 1:
        tray.Execute() # Process all frames
    else:
        tray.Execute(args.num_events + 1) # Number of frames to process is num_events DAQ frames plus one Geometry frame

