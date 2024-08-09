from icecube import icetray, dataclasses
from icecube.icetray import I3ConditionalModule

class BCMCuts(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("CCMGeometryKey", "Name of the key for the CCMGeometry object", "CCMGeometry")
        self.AddParameter("CCMBCMSummaryKey", "Name of the key for the CCMBCMSummary object", "BCMSummary")
        self.AddParameter("NIMPulsesKey", "Name of the key for the NIMPulses object", "NIMPulses")
        self.AddParameter("MinStartTime", "Minimum start time of the BCM", -400)
        self.AddParameter("MaxStartTime", "Maximum start time of the BCM", -275)
        self.AddParameter("MinEndTime", "Minimum end time of the BCM", 150)
        self.AddParameter("MaxEndTime", "Maximum end time of the BCM", 275)
        self.AddParameter("MinDuration", "Minimum duration of the BCM", 450)
        self.AddParameter("MaxDuration", "Maximum duration of the BCM", 650)
        self.AddParameter("MinPeakTime", "Minimum peak time of the BCM", -200)
        self.AddParameter("MaxPeakTime", "Maximum peak time of the BCM", -100)
        self.AddParameter("MinPeakValue", "Minimum peak value of the BCM", 7000)
        self.AddParameter("MaxPeakValue", "Maximum peak value of the BCM", 12000)
        self.AddParameter("ApplyCut", "Apply the cut", True)
        self.AddParameter("OutputKey", "Name of the key for the output object", "PassedBCMCut")

    def Configure(self):
        self.geometry_key = self.GetParameter("CCMGeometryKey")
        self.bcm_key = self.GetParameter("CCMBCMSummaryKey")
        self.nim_key = self.GetParameter("NIMPulsesKey")
        self.min_start_time = self.GetParameter("MinStartTime")
        self.max_start_time = self.GetParameter("MaxStartTime")
        self.min_end_time = self.GetParameter("MinEndTime")
        self.max_end_time = self.GetParameter("MaxEndTime")
        self.min_duration = self.GetParameter("MinDuration")
        self.max_duration = self.GetParameter("MaxDuration")
        self.min_peak_time = self.GetParameter("MinPeakTime")
        self.max_peak_time = self.GetParameter("MaxPeakTime")
        self.min_peak_value = self.GetParameter("MinPeakValue")
        self.max_peak_value = self.GetParameter("MaxPeakValue")
        self.apply_cut = self.GetParameter("ApplyCut")
        self.output_key = self.GetParameter("OutputKey")

    def Geometry(self, frame):
        bcm_key = icetray.CCMPMTKey(10, 1, 0)
        geo = frame[self.geometry_key]
        self.bcm_trigger_key = geo.trigger_copy_map[bcm_key]
        self.PushFrame(frame)

    def DAQ(self, frame):
        if not frame.Has(self.bcm_key):
            self.PushFrame(frame)
            return

        if not frame.Has(self.nim_key):
            icetray.logging.log_fatal(f"NIMPulses object \"{self.nim_key}\" not found in frame")
            self.PushFrame(frame)
            return

        bcm = frame[self.bcm_key]
        nim_pulses = frame[self.nim_key]

        trigger_nim_pulses = nim_pulses[self.bcm_trigger_key]
        if len(trigger_nim_pulses) == 0:
            self.PushFrame(frame)
            return
        trigger_time = trigger_nim_pulses[0].time

        start_time = bcm.start_time - trigger_time
        end_time = bcm.end_time - trigger_time
        peak_time = bcm.peak_time - trigger_time
        peak_value = bcm.peak_value
        duration = end_time - start_time

        bad_bcm = (
                start_time < self.min_start_time or
                start_time > self.max_start_time or
                end_time < self.min_end_time or
                end_time > self.max_end_time or
                duration < self.min_duration or
                duration > self.max_duration or
                peak_time < self.min_peak_time or
                peak_time > self.max_peak_time or
                peak_value < self.min_peak_value or
                peak_value > self.max_peak_value
        )

        if self.apply_cut and bad_bcm:
            return

        passed_bcm = dataclasses.I3UInt32(not bad_bcm)
        frame[self.output_key] = passed_bcm

        self.PushFrame(frame)

