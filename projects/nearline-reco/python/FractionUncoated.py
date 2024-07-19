import numpy as np

from icecube import icetray, dataclasses
from icecube.icetray import I3ConditionalModule

class FractionUncoated(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("TimeWindow", "Time window in ns", 10)
    def Configure(self):
        self.time_window = self.GetParameter("TimeWindow")
    def Physics(self, frame):
        uncoated = frame[f"UncoatedEventCharge{self.time_window}NS"].value
        coated = frame[f"CoatedEventCharge{self.time_window}NS"].value
        s = uncoated + coated
        if np.isnan(s) or s == 0:
            frame[f"FracUncoated{self.time_window}NS"] = dataclasses.I3Double(0)
        else:
            frame[f"FracUncoated{self.time_window}NS"] = dataclasses.I3Double(uncoated / s)

        self.PushFrame(frame)

