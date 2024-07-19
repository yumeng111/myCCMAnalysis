from icecube import icetray, dataclasses
from icecube.icetray import I3ConditionalModule

class ChargeRatio(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("TimeWindowNumerator", "Time window in ns for numerator", 20)
        self.AddParameter("TimeWindowDenominator", "Time window in ns for denominator", 40)

    def Configure(self):
        self.time_window_num = int(self.GetParameter("TimeWindowNumerator"))
        self.time_window_den = int(self.GetParameter("TimeWindowDenominator"))

    def Physics(self, frame):
        coated_num = frame[f"CoatedEventCharge{self.time_window_num}NS"].value
        uncoated_num = frame[f"UncoatedEventCharge{self.time_window_num}NS"].value

        coated_den = frame[f"CoatedEventCharge{self.time_window_den}NS"].value
        uncoated_den = frame[f"UncoatedEventCharge{self.time_window_den}NS"].value

        tot_num = coated_num + uncoated_num
        tot_den = coated_den + uncoated_den

        if tot_den == 0:
            prop = 0
        else:
            prop = tot_num / tot_den

        frame[f"ChargeRatio{self.time_window_num}v{self.time_window_den}"] = dataclasses.I3Double(prop)

        self.PushFrame(frame)

