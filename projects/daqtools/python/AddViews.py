from icecube import dataclasses
from icecube.icetray import I3Module

class AddViews(I3Module): # module definitions are classes that inherit from I3Module or I3ConditionalModule
    def __init__(self, context): # __init__ is the constructor
        I3Module.__init__(self, context) # call the constructor of the base class (this boilerplate)

    def DAQ(self, frame): # DAQ is a method that is called for every DAQ frame
        if "TriggerTimePulses" not in frame:
            frame["TriggerTimePulses"] = dataclasses.CCMRecoPulseSeriesMapApplySPECalPlusTriggerTime("WavedeformPulses", "CCMCalibration", "NIMPulses", "CCMGeometry")
        if "BCMSummary" in frame:
            if "BeamTimePulses" not in frame:
                frame["BeamTimePulses"] = dataclasses.CCMRecoPulseSeriesMapApplySPECalPlusBeamTime("WavedeformPulses", "CCMCalibration", "NIMPulses", "CCMGeometry", "BCMSummary")
        self.PushFrame(frame)
