from icecube.icetray import (
    CCMTriggerKey,
    I3ConditionalModule,
    I3Frame,
    traysegment,
)

class TriggerTypeFilterModule(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("DesiredTriggerTypes", "Trigger type you want to select for and continue processing data", [])
        self.AddParameter("FrameTypesToFilter", "The frame types to filter on", [I3Frame.DAQ, I3Frame.Physics])
        self.n_frames = 0

    def Configure(self):
        self.desired_trigger_types = self.GetParameter("DesiredTriggerTypes")
        self.frame_types = self.GetParameter("FrameTypesToFilter")


    def IsCorrectTriggerType(self, frame):
        ### let's get the nim pulse
        nim_pulses = frame["NIMPulses"]
        for trigger_key, pulses in nim_pulses.items():
            if trigger_key.type in self.desired_trigger_types and len(pulses) > 0:
                return True
        return False

    def Process(self):
        frame = self.PopFrame()
        if frame.Stop in self.frame_types:
            if self.IsCorrectTriggerType(frame):
                self.PushFrame(frame)
            else:
                return
        else:
            self.PushFrame(frame)

    def Finish(self):
        pass

@traysegment
def TriggerTypeFilter(tray, name, DesiredTriggerTypes=None, FrameTypesToFilter=None):
    if DesiredTriggerTypes is None:
        return
    if FrameTypesToFilter is None:
        FrameTypesToFilter = [I3Frame.DAQ, I3Frame.Physics]

    trigger_types_map = {
        "strobe": [CCMTriggerKey.TriggerType.StrobeTrigger],
        "ledtop": [CCMTriggerKey.TriggerType.LEDTopTrigger],
        "ledbottom": [CCMTriggerKey.TriggerType.LEDBottomTrigger],
        "led": [CCMTriggerKey.TriggerType.LEDTopTrigger, CCMTriggerKey.TriggerType.LEDBottomTrigger],
        "beam": [CCMTriggerKey.TriggerType.BeamTrigger],
        "cosmic": [CCMTriggerKey.TriggerType.CosmicTrigger],
        "laser": [CCMTriggerKey.TriggerType.LaserTrigger],
        "all": [CCMTriggerKey.TriggerType.StrobeTrigger, CCMTriggerKey.TriggerType.LEDTopTrigger, CCMTriggerKey.TriggerType.LEDBottomTrigger, CCMTriggerKey.TriggerType.BeamTrigger, CCMTriggerKey.TriggerType.CosmicTrigger, CCMTriggerKey.TriggerType.LaserTrigger],
    }

    trigger_types = []
    for s in DesiredTriggerTypes:
        if type(s) == CCMTriggerKey.TriggerType:
            trigger_types.append(s)
        elif s.lower() in trigger_types_map:
            trigger_types.extend(trigger_types_map[s.lower()])
    DesiredTriggerTypes = list(sorted(list(set(trigger_types))))

    tray.AddModule(TriggerTypeFilterModule, "TriggerTypeFilter", DesiredTriggerTypes=DesiredTriggerTypes, FrameTypesToFilter=FrameTypesToFilter)

