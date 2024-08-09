import collections
from icecube.icetray import (
    CCMTriggerKey,
    I3ConditionalModule,
)
from icecube import dataclasses

class TriggerType(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("NIMPulsesKey", "Name of the NIMPulses key", "NIMPulses")

    def Configure(self):
        self.nim_pulses_key = self.GetParameter("NIMPulsesKey")
        self.possible_trigger_types = set(
            [
                CCMTriggerKey.TriggerType.StrobeTrigger,
                CCMTriggerKey.TriggerType.LEDTopTrigger,
                CCMTriggerKey.TriggerType.LEDBottomTrigger,
                CCMTriggerKey.TriggerType.BeamTrigger,
                CCMTriggerKey.TriggerType.CosmicTrigger,
                CCMTriggerKey.TriggerType.LaserTrigger,
            ]
        )

    def DAQ(self, frame):
        if not frame.Has(self.nim_pulses_key):
            self.PushFrame(frame)
            return

        seen_triggers = collections.defaultdict(int)

        nim_pulses = frame[self.nim_pulses_key]

        for trigger_key, pulses in nim_pulses.items():
            if trigger_key.type in self.possible_trigger_types and len(pulses) > 0:
                seen_triggers[trigger_key.type] += 1

        if len(seen_triggers) == 1:
            frame["TriggerType"] = dataclasses.I3Int32(
                int(list(seen_triggers.keys())[0])
            )
        else:
            frame["TriggerType"] = dataclasses.I3Int32(
                int(CCMTriggerKey.TriggerType.UnknownType)
            )

        self.PushFrame(frame)
