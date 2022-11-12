from icecube.icetray import I3ConditionalModule
from icecube import dataclasses

class DOMFlagger(I3ConditionalModule):
	def __init__(self, context):
		I3ConditionalModule.__init__(self, context)
		self.AddParameter("Key", "Name of the DOMLaunches/Waveforms/Pulses in the frame", None)
		self.AddParameter("Flag", "Name of the flag vector to put in the frame", None)
		self.AddParameter("Condition", "Function that returns True if a series of DOMLaunches/Waveforms/Pulses should be transmitted in raw form, False otherwise", lambda vec: False)
	
	def Configure(self):
		self.key = self.GetParameter("Key")
		self.flag = self.GetParameter("Flag")
		self.callable = self.GetParameter("Condition")
	
	def DAQ(self, frame):
		
		if not self.key in frame:
			self.PushFrame(frame)
			return
		
		frameobject = frame[self.key]
		
		flags = dataclasses.I3VectorOMKey()
		for om, vector in frameobject.items():
			if self.callable(vector):
				flags.append(om)
		
		if len(flags) > 0:
			frame[self.flag] = flags
		
		self.PushFrame(frame)
		return

		
