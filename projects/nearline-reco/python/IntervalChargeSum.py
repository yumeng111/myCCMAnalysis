import icecube.icetray

@icecube.icetray.traysegment
def IntervalChargeSum(tray, name,
        CCMGeometryName=None,
        PMTTypes=None,
        TimeWindows=None,
        InputPulsesMaskName=None,
        InputRawPulsesName=None,
        InputEventPrefix=None,
        OutputPrefix=None,
        ProcessDAQFrames=True,
        ProcessPhysicsFrames=True,
        ):
    kwargs = {
            "CCMGeometryName": CCMGeometryName,
            "PMTTypes": PMTTypes,
            "TimeWindows": TimeWindows,
            "InputPulsesMaskName": InputPulsesMaskName,
            "InputRawPulsesName": InputRawPulsesName,
            "InputEventPrefix": InputEventPrefix,
            "OutputPrefix": OutputPrefix,
            }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    if not ProcessDAQFrames and not ProcessPhysicsFrames:
        icetray.logging.log_fatal("You must process at least one frame type")
    if ProcessDAQFrames:
        tray.AddModule("IntervalChargeSumQ", name + "_DAQ", **kwargs)
    if ProcessPhysicsFrames:
        tray.AddModule("IntervalChargeSumP", name + "_Physics", **kwargs)

