import icecube.icetray

@icecube.icetray.traysegment
def IntervalChargeSum(tray, name,
        CCMGeometryName=None,
        PMTTypes=None,
        TimeWindows=None,
        InputRawPulsesName=None,
        InputEventPrefix=None,
        OutputPrefix=None,
        BadKeysNames=None,
        ProcessDAQFrames=True,
        ProcessPhysicsFrames=True,
        If=None):
    kwargs = {
            "CCMGeometryName": CCMGeometryName,
            "PMTTypes": PMTTypes,
            "TimeWindows": TimeWindows,
            "InputRawPulsesName": InputRawPulsesName,
            "InputEventPrefix": InputEventPrefix,
            "OutputPrefix": OutputPrefix,
            "BadKeysNames": BadKeysNames,
            "If": If,
            }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    if not ProcessDAQFrames and not ProcessPhysicsFrames:
        icecube.icetray.logging.log_fatal("You must process at least one frame type")
    if ProcessDAQFrames:
        tray.AddModule("IntervalChargeSumQ", name + "_DAQ", **kwargs)
    if ProcessPhysicsFrames:
        tray.AddModule("IntervalChargeSumP", name + "_Physics", **kwargs)

