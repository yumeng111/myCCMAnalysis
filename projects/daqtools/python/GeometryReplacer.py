#!/usr/bin/env python3
import os, numpy, glob
from icecube import icetray, dataclasses, dataio

class GeometryReplacer(icetray.I3Module):
    """
    Replaces information in the geometry frame
    """
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("GeometryFile", "Geometry file to use as replacement", None)
        self.AddParameter("Keys", "Keys to replace", ["CCMGeometry", "CCMDAQConfig"])

    def Configure(self):
        self.geometry_file = self.GetParameter("GeometryFile")
        self.keys = self.GetParameter("Keys")
        i3file = dataio.I3File(self.geometry_file)
        self.geo_frame = None
        while i3file.more():
            geo_frame = i3file.pop_frame()
            if geo_frame.Stop == icetray.I3Frame.Geometry:
                self.geo_frame = geo_frame
                break
        if self.geo_frame is None:
            raise RuntimeError("Could not find geometry frame in supplied file")
        for key in self.keys:
            if not self.geo_frame.Has(key):
                raise RuntimeError("Supplied geometry frame does not have specified key")


    def Geometry(self, frame):
        for key in self.keys:
            if frame.Has(key):
                del frame[key]
            frame[key] = self.geo_frame[key]
        self.PushFrame(frame)

