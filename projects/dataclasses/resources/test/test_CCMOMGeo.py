#!/usr/bin/env python3

import unittest
from icecube import dataclasses
from I3Tray import I3Units

class TestCCMModuleGeo(unittest.TestCase):

    def test_CCMOMGeo(self):
        omgeo1 = dataclasses.CCMOMGeo()
        omgeo1.omtype = dataclasses.CCMOMGeo.OMType.IceCube
        omgeo1.position = dataclasses.I3Position(0,0,0)
        omgeo1.orientation = dataclasses.I3Orientation(1,0,0,0,1,0)
        omgeo1.area = 42.

        omgeo2 = dataclasses.CCMOMGeo()
        omgeo2.omtype = dataclasses.CCMOMGeo.OMType.IceCube
        omgeo2.position = dataclasses.I3Position(0,0,0)
        omgeo2.orientation = dataclasses.I3Orientation(1,0,0,0,1,0)
        omgeo2.area = 42.

        self.assertEqual(omgeo1, omgeo2, "these should be the same.")


unittest.main()
