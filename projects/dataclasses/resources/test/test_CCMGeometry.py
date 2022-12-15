#!/usr/bin/env python3

import unittest
from icecube import dataclasses
from I3Tray import I3Units

class TestCCMGeometry(unittest.TestCase):

    def test_CCMGeometry_equality(self):
        geo1 = dataclasses.CCMGeometry()
        geo1.omgeo = dataclasses.I3OMGeoMap()
        geo1.start_time = dataclasses.I3Time()
        geo1.end_time = dataclasses.I3Time()

        geo2 = dataclasses.CCMGeometry()
        geo2.omgeo = dataclasses.I3OMGeoMap()
        geo2.start_time = dataclasses.I3Time()
        geo2.end_time = dataclasses.I3Time()

        self.assertEqual(geo1, geo2, "these should be the same.")


unittest.main()
