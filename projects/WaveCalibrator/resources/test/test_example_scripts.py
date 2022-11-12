#!/usr/bin/env python3

import os
import unittest
from tempfile import NamedTemporaryFile
from icecube.icetray import I3Test

class TestSimpleExample(I3Test.TestExampleScripts):

    project_name = "WaveCalibrator"

    def test_simple_example(self):
        '''
        Test that runs the one example script in this project.
        '''
        gcd_file = self.I3_TESTDATA + "/GCD/GeoCalibDetectorStatus_2013.56429_V1.i3.gz"
        input_file = self.I3_TESTDATA + "/sim/Level2_IC86.2011_corsika.010281.001664.00.i3.bz2"
        with NamedTemporaryFile(dir=os.getcwd(), suffix=".gz") as output_file:
            self.run_example('simple_example.py', gcd_file, input_file, output_file.name)

unittest.main()
