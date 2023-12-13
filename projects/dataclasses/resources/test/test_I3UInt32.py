#!/usr/bin/env python3

import unittest
from icecube.dataclasses import I3UInt32

class TestI3UInt32(unittest.TestCase):
    def test_bool(self):
        self.assertFalse(bool(I3UInt32(0)), "this should be false.")
        self.assertTrue(bool(I3UInt32(1)), "this should be false.")        

if __name__ == "__main__":
    unittest.main()
