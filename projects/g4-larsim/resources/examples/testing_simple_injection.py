#!/usr/bin/env python3
import argparse
from icecube.icetray import load
from I3Tray import I3Tray
from icecube import CCMBinary, icetray, dataio, dataclasses, phys_services, daqtools, simclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.icetray import I3Units
from icecube.dataclasses import CCMOMGeo, I3VectorDouble
from icecube.icetray import CCMTriggerKey, CCMPMTKey
load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)
load("g4-larsim")

parser = argparse.ArgumentParser()
parser.add_argument("-p", dest="particle", default="e-", help="particle type")
parser.add_argument("-e", dest="energy", default=1.0*I3Units.MeV, help="particle energy")
parser.add_argument("-l", dest="location", default=I3VectorDouble([0.0, 0.0, 0.0]), help="particle location")
parser.add_argument("-d", dest="direction", default=I3VectorDouble([0.0, 0.0, 1.0]), help="particle direction")
parser.add_argument("-n", dest="Nevents", default=1, help = "number of events to process")
parser.add_argument("-pmt", dest="PMTSDStatus", default=True, help = "true for recording hits in PMT")
parser.add_argument("-lar", dest="LArSDStatus", default=True, help = "true for recording hits in LAr")
args = parser.parse_args()

tray = I3Tray()
tray.AddModule("I3InfiniteSource", "source", Stream=icetray.I3Frame.DAQ)
tray.AddService("CCMSimpleInjectorFactory", "injector", ParticleType=args.particle, ParticleEnergy=args.energy, ParticleLocation=args.location, ParticleDirection=args.direction)
tray.AddService("CCM200ResponseFactory", "response", PMTSDStatus=args.PMTSDStatus, LArSDStatus=args.LArSDStatus)
tray.AddModule("CCMSimulator", "ccm simulator", InjectorServiceName="injector", ResponseServiceName="response")
tray.AddModule("I3Writer", "i3-writer", filename = "g4-larsim-test.i3", streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Simulation])
tray.Execute(args.Nevents)

