from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from I3Tray import I3Tray
from icecube.icetray import load
import random
import numpy as np
import bisect
import sys

np.set_printoptions(threshold=sys.maxsize)
np.seterr(all='raise')

class MuonInjector(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("MuonDistance", "The distance from the detector to place the muon in meters", 0.0)
        self.AddParameter("HalfPlaneSize", "The size of the injection half plane in meters", 1.0)
        self.seen_s_frame_ = False


    def Configure(self):
        self.muon_distance = self.GetParameter("MuonDistance")
        self.half_plane_size = self.GetParameter("HalfPlaneSize")
    '''
    def ThetaValue(self):
        x = 0.0
        y = 1000
        while y > x**2:
            x = np.random.uniform(0.0, 1.0)
            y = np.random.uniform(0.0, 1.0)
        th = np.arccos(x)
        return th

    def SmithDuller(self, Eu, zenith):
        Emu = np.power(10, Eu)
        Au = 2e9
        k = 2.66
        r = 0.76
        a = 0.0025 #GeV/gcm^{2}#
        y0 = 1000.0 #g/cm^{2}#
        bmu = 0.80
        c = 2.99e10  #cm/s#
        tau_muon = 2.2e-6 #s#
        rho_o = 0.001205  #g/cm^{3}#
        mass_muon = 0.1056 #GeV/c^{2}
        Bmu = (bmu*mass_muon*y0)/(c*tau_muon*rho_o)  #GeV #
        lamb = 120.0 #120 g/cm^{2} #
        b = 0.771
        mass_pio = 0.139 #GeV/c^{2}#
        tau_pio = 2.6e-8 #s#
        jpio = (y0*mass_pio)/(c*tau_pio*rho_o) # GeV #
        Epio = (Emu + a*y0*((1.0/np.cos(zenith))-0.100))/(r) #GeV#
        exp = (Bmu)/((r*Epio+100*a)*np.cos(zenith)) # 1 #
        in1 = a*((y0)/(np.cos(zenith))-100)/(r*Epio)  # 1  #
        in1 = np.amin([in1, np.full_like(in1, 1.0)], axis=0)
        ins = 0.100*np.cos(zenith)*(1.0-in1)
        Pmu = np.power(ins,exp) #1#
        p = (Au*(np.power(Epio,-k))*Pmu*lamb*b*jpio)/(Epio*np.cos(zenith)+b*jpio) # 1/GeV*str*cm^{2}*s #
        return p

    def sample_energy(self, cdf, x):
        r = np.random.uniform()
        cdf_idx = bisect.bisect_left(cdf, r)
        low_x = x[cdf_idx - 1]
        high_x = x[cdf_idx]
        low_cdf = cdf[cdf_idx - 1]
        high_cdf = cdf[cdf_idx]
        x_interp = (r - low_cdf) / (high_cdf - low_cdf) * (high_x - low_x) + low_x
        return x_interp
    '''
    def FillSimulationFrame(self, frame):
        # Fill the simulation frame with the event information
        config = dataclasses.I3MapStringDouble()
        #config[""] = 
        frame["MuonInjectorConfiguration"] = config

    def Simulation(self, frame):
        self.seen_s_frame_ = True
        self.FillSimulationFrame(frame)
        self.PushFrame(frame)

    def DAQ(self, frame):
        if not self.seen_s_frame_:
            sim_frame = icetray.I3Frame(icetray.I3Frame.Simulation)
            self.FillSimulationFrame(sim_frame)
            self.PushFrame(sim_frame)
            self.seen_s_frame_ = True

      # 1. choose the properties of the muon
        ## Zenith and Azimuth
        azimuth = np.random.uniform(0.0, 2*np.pi)
        #zenith = self.ThetaValue()
        zenith = np.random.uniform(0.0, np.pi)
        ## x, y and z position
        R = self.muon_distance + self.half_plane_size*np.random.uniform(0.0,1.0)  #radius of the hemisphere
        X = R*np.sin(zenith)*np.cos(azimuth)
        Y = R*np.sin(zenith)*np.sin(azimuth)
        Z = R*np.cos(zenith)
        ## random positions on the plane
        px = self.half_plane_size
        py = self.half_plane_size
        u = np.random.uniform(-px, px)
        v = np.random.uniform(-py, py)
        #########################  muon coordinates
        x = X + u*np.cos(zenith)*np.cos(azimuth)-v*np.sin(azimuth)
        y = Y + u*np.cos(zenith)*np.sin(azimuth)+v*np.cos(azimuth)
        z = Z - u*np.sin(zenith)

        ## energies
        # range of energies

        ### energy definition
        '''
        Emu_values = np.linspace(-5, -3, int(1e5)) # Emu in GeV (actually is the exponent of base 10 in the log scale E = 10^{Emu_values}) 
        pdf = self.SmithDuller(Emu_values,zenith)
        cdf = np.cumsum(pdf)
        cdf /= cdf[-1]

        exponent_kinetic_energy = self.sample_energy(cdf, Emu_values)
        kinetic_energy = np.power(10, exponent_kinetic_energy)
        mu_plus_minus_ratio = 1.2
        mu_threshold = mu_plus_minus_ratio / (1.0 + mu_plus_minus_ratio)
        r = np.random.uniform(0.0, 1.0)
        if r < mu_threshold:
            particle_type = dataclasses.I3Particle.MuPlus
        else:
            particle_type = dataclasses.I3Particle.MuMinus

        position = dataclasses.I3Position(x,y,z)
        direction = dataclasses.I3Direction(zenith, azimuth)
        '''
        # 2. Create and I3Particle object that represents the muon
        muon = dataclasses.I3Particle()
        muon.energy = 0.01 # in GeV
        muon.pos = dataclasses.I3Position(x,y,z)
        muon.dir = dataclasses.I3Direction(zenith, azimuth)
        muon.type = dataclasses.I3Particle.MuMinus

        #3. Create an I3MCTree object that represents all the particles we want to simulate
        monte_carlo_tree = dataclasses.I3MCTree()
        monte_carlo_tree.add_primary(muon)

        #4. put the I3MCTree
        frame["I3MCTree"] = monte_carlo_tree
        #5. Pass the frame to the next
        self.PushFrame(frame)



