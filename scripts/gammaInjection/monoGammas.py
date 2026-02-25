from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
#change
#from I3Tray import I3Tray
#from icecube.icetray import load
#import random
import numpy as np
import bisect
import sys

#change
#np.set_printoptions(threshold=sys.maxsize)
#np.seterr(all='raise')

# class MuonInjector(I3ConditionalModule):
#change
class GammaInjector(I3ConditionalModule):
    def __init__(self, context):
        # I3ConditionalModule.__init__(self, context)
        #change
        super().__init__(context)

        # change:geometry knobs for plane source injection
        #self.AddParameter("MuonDistance", "The distance from the detector to place the muon in meters", 0.0)
        self.AddParameter("GammaDistance", "The distance from the detector to place the gamma in meters", 0.0)
        self.AddParameter("HalfPlaneSize", "The size of the injection half plane in meters", 1.0)
        

        #change
        self.AddParameter("GammaEnergyMeV", "Gamma energy in MeV", 10.0)
        self.AddParameter("UsePointSource", "If True, inject from point; else use surface geometry", True)
        self.AddParameter("DecayX", "Point source x (meters)", 0.0)
        self.AddParameter("DecayY", "Point source y (meters)", 0.0)
        self.AddParameter("DecayZ", "Point source z (meters)", 0.0)
        self.AddParameter("ArgonInjectionRadiusM", "fOuterLAr radius (m)", 1.20)
        self.AddParameter("ArgonInjectionHalfHeightM", "fOuterLAr half-height (m)", 1.15)
        self.AddParameter("Isotropic", "If True, sample isotropic direction", False)

        self.seen_s_frame_ = False

    def Configure(self):
        # self.muon_distance = self.GetParameter("MuonDistance")
        #change
        self.gamma_distance = float(self.GetParameter("GammaDistance"))
        self.half_plane_size = float(self.GetParameter("HalfPlaneSize"))
        #change
        self.gamma_energy_mev = float(self.GetParameter("GammaEnergyMeV"))
        self.use_point = bool(self.GetParameter("UsePointSource"))
        self.decay_x = float(self.GetParameter("DecayX"))
        self.decay_y = float(self.GetParameter("DecayY"))
        self.decay_z = float(self.GetParameter("DecayZ"))
        self.argon_injection_radius = float(self.GetParameter("ArgonInjectionRadiusM"))
        self.argon_injection_half_height = float(self.GetParameter("ArgonInjectionHalfHeightM"))
        self.isotropic = bool(self.GetParameter("Isotropic"))

    def FillSimulationFrame(self, frame):
        # Fill the simulation frame with the event information
        config = dataclasses.I3MapStringDouble()
        #config[""] = 
        # frame["MuonInjectorConfiguration"] = config
        #change
        config["GammaEnergyMeV"] = self.gamma_energy_mev
        config["UsePointSource"] = 1.0 if self.use_point else 0.0
        config["DecayX"] = self.decay_x
        config["DecayY"] = self.decay_y
        config["DecayZ"] = self.decay_z
        config["Isotropic"] = 1.0 if self.isotropic else 0.0
        frame["GammaInjectorConfiguration"] = config

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

      # 1. choose the properties of the muon/gamma
        #change (fix isotropic direction sampling?: uniform zenithis not isotropic, but will oversample the poles)
        if self.isotropic:
            u = np.random.uniform(-1.0, 1.0)
            zenith = np.arccos(u)
            azimuth = np.random.uniform(0.0, 2.0 * np.pi)
        else:
        #   zenith = self.ThetaValue()
        #   zenith = 0.0
        #   azimuth = 0.0
           azimuth = np.random.uniform(0.0, 2*np.pi)
           zenith = np.random.uniform(0.0, np.pi)
 
        if self.use_point:
            #x, y, z = self.decay_x, self.decay_y, self.decay_z
            pos_phi = np.random.uniform(0.0, 2*np.pi)
            pos_r   = np.random.uniform(0.0, self.argon_injection_radius)
            pos_z   = np.random.uniform(-self.argon_injection_half_height, self.argon_injection_half_height)

            x = pos_r * np.cos(pos_phi)
            y = pos_r * np.sin(pos_phi)
            z = pos_z
        
        else:
            R = self.gamma_distance + self.half_plane_size * np.random.uniform(0.0, 1.0) #radius of the hemisphere
            X = R * np.sin(zenith) * np.cos(azimuth)
            Y = R * np.sin(zenith) * np.sin(azimuth)
            Z = R * np.cos(zenith)

            # ## random positions on the plane
            px = self.half_plane_size
            py = self.half_plane_size
            uu = np.random.uniform(-px, px)
            vv = np.random.uniform(-py, py)

            # #########################  muon coordinates
            x = X + uu * np.cos(zenith) * np.cos(azimuth) - vv * np.sin(azimuth)
            y = Y + uu * np.cos(zenith) * np.sin(azimuth) + vv * np.cos(azimuth)
            z = Z - uu * np.sin(zenith)

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

        #change (use parameter and convert MeV toGeV)
        energy_gev = self.gamma_energy_mev * 1e-3

        # 2. Create and I3Particle object that represents the muon/gamma
        #change
        # muon = dataclasses.I3Particle()
        # muon.energy = 0.01 # in GeV
        # muon.pos = dataclasses.I3Position(x,y,z)
        # muon.dir = dataclasses.I3Direction(zenith, azimuth)
        # muon.type = dataclasses.I3Particle.MuMinus
        
        gamma = dataclasses.I3Particle()
        gamma.energy = energy_gev
        gamma.pos = dataclasses.I3Position(x, y, z)
        gamma.dir = dataclasses.I3Direction(zenith, azimuth)
        gamma.type = dataclasses.I3Particle.Gamma
        #change
        gamma.time = 0.0

        #3. Create an I3MCTree object that represents all the particles we want to simulate
        # monte_carlo_tree = dataclasses.I3MCTree()
        # monte_carlo_tree.add_primary(muon)
        #change
        monte_carlo_tree = dataclasses.I3MCTree()
        monte_carlo_tree.add_primary(gamma)

        #4. put the I3MCTree
        frame["I3MCTree"] = monte_carlo_tree
        #5. Pass the frame to the next
        self.PushFrame(frame)
