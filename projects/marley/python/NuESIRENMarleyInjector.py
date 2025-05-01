import siren
from siren import utilities
import pickle
import argparse
import os
import numpy as np
from icecube import icetray, dataclasses
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule

class NuESIRENMarleyInjector(I3ConditionalModule):
    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter("Target", "0 for upper tungsten target, 1 for lower tungsten target", 0)
        self.AddParameter("EventsToInject", "Number of events to inject", 10000)
        self.seen_s_frame_ = False

    def Configure(self):
        self.target = self.GetParameter("Target")
        self.events_to_inject = self.GetParameter("EventsToInject")

        #Load experiment and detector
        self.experiment = "CCM"
        self.detector_model = utilities.load_detector(self.experiment)
        self.fiducial_volume = utilities.get_fiducial_volume(self.experiment)

        #Define primary particle
        self.primary_type = siren.dataclasses.Particle.ParticleType.NuE

        #Load Marley process
        self.primary_processes, self.secondary_processes = siren.resources.load_processes(
                    "MarleyCrossSection",
                    primary_types=[self.primary_type],
                    process_types=["CC"],
                )

        #Primary distributions
        self.mass_ddist = siren.distributions.PrimaryMass(0) #mass
        self.edist = siren.distributions.PiDARNuEDistribution() #energy

        # Cone direction distribution from tungsten target
        #self.opening_angle = np.arcsin(1.21 / 23.0)
        self.opening_angle = np.arcsin(2.0 / 23.0)

        # Fixed tungsten target
        if self.target == 0:
            self.target_origin = siren.math.Vector3D(0, 0, 0.1375)
        else:
            self.target_origin = siren.math.Vector3D(0, 0, -0.241)

        self.detector_origin = siren.math.Vector3D(23, 0, -0.65)
        self.target_dir = self.detector_origin - self.target_origin
        self.target_dir.normalize()
        self.target_inj_ddist = siren.distributions.Cone(self.target_dir, self.opening_angle)

        # Physical direction distribution
        self.phys_ddist = siren.distributions.IsotropicDirection()

        #Position distribution
        self.max_dist = 25.0
        # Position distribution: use fixed target distribution that considers physical volume of the target
        if self.target == 0:
            # upper target
            self.target_cylinder = siren.geometry.Cylinder(siren.geometry.Placement(self.target_origin - self.detector_origin), 0.05, 0.0, 0.091)
        else:
            # lower target
            self.target_cylinder = siren.geometry.Cylinder(siren.geometry.Placement(self.target_origin - self.detector_origin), 0.05, 0.0, 0.298)

        self.target_pos_dist = siren.distributions.FixedTargetPositionDistribution(self.target_cylinder, self.fiducial_volume, self.max_dist)

        self.target_area_dist = siren.distributions.FixedTargetAreaDistribution(self.target_cylinder)

        # Primary injection distributions
        self.primary_injection_distributions = [
            self.mass_ddist, # Mass distribution
            self.edist, # Energy distribution
            self.target_inj_ddist, # Direction distribution (cone)
            self.target_pos_dist, # Position distribution (fixed target position)
            ]

        self.primary_physical_distributions = [
            self.edist,  # Energy distribution
            self.phys_ddist,  # Direction distribution (isotropic)
            self.target_area_dist,  # Area distribution (fixed target area)
            ]

        # Setup injector
        self.injector = siren.injection.Injector()
        self.injector.number_of_events = self.events_to_inject
        self.injector.detector_model = self.detector_model
        self.injector.primary_type = self.primary_type
        self.injector.primary_interactions = self.primary_processes[self.primary_type]
        self.injector.primary_injection_distributions = self.primary_injection_distributions
        self.injector._Injector__initialize_injector()

        # Weighting (from DIS_ATLAS example from SIREN)
        self.weighter = siren.injection.Weighter()
        self.weighter.injectors = [self.injector]
        self.weighter.detector_model = self.detector_model
        self.weighter.primary_type = self.primary_type
        self.weighter.primary_interactions = self.primary_processes[self.primary_type]
        self.weighter.primary_physical_distributions = self.primary_physical_distributions
        self.weighter.secondary_interactions = {}
        self.weighter.secondary_physical_distributions = {}


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

        event = self.injector.generate_event()

        weight = self.weighter(event)

        #Create tree
        tree = dataclasses.I3MCTree()
        particle = dataclasses.I3Particle() #Create a particle object
        interaction_record = siren.dataclasses.InteractionRecord()

        primaries = []

        for datum in event.tree:
            if datum.parent is None:
                record = datum.record
                primaries.append(record)

        secondaries = []

        for prim in primaries:
            record = prim
            #Extract the information from the SIREN interaction record and pass them to the I3 particle
            major_ID = record.primary_id.major_id
            minor_ID = record.primary_id.minor_id
            p_type = record.signature.primary_type
            prim_position = record.primary_initial_position
            vertex = record.interaction_vertex
            length = np.sqrt(np.sum((np.array(vertex)- np.array(prim_position))**2))

            px = record.primary_momentum[1]
            py = record.primary_momentum[2]
            pz = record.primary_momentum[3]
            momentum_magnitude = np.sqrt(px**2 + py**2 + pz**2)
            dir_x = px/momentum_magnitude
            dir_y = py/momentum_magnitude
            dir_z = pz/momentum_magnitude
            kinetic_energy = momentum_magnitude
            if record.primary_mass == 0:
                beta = 1
            else:
                gamma = record.primary_momentum[0]/record.primary_mass
                beta = np.sqrt(1 - 1/gamma**2)
            speed = beta * dataclasses.I3Constants.c
            time = 0

            #And now pass the information to the I3Particle
            particle.type = dataclasses.I3Particle.ParticleType(int(p_type))
            particle.id.majorID = major_ID
            particle.id.minorID = minor_ID
            particle.pos = dataclasses.I3Position(*prim_position)
            particle.dir = dataclasses.I3Direction(dir_x, dir_y, dir_z)
            particle.time = time #Time for the neutrino in the W target
            particle.energy = kinetic_energy
            particle.length = length
            particle.speed = speed

            tree.add_primary(particle)

            time = time + length/speed

            for i in range(len(record.secondary_ids)):
                secondaries.append([record.primary_id, particle.id, record, time, i])

        while len(secondaries) > 0:
            s_parent_id, parent_id, parent_record, parent_time, index = secondaries[0]
            secondaries.pop(0)

            particle = dataclasses.I3Particle()
            s_type = parent_record.signature.secondary_types
            particle.type = dataclasses.I3Particle.ParticleType(int(s_type[index]))
            major_ID = parent_record.secondary_ids[index].major_id
            minor_ID = parent_record.secondary_ids[index].minor_id
            sec_position = parent_record.interaction_vertex
            if record.secondary_masses[index] == 0:
                beta = 1
            else:
                gamma = record.secondary_momenta[index][0]/record.secondary_masses[index]
                beta = np.sqrt(1 - 1/gamma**2)
            speed = beta * dataclasses.I3Constants.c

            px = parent_record.secondary_momenta[index][1]
            py = parent_record.secondary_momenta[index][2]
            pz = parent_record.secondary_momenta[index][3]
            momentum_magnitude = np.sqrt(px**2 + py**2 + pz**2)
            dir_x = px/momentum_magnitude
            dir_y = py/momentum_magnitude
            dir_z = pz/momentum_magnitude
            kinetic_energy = momentum_magnitude

            particle.id.majorID = major_ID
            particle.id.minorID = minor_ID
            particle.pos = dataclasses.I3Position(*sec_position)
            particle.dir = dataclasses.I3Direction(dir_x, dir_y, dir_z)
            particle.time = parent_time
            particle.energy = kinetic_energy
            particle.speed = speed

            records = [datum.record for datum in event.tree if datum.record.primary_id == parent_record.secondary_ids[index]]
            if len(records) > 0:
                record = records[0]
                length = np.sqrt(np.sum((np.array(record.interaction_vertex) - np.array(sec_position))**2))
                time = parent_time + length/speed
                for i in range(len(record.secondary_ids)):
                    secondaries.append([record.primary_id, particle.id, record, time, i])
            else:
                length = np.nan

            particle.length = length

            tree.append_child(parent_id, particle)


        #Add the mctree to the frame
        frame["SIRENMarleyInjectionTree"] = tree

        #Add the weight to the frame
        frame["SIRENInjectionWeight"] = dataclasses.I3Double(weight)
        self.PushFrame(frame)

    def FillSimulationFrame(self, frame):
        # Fill the simulation frame with the event information
        frame["SIRENConfiguration"] = dataclasses.I3Int32(1) #SIREN configuration empty for now

