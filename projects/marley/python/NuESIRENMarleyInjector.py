import siren
from siren import utilities
import numpy as np
import traceback
from icecube import icetray, dataclasses
from icecube.icetray import I3ConditionalModule, logging

if siren.dataclasses.ParticleID.__hash__ is None:

    def pid_hash(self):
        return hash((self.is_set(), self.major_id, self.minor_id))

    siren.dataclasses.ParticleID.__hash__ = pid_hash


def siren_primary_to_i3_particle(record):
    particle = dataclasses.I3Particle()  # Create a particle object
    # Extract the information from the SIREN interaction record and pass them to the I3 particle
    major_ID = record.primary_id.major_id
    minor_ID = record.primary_id.minor_id
    p_type = record.signature.primary_type
    position = record.primary_initial_position
    vertex = record.interaction_vertex
    length = np.sqrt(np.sum((np.array(vertex) - np.array(position)) ** 2))

    momentum = np.array(record.primary_momentum[1:])
    momentum_magnitude = np.sqrt(np.sum(momentum**2))
    direction = momentum / momentum_magnitude
    mass = record.primary_mass
    energy = record.primary_momentum[0]
    kinetic_energy = energy - mass
    if mass == 0:
        beta = 1
    else:
        gamma = energy / mass
        beta = np.sqrt(1 - 1 / gamma**2)
    speed = beta * dataclasses.I3Constants.c
    time = 0

    # And now pass the information to the I3Particle
    particle.type = dataclasses.I3Particle.ParticleType(int(p_type))
    particle.id.majorID = major_ID
    particle.id.minorID = minor_ID
    particle.pos = dataclasses.I3Position(*position)
    particle.dir = dataclasses.I3Direction(*direction)
    particle.time = time  # Time for the neutrino in the W target
    particle.energy = kinetic_energy
    particle.length = length
    particle.speed = speed

    return particle


def siren_secondary_to_i3_particle(primary_particle, record, index):
    particle = dataclasses.I3Particle()

    secondary_type = record.signature.secondary_types[index]
    particle.type = dataclasses.I3Particle.ParticleType(int(secondary_type))
    major_ID = record.secondary_ids[index].major_id
    minor_ID = record.secondary_ids[index].minor_id
    position = record.interaction_vertex
    mass = record.secondary_masses[index]
    energy = record.secondary_momenta[index][0]
    if mass:
        beta = 1
    else:
        gamma = energy / mass
        beta = np.sqrt(1 - 1 / gamma**2)
    speed = beta * dataclasses.I3Constants.c
    time = primary_particle.time + primary_particle.length / primary_particle.speed

    momentum = np.array(record.secondary_momenta[index][1:])
    momentum_magnitude = np.sqrt(np.sum(momentum**2))

    direction = momentum / momentum_magnitude
    kinetic_energy = energy - mass

    particle.id.majorID = major_ID
    particle.id.minorID = minor_ID
    particle.pos = dataclasses.I3Position(*position)
    particle.dir = dataclasses.I3Direction(*direction)
    particle.time = time
    particle.energy = kinetic_energy
    particle.length = np.nan
    particle.speed = speed

    return particle


class NuESIRENMarleyInjector(I3ConditionalModule):
    """
    A module that injects NuE CC events into CCM using SIREN and the total cross section from MARLEY
    """

    def __init__(self, context):
        I3ConditionalModule.__init__(self, context)
        self.AddParameter(
            "Target", "0 for upper tungsten target, 1 for lower tungsten target", 0
        )
        self.AddParameter("EventsToInject", "Number of events to inject", 10000)
        self.seen_s_frame_ = False

    def Configure(self):
        self.target = self.GetParameter("Target")
        self.events_to_inject = self.GetParameter("EventsToInject")
        if self.target not in [0, 1]:
            raise ValueError("Target option must be 0 or 1")
        if self.events_to_inject <= 0:
            raise ValueError("EventsToInject must be a positive integer")

        # Load experiment and detector
        self.experiment = "CCM"
        self.detector_model = utilities.load_detector(self.experiment)
        self.fiducial_volume = utilities.get_fiducial_volume(self.experiment)

        # Define primary particle
        self.primary_type = siren.dataclasses.Particle.ParticleType.NuE

        # Load Marley process
        try:
            self.primary_processes, self.secondary_processes = (
                siren.resources.load_processes(
                    "MarleyCrossSection",
                    primary_types=[self.primary_type],
                    process_types=["CC"],
                )
            )
        except Exception as e:
            logging.log_info(
                "Caught exception when loading MarleyCrossSection siren process"
            )
            logging.log_info(traceback.format_exc())
            logging.log_fatal(str(e))
            raise

        # Primary distributions
        self.mass_ddist = siren.distributions.PrimaryMass(0)  # mass
        self.edist = siren.distributions.PiDARNuEDistribution()  # energy

        # Cone direction distribution from tungsten target
        # self.opening_angle = np.arcsin(1.21 / 23.0)
        self.opening_angle = np.arcsin(2.0 / 23.0)

        # Fixed tungsten target
        if self.target == 0:
            self.target_origin = siren.math.Vector3D(0, 0, 0.1375)
        elif self.target == 1:
            self.target_origin = siren.math.Vector3D(0, 0, -0.241)

        self.detector_origin = siren.math.Vector3D(23, 0, -0.65)
        self.target_dir = self.detector_origin - self.target_origin
        self.target_dir.normalize()
        self.target_inj_ddist = siren.distributions.Cone(
            self.target_dir, self.opening_angle
        )

        # Physical direction distribution
        self.phys_ddist = siren.distributions.IsotropicDirection()

        # Position distribution
        self.max_dist = 25.0
        # Position distribution: use fixed target distribution that considers physical volume of the target
        if self.target == 0:
            # upper target
            self.target_cylinder = siren.geometry.Cylinder(
                siren.geometry.Placement(self.target_origin - self.detector_origin),
                0.05,
                0.0,
                0.091,
            )
        else:
            # lower target
            self.target_cylinder = siren.geometry.Cylinder(
                siren.geometry.Placement(self.target_origin - self.detector_origin),
                0.05,
                0.0,
                0.298,
            )

        self.target_pos_dist = siren.distributions.FixedTargetPositionDistribution(
            self.target_cylinder, self.fiducial_volume, self.max_dist
        )

        self.target_area_dist = siren.distributions.FixedTargetAreaDistribution(
            self.target_cylinder
        )

        # Primary injection distributions
        self.primary_injection_distributions = [
            self.mass_ddist,  # Mass distribution
            self.edist,  # Energy distribution
            self.target_inj_ddist,  # Direction distribution (cone)
            self.target_pos_dist,  # Position distribution (fixed target position)
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
        self.injector.primary_injection_distributions = (
            self.primary_injection_distributions
        )
        self.injector._Injector__initialize_injector()

        # Weighting (from DIS_ATLAS example from SIREN)
        self.weighter = siren.injection.Weighter()
        self.weighter.injectors = [self.injector]
        self.weighter.detector_model = self.detector_model
        self.weighter.primary_type = self.primary_type
        self.weighter.primary_interactions = self.primary_processes[self.primary_type]
        self.weighter.primary_physical_distributions = (
            self.primary_physical_distributions
        )
        self.weighter.secondary_interactions = {}
        self.weighter.secondary_physical_distributions = {}

    def ProcessSimulationFrame(self, frame):
        self.seen_s_frame_ = True
        self.FillSimulationFrame(frame)
        self.PushFrame(frame)

    def Simulation(self, frame):
        self.ProcessSimulationFrame(frame)

    def DAQ(self, frame):

        if not self.seen_s_frame_:
            sim_frame = icetray.I3Frame(icetray.I3Frame.Simulation)
            self.ProcessSimulationFrame(sim_frame)

        event = self.injector.generate_event()
        weight = self.weighter(event)

        all_records = {datum.record.primary_id: datum.record for datum in event.tree}
        primaries = [datum.record for datum in event.tree if datum.parent is None]

        # Create tree
        tree = dataclasses.I3MCTree()

        secondaries = []

        for record in primaries:
            particle = siren_primary_to_i3_particle(record)
            tree.add_primary(particle)
            secondaries.extend(
                [(particle, record, i) for i in range(len(record.secondary_ids))]
            )

        while len(secondaries) > 0:
            primary_particle, record, index = secondaries[0]
            secondary_particle = siren_secondary_to_i3_particle(*secondaries[0])
            secondaries.pop(0)

            if record.secondary_ids[index] in all_records:
                secondary_record = all_records[record.secondary_ids[index]]
                position = secondary_record.primary_initial_position
                vertex = secondary_record.interaction_vertex
                length = np.sqrt(np.sum((np.array(vertex) - np.array(position)) ** 2))
                secondary_particle.length = length
                secondaries.extend(
                    [
                        (secondary_particle, secondary_record, i)
                        for i in range(len(secondary_record.secondary_ids))
                    ]
                )

            tree.append_child(primary_particle.id, secondary_particle)

        # Add the mctree to the frame
        frame["SIRENMarleyInjectionTree"] = tree

        # Add the weight to the frame
        frame["SIRENInjectionWeight"] = dataclasses.I3Double(weight)
        self.PushFrame(frame)

    def FillSimulationFrame(self, frame):
        # Fill the simulation frame with the event information
        frame["SIRENConfiguration"] = dataclasses.I3Int32(
            1
        )  # SIREN configuration empty for now
