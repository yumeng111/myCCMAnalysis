#!/usr/bin/env python3
import I3Tray
import os
import argparse
import uuid
import logging #
import sys
import json
import numpy as np 
from icecube import icetray, dataclasses, dataio, phys_services, hdfwriter, daqtools
from icecube.icetray import load
from icecube.icetray import I3Frame, I3Module, I3ConditionalModule
from icecube.icetray import I3Units
import monoMuons as injectorMuons
from icecube.icetray import CCMTriggerKey, CCMPMTKey
from icecube.dataclasses import CCMOMGeo, I3VectorDouble, I3Double

load("CCMBinary", False)
load("daqtools", False)
load("dataclasses", False)
load("nearline-reco", False)
load("g4-larsim")
load("simutils", False)

#icetray.set_log_level(icetray.I3LogLevel.LOG_DEBUG) 

def compute_total_deposited_energy_LAr(frame):
    tree = frame["InnerLArMCTreeVector"]
    total_deposited_energy_LAr = np.sum([p.energy for p in tree])
    total_deposited_energy_LAr = dataclasses.I3Double(total_deposited_energy_LAr)
    frame["TotalDepositedEnergyLAr"] = total_deposited_energy_LAr

def compute_total_muon_length_LAr(frame):
    tree = frame["InnerLArMCTreeVector"]
    total_muon_length_LAr = np.sum([q.length for q in tree if (q.type == 13 or q.type ==-13)])
    total_muon_length_LAr = dataclasses.I3Double(total_muon_length_LAr)
    frame["TotalMuonLengthLAr"] = total_muon_length_LAr

if __name__ == "__main__":
    return_code = 0

    parser = argparse.ArgumentParser()
    parser.add_argument("--n", default=10, help="Number of events to simulate", type=int)
    parser.add_argument("--file_name", default = "Muons", help = "name of file")
    args=parser.parse_args()

    job_file_name = str(args.file_name)
    # Config logging
    log_filename = "MuonsG4.log"
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
        logging.FileHandler(log_filename), # log to file
        logging.StreamHandler(sys.stdout), # log to console
        ],
    )

    logger = logging.getLogger()
    logger.info(f"Starting Injector+G4LAr simulation for {args.n} events")

    # Add our random service
    seed = int(int(uuid.uuid4()) % 1e16)
    print(f"seed = {seed}")
    logger.info(f"Using random seed: {seed}")
    randomService = phys_services.I3GSLRandomService(seed)

    tray = I3Tray.I3Tray()
    tray.context["I3RandomService"] = randomService
    #Empty DAQ frames
    tray.Add("I3InfiniteSource", Prefix="sim_gcd_combined.i3.zst")

    tray.Add(injectorMuons.MuonInjector)
    tray.Add("Dump")

    tray.AddService(
        "CCM200ResponseFactory",
        "response",
        SaveAllEnergyLossesTree = False, # Save ALL energy losses to a tree "LArMCTree"
        VetoSDSaveEnergyLossesVector = False, # Save veto energy losses to a vector "VetoLArMCTreeVector"
        VetoSDSaveEnergyLossesTree = True, # Save veto energy losses to a tree "VetoLArMCTree"
        VetoSDPruneTree = False, # Remove particles from the veto tree that do not lead to an energy loss in the veto

        InteriorSDSaveEnergyLossesVector = True, # Save interior energy losses to a vector "InnerLArMCTreeVector"
        InteriorSDSaveEnergyLossesTree = False, # Save interior energy losses to a tree "InnerLArMCTree"
        InteriorSDPruneTree = False, # Remove particles from the interior tree that do not lead to an energy loss in the interior
   
        KillNeutrinos = False, # Kill neutrinos upon first interaction
        KillPhotons = False, # Kill photons upon first interaction
        KillScintillation = False, # Kill scintillation photons upon first interaction
        KillCherenkov = False, # Kill Cherenkov photons upon first interaction

        TimeCut = False, # Cut all hits that are not in a 200ns time window (this is fine to do if modeling scintillation timing via the PMT response module)
        DetailedPhotonTracking = False, # Track photon history in the simulation (this is not recommended)
        TrackParticles = False, # Track all particles in the simulation (this is not needed for production simulation but can be useful)
        TrackEnergyLosses = False, # Track energy losses in the simulation (this is not needed for production simulation but can be useful)

        SimulateNuclearRecoils = False, # Sets the production threshold of protons to zero
        G4RangeCut = 0.7 * I3Units.mm, # Range cut for all particles
        G4EDepMin = 0.0 * I3Units.keV, # The minimum energy deposit to be recorded (sub-threshold energy deposits are aggregated, but the energy loss type is no longer valid)
        G4ETrackingMin = 0.0 * I3Units.keV, # The minimum kinetic energy a particle needs to have to be tracked (energy losses of sub-threshold particles are transferred to their parent)

        RecordHits = False, # Record photon hits on PMTs (this is needed for production simulation)
        SourceRodIn = False, # Set to true if the source rod is inside the detector
        SourceRodLocation = 0.0 * I3Units.cm, # The location of the source rod in the detector, 0cm is the center of the detector

        CobaltSourceRun = False, # Set to true if the source is a cobalt source
        SodiumSourceRun = False, # Set to true if the source is a sodium source
        TrainingSource = False, # Set to true if the source is a training source (this is not the "blank source")
        DecayX = 0.0 * I3Units.cm, # The x position of the source in the detector (TrainingSource must be true)
        DecayY = 0.0 * I3Units.cm, # The y position of the source in the detector (TrainingSource must be true)
        DecayZ = 0.0 * I3Units.cm, # The z position of the source in the detector (TrainingSource must be true)

        Rayleigh128Length = 99.1 * I3Units.cm, # Rayleigh scattering length for 128nm photons
        #EnableUVAbsorption = False,
        UVAbsA = 0.0618 , # f_(Transmission) = 1- exp(-(a * lambda - b)) [1/nm]
        UVAbsB = 113.0, # [nm]
        UVAbsD = 5.8, # abs_length = -d / log(f_T) [cm]
        UVAbsScaling = 3.154, # abs_length *= UVAbsScaling

        # Poisson mean number of photons produced per wavelength shift in different regions of TPB
        WLSNPhotonsEndCapFoil = 0.605,
        WLSNPhotonsSideFoil = 0.605,
        WLSNPhotonsPMT = 0.605,

        EndCapFoilTPBThickness = 0.00278035 * I3Units.mm,
        SideFoilTPBThickness = 0.00278035 * I3Units.mm,
        PMTTPBThickness = 0.00203892 * I3Units.mm,

        # TPB absorption is defined by a table, below ~380 nm
        # L = Scale * Table(lambda)
        # Above this value it is defined by an exponential distribution
        # L = Scale * Norm * exp(lambda / Tau)
        TPBAbsorptionTau = 0.13457,
        TPBAbsorptionNorm = 8.13914e-21,
        TPBAbsorptionScale = 1.0,

        MieGG = 0.84, # Forward and backward anisotropy of the Mie scattering in the TPB (MIEHG_FORWARD, MIEHG_BACKWARD)
        MieRatio = 0.90, # Ratio of the Mie scattering in the TPB (MIEHG_FORWARD_RATIO)

        Normalization = 0.778, # Normalization factor for scintillation light production (relative to 4e4 photons/MeV)
        PhotonSampling = 0.0, # Sampling factor for photons (below 1 is under sampling, above 1 is over sampling). This is corrected in the PMT response module, but can improve speed for high energy events
        RandomSeed=seed,
        )

    tray.AddModule(
        "CCMSimulator",
        "ccm simulator",
        ResponseServiceName="response",
        #ConfigurationName="DetectorResponseConfig",
        InputMCTreeName="I3MCTree",
        PMTHitSeriesName="PMTMCHitsMap",
        LArMCTreeName="LArMCTree",
        #PhotonSummarySeriesName="PhotonSummarySeries",
        Multithreaded = True,
        MultiParticleOutput = True,
        PerEventOutput = False,
        NumberOfThreads=2,
        BatchSize = 1
        )
    tray.Add("Dump")

    tray.Add(
        "PMTResponse",
        InputHitsMapName="PMTMCHitsMapMultiParticle",
        OutputRecoPulseName="MCRecoPulses",
        OutputTimeOffsetsName="SimulatedBoardTimeOffsets",
        OutputTrueEventTimeName="TrueEventTime",
        DetectorConfigurationName="DetectorResponseConfig",
        RandomServiceName="I3RandomService", # Use default random service
        #QEWavelengths=[], # Reasonable default values are hardcoded in the module
        #QEValues=[],
        RemoveCherenkov=False, # Remove Cherenkov photons
        FlatEfficiency=False, # Use flat efficiency
        WeightUVAbsorption=True, # UV absorption
    )
    tray.Add(compute_total_deposited_energy_LAr, Streams = [icetray.I3Frame.DAQ])
    #tray.Add(compute_total_charge, streams = [icetray.I3Frame.DAQ]) 
    tray.Add(compute_total_muon_length_LAr, Streams = [icetray.I3Frame.DAQ])

    charge_threshold = 10.0
    time_window = 10.0 * I3Units.ns
    event_prefix = f"RecoQ{int(charge_threshold):d}T{int(time_window):d}"
    event_stream_name = f"EventFinder{int(charge_threshold):d}PE{int(time_window):d}ns"

    tray.Add("Dump")
    
    tray.Add(
    "EventFinder",
    EventChargeThreshold=charge_threshold, #This is only one value (10PE in this case)
    TimeWindow=time_window,
    PMTTypes = [dataclasses.CCMOMGeo.OMType.CCM8inUncoated, dataclasses.CCMOMGeo.OMType.CCM8inCoated],
    BadKeysNames = ["KeysExcludedFromFit", "KeysIncompleteSPE", "KeysHighNoise"],
    Output=event_prefix,
    Pulses="MCRecoPulses",
    )
    
    tray.Add(
    "CCMEventSplitter",
    InputPrefix=event_prefix,
    OutputPrefix=event_prefix,
    SubEventStreamName=event_stream_name,
    #If=has_pulses_for_splitter #I added this, it was exploding
    )

    tray.Add(
        hdfwriter.CCMSimHDFWriter,
        Output="test.h5",
        Keys=[
           "CCMEventHeader", "TotalDepositedEnergyLAr", "I3MCTree", "TotalMuonLengthLAr"]
    )

    tray.AddModule("I3Writer", 
        "writer",
        FileName ="test.i3.zst",
        Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Simulation, icetray.I3Frame.Physics])

    try:
        tray.Execute(args.n + 2)
        logger.info(f"Muons G4 simulation finished successfully for {args.n} events")
    except KeyboardInterrupt as e:
        return_code = 2 ### unprocessed
    except Exception as e:
        logger.error(f"Error during execution: {e}", exc_info=True)
        return_code = 1 ### failure
    else:
        return_code = 0 ### processed
    exit(return_code)

