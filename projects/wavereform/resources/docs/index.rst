.. 
.. copyright  (C) 2010
.. The Icecube Collaboration
.. 
.. $Id: index.rst 66546 2010-08-31 18:09:44Z kislat $
.. 
.. @version $Revision: 66546 $
.. @date $LastChangedDate: 2010-08-31 13:09:44 -0500 (Tue, 31 Aug 2010) $
.. @author Jakob van Santen <vansanten@wisc.edu> $LastChangedBy: kislat $

.. highlight:: python

.. _wavereform-main:

wavereform
================

The wavereform project contains tools for convolving I3RecoPulses with a pulse template to reconstruct the
waveform they represent and checking this reconstruction against the input waveform with a simple
:math:`\chi^2`-like criterion. This criterion can be used to identify failures in calibration and feature
extraction.

.. toctree::
   :maxdepth: 1
   
   release_notes


Modules
___________


I3Wavereform
--------------------

The I3Wavereform module takes an I3WaveformSeriesMap and an I3RecoPulseSeries map as input. For each
I3Waveform, it convolves the corresponding I3RecoPulseSeries with the appropriate pulse template to obtain
the waveform the pulses represent and calculates a simple :math:`\chi^2`-like quantity between the
reconstructed and measured waveforms. This is defined as

.. math::
    
    \chi^2 \equiv \sum_i \left\{\frac{{\rm waveform}_i - {\rm refolded}_i}{{\rm step}_i} \right\}^2

where the sum goes over all digitizer bins and the step is the digitizer's voltage discretization step.

Optionally the module can write out a list of OMKeys where the :math:`\chi^2` for at least one waveform
exceeded a configured threshold. This can be used to flag those readouts for special handling (e.g.
transmission over the satellite in raw rather than compressed form).

Parameters
^^^^^^^^^^

**Waveforms**:

  Name of the I3WaveformSeriesMap to get from the frame, e.g. "CalibratedATWD".
  
**RecoPulses**:

  Name of the I3RecoPulseSeriesMap to get from the frame, e.g. "DeformedPulses".
  
**Chi**:

  Name of the I3MapOMKeyVectorDouble to put in the frame, e.g. "ATWDChi_DeformedPulses".
  
**ChiThreshold**:

  Threshold at which to mark a waveform as badly handled.
  
**Flag**:

  Name of the I3VectorOMKey to put in the frame containing those keys for which at least one waveform
  was marked as badly handled.

OutBoxes
^^^^^^^^
One.

Example
^^^^^^^^

::
    
    from icecube import icetray, dataio, WaveCalibrator
    import I3Tray
    
    tray = I3Tray.I3Tray()
    
    tray.AddModule("I3Reader","reader",
        Filename="data.i3.gz"
        )
    
    tray.AddModule("I3WaveCalibrator", "sedan",
        Launches="InIceRawData",
        Waveforms="CalibratedWaveforms",
        Errata="BorkedOMs",
        ATWDMode=WaveCalibrator.CalibrationMode.CALIBRATE_UNSATURATED,
        ATWDSaturationMargin=123,
        FADCSaturationMargin=0,
        )
        
    # Extract pulses.
    tray.AddModule('I3Wavedeform', 'deform')

    # Perform some sanity checks.
    tray.AddModule("I3Wavereform", "wavereform_atwd",
    	Waveforms="CalibratedWaveforms",
    	Pulses="WavedeformPulses",
    	Chi="Chi_HLCPulses_CalibratedATWD",
    	ChiThreshold=1e4, # flag ~ 0.2% of waveforms
    	Flag="Borked_ATWDs",
    	)
    
    
    tray.Execute()
    


I3LaunchSelector
--------------------

The I3LaunchSelector module takes a number of I3VectorOMKeys and an I3DOMLaunchSeriesMap from the frame
and outputs an I3DOMLaunchSeriesMap containing only entries for the keys that appear in one of the input
lists.

Parameters
^^^^^^^^^^

**Launches**:

  Name of the I3DOMLaunchSeriesMap to get from the frame, e.g. "InIceRawData".

**Flag**:

  List of names of I3VectorOMKeys to get from the frame.

**Output**:

  Name of the I3DOMLaunchSeriesMap to put in the frame, e.g. "InIceRawDataErrata".

OutBoxes
^^^^^^^^
One.

Example
^^^^^^^^

::

    from icecube import icetray, dataio, WaveCalibrator
    import I3Tray

    tray = I3Tray.I3Tray()

    tray.AddModule("I3Reader","reader",
        Filename="data.i3.gz"
        )

    tray.AddModule("I3WaveCalibrator", "sedan",
        Launches="InIceRawData",
        Waveforms="CalibratedWaveforms",
        Errata="BorkedOMs",
        ATWDMode=WaveCalibrator.CalibrationMode.CALIBRATE_UNSATURATED,
        ATWDSaturationMargin=123,
        FADCSaturationMargin=0,
        )
        
    # Extract pulses.
    tray.AddModule('I3Wavedeform', 'deform')

    # Perform some sanity checks.
    tray.AddModule("I3Wavereform", "wavereform_atwd",
    	Waveforms="CalibratedWaveforms",
    	Pulses="WavedeformPulses",
    	Chi="Chi_HLCPulses_CalibratedATWD",
    	ChiThreshold=1e4, # flag ~ 0.2% of waveforms
    	Flag="Borked_ATWDs",
    	)
    	
    # Save raw data for those OMs where calibration or feature extraction
    # failed in an obvious way. Two keys are added to the frame:
    # InIceErrata => raw data for screwy DOMs (keep only this if discarding the DAQ payload)
    # InIceErrataKeys => compact list of screwy keys (keep only this if transmitting the DAQ payload) 
    tray.AddModule("I3LaunchSelector", "seatbelt",
    	Launches="InIceRawData",
    	Flags=["CalibrationErrata", "Borked_ATWDs"],
    	Output="InIceErrata",
    	)

    
    tray.Execute()
    


