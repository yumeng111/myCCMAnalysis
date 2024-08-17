set(PROJECT_NAME daqtools)

i3_project(${PROJECT_NAME}
    DOCS_DIR resources/docs
    PYTHON_DIR python
    PYTHON_DEST icecube/${PROJECT_NAME}
)

i3_add_library(${PROJECT_NAME}
    private/daqtools/CCMDAQAlignment.cxx
    private/daqtools/CCMGeometryGenerator.cxx
    private/daqtools/CCMTriggerMerger.cxx
    private/daqtools/CCMPeakFinder.cxx
    private/daqtools/CCMIDSPERegions.cxx
    private/daqtools/CCMMichelElectrons.cxx
    private/daqtools/ElectronicsCorrection.cxx
    private/daqtools/CountCosmicFrames.cxx
    private/daqtools/TestWaveformAccum.cxx
    private/daqtools/FitSummedPulses.cxx
    private/daqtools/BaselineEstimator.cxx
    private/daqtools/SumWaveforms.cxx
    private/daqtools/BeamCurrentMonitorSummary.cxx
    private/daqtools/FlightPathThreeSummary.cxx
    private/daqtools/NIMLogicPulseFinder.cxx
    private/daqtools/PulseCollector.cxx
    private/daqtools/SumPulses.cxx
    private/daqtools/FirstPulseFinder.cxx
    private/daqtools/CCMDerivativePulseFinder.cxx
    private/daqtools/AccumulateSummedWaveforms.cxx
    private/daqtools/AccumulateIndividualPMTWaveforms.cxx
    private/daqtools/AccumulateIndividualChannelWaveforms.cxx
    private/daqtools/SumPulsesSPEShapeCalibration.cxx
    IWYU
    USE_TOOLS boost python nlopt
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services
)

i3_add_pybindings(${PROJECT_NAME}
    private/pybindings/module.cxx
    private/pybindings/*.cxx
    USE_TOOLS boost python
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services daqtools
)
