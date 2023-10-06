set(PROJECT_NAME daqtools)

i3_project(${PROJECT_NAME}
    PYTHON_DIR python
    PYTHON_DEST icecube/${PROJECT_NAME}
)

i3_add_library(${PROJECT_NAME}
    private/CCMDAQAlignment.cxx
    private/CCMGeometryGenerator.cxx
    private/CCMTriggerMerger.cxx
    private/CCMPeakFinder.cxx
    private/CCMIDSPERegions.cxx
    private/CCMMichelElectrons.cxx
    private/ElectronicsCorrection.cxx
    private/CountCosmicFrames.cxx
    private/TestWaveformAccum.cxx
    private/FitSummedPulses.cxx
    private/BaselineEstimator.cxx
    private/SumWaveforms.cxx
    private/BeamCurrentMonitorSummary.cxx
    private/NIMLogicPulseFinder.cxx
    private/PulseCollector.cxx
    private/SumPulses.cxx
    private/FirstPulseFinder.cxx
    private/CCMDerivativePulseFinder.cxx
    private/EventFinder.cxx
    private/AccumulateWaveforms.cxx
    IWYU
    USE_TOOLS boost python nlopt
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services
)
