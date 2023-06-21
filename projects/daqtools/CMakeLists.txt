set(PROJECT_NAME daqtools)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMDAQAlignment.cxx
    private/CCMGeometryGenerator.cxx
    private/CCMTriggerMerger.cxx
    private/CCMBaselineAnalyzer.cxx
    private/CCMPeakFinder.cxx
    private/CCMNewBaselines.cxx
    private/CCMAverageBaselines.cxx
    private/CCMIDSPERegions.cxx
    private/CCMMichelElectrons.cxx
    private/CountCosmicFrames.cxx
    private/BeamCurrentMonitorSummary.cxx
    private/NIMLogicPulseFinder.cxx
    private/PulseCollector.cxx
    private/SumPulses.cxx
    private/FirstPulseFinder.cxx
    private/BaselineEstimator.cxx
    IWYU
    USE_TOOLS boost python nlopt
    USE_PROJECTS icetray serialization dataio dataclasses CCMBinary phys-services
)
