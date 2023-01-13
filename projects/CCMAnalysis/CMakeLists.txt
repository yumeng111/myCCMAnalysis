i3_project(CCMAnalysis)

i3_executable(CCMAnalysis
    apps/CCMAnalysis.cxx
    USE_TOOLS root xercesc argagg
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMIO CCMFramework CCMReco CCMNearline CCMPlots CCMSPE CCMSimulationUtils
    NO_PREFIX
)
