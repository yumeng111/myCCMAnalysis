set(PROJECT_NAME CCMPlots)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMSumWaveforms.cxx
    private/CCMProcessingPlots.cxx
    private/CCMNuDMPlots.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMIO CCMFramework
    ROOTCINT
    public/CCMAnalysis/CCMPlots/CCMSumWaveforms.h
    public/CCMAnalysis/CCMPlots/CCMProcessingPlots.h
    public/CCMAnalysis/CCMPlots/CCMNuDMPlots.h
)
