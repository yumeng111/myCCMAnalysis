set(PROJECT_NAME CCMSimulationUtils)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMPMTResponse.cxx
    private/CCMStrobeOverlay.cxx
    private/CCMTimeShift.cxx
    private/CCMQuenchingFactor.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMFramework CCMIO
    ROOTCINT
    public/CCMAnalysis/CCMSimulationUtils/CCMPMTResponse.h
    public/CCMAnalysis/CCMSimulationUtils/CCMStrobeOverlay.h
    public/CCMAnalysis/CCMSimulationUtils/CCMTimeShift.h
    public/CCMAnalysis/CCMSimulationUtils/CCMQuenchingFactor.h
)
