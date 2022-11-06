set(PROJECT_NAME CCMReco)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMFindEvents.cxx
    private/CCMNa22Cuts.cxx
    private/CCMBuildAccumWaveform.cxx
    private/CCMSingletTriplet.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMIO CCMFramework
    ROOTCINT
    public/CCMAnalysis/CCMReco/CCMFindEvents.h
    public/CCMAnalysis/CCMReco/CCMNa22Cuts.h
    public/CCMAnalysis/CCMReco/CCMBuildAccumWaveform.h
    public/CCMAnalysis/CCMReco/CCMSingletTriplet.h
)
