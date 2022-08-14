set(PROJECT_NAME CCMNearline)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/CCMFindPulses.cxx
    private/CCMSPECalc.cxx
    private/CCMRateCalc.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMSPE CCMDataStructures CCMFramework
    ROOTCINT
    public/CCMAnalysis/CCMNearline/CCMFindPulses.h
    public/CCMAnalysis/CCMNearline/CCMSPECalc.h
    public/CCMAnalysis/CCMNearline/CCMRateCalc.h
)
