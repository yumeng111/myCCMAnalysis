ccm_project(CCMAnalysis)

ccm_executable(CCMAnalysis
    apps/CCMAnalysis.cxx
    USE_TOOLS root xercesc argagg
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures CCMIO CCMFramework CCMReco CCMNearline
    NO_PREFIX
)
