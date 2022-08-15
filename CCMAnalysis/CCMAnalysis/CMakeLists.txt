ccm_project(CCMAnalysis)

ccm_executable(CCMAnalysis
    apps/CCMAnalysis.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMIO CCMFramework CCMReco CCMNearline
    NO_PREFIX
)
