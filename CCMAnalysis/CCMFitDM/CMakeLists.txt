ccm_project(CCMFitDM)

ccm_executable(CCMFitDM
    apps/CCMFitDM.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMIO CCMFramework
    NO_PREFIX
)
