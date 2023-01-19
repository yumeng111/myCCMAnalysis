i3_project(CCMFitDM)

i3_executable(CCMFitDM
    apps/CCMFitDM.cxx
    USE_TOOLS root xercesc
    USE_PROJECTS CCMBinary CCMUtils CCMIO CCMFramework icetray serialization
    NO_PREFIX
)
