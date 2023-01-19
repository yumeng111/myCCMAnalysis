set(PROJECT_NAME CCMIO)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMRootIO.cxx
    private/CCMRawIO.cxx
    private/CCMEventTreeHandle.cxx
    private/IOUtils.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures icetray serialization
    ROOTCINT
    public/CCMAnalysis/CCMIO/CCMRootIO.h
    public/CCMAnalysis/CCMIO/CCMRawIO.h
    public/CCMAnalysis/CCMIO/CCMEventTreeHandle.h
    public/CCMAnalysis/CCMIO/IOUtils.h
)
