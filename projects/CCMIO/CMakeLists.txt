set(PROJECT_NAME CCMIO)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/CCMI3RootReader.cxx
    private/CCMRootIO.cxx
    private/CCMRootHandle.cxx
    private/CCMRawIO.cxx
    private/CCMEventTreeHandle.cxx
    private/IOUtils.cxx
    USE_TOOLS root boost python
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures icetray serialization dataio
    ROOTCINT
    ../icetray/public/icetray/I3FrameObject.h
    public/CCMAnalysis/CCMIO/CCMRootIO.h
    public/CCMAnalysis/CCMIO/CCMRawIO.h
    public/CCMAnalysis/CCMIO/CCMEventTreeHandle.h
    public/CCMAnalysis/CCMIO/IOUtils.h
)
