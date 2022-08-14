set(PROJECT_NAME CCMIO)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/CCMRootIO.cxx
    private/CCMRawIO.cxx
    private/CCMEventTreeHandle.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures
    ROOTCINT
    public/CCMAnalysis/CCMIO/CCMRootIO.h
    public/CCMAnalysis/CCMIO/CCMRawIO.h
    public/CCMAnalysis/CCMIO/CCMEventTreeHandle.h
)
