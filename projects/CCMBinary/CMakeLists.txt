set(PROJECT_NAME CCMBinary)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME} INTERFACE
    ROOTCINT
    public/CCMAnalysis/CCMBinary/data_structures.hh)
