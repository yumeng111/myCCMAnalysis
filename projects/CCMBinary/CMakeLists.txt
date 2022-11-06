set(PROJECT_NAME CCMBinary)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME} INTERFACE
    ROOTCINT
    public/CCMAnalysis/CCMBinary/data_structures.hh)
