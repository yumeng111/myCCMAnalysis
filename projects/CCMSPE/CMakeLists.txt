set(PROJECT_NAME CCMSPE)

i3_project(${PROJECT_NAME})

i3_add_library(${PROJECT_NAME}
    private/SPECalibrationVariables.cxx
    private/MakeSPECalibrationHists.cxx
    private/NearlineSPEDiag.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures icetray serialization
    ROOTCINT
    public/CCMAnalysis/CCMSPE/SPECalibrationVariables.h
    public/CCMAnalysis/CCMSPE/MakeSPECalibrationHists.h
    public/CCMAnalysis/CCMSPE/NearlineSPEDiag.h
)
