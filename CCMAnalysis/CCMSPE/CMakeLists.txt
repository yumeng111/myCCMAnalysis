set(PROJECT_NAME CCMSPE)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/SPECalibrationVariables.cxx
    private/MakeSPECalibrationHists.cxx
    private/NearlineSPEDiag.cxx
    USE_TOOLS root
    USE_PROJECTS CCMBinary CCMUtils CCMDataStructures
    ROOTCINT
    public/CCMAnalysis/CCMSPE/SPECalibrationVariables.h
    public/CCMAnalysis/CCMSPE/MakeSPECalibrationHists.h
    public/CCMAnalysis/CCMSPE/NearlineSPEDiag.h
)
