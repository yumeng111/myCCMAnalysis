set(PROJECT_NAME CCMUtils)

ccm_project(${PROJECT_NAME})

ccm_add_library(${PROJECT_NAME}
    private/MsgLog.cxx
    private/PMTInformation.cxx
    private/PMTInfoMap.cxx
    private/CCMBeamInfo.cxx
    private/Utility.cxx
    USE_TOOLS ROOT
    USE_PROJECTS CCMBinary
    ROOTCINT
    public/CCMAnalysis/CCMUtils/MsgLog.h
    public/CCMAnalysis/CCMUtils/PMTInformation.h
    public/CCMAnalysis/CCMUtils/PMTInfoMap.h
    public/CCMAnalysis/CCMUtils/PMTInfoMap.tcc
    public/CCMAnalysis/CCMUtils/CCMBeamInfo.h
    public/CCMAnalysis/CCMUtils/Utility.h
)
