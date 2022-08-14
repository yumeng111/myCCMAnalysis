set(PROJECT_NAME utils)

ccm_project(${PROJECT_NAME})

ccm_add_library(utils
    private/MsgLog.cxx
    private/PMTInformation.cxx
    private/PMTInfoMap.cxx
    private/CCMBeamInfo.cxx
    private/Utility.cxx
    USE_TOOLS ROOT
)

# set(LOCAL_HEADERS
#   MsgLog.h
#   PMTInformation.h
#   PMTInfoMap.h
#   PMTInfoMap.tcc
#   CCMBeamInfo.h
#   Utility.h
# )
# 
# set(LOCAL_SRCS
#   MsgLog.cxx
#   PMTInformation.cxx
#   PMTInfoMap.cxx
#   CCMBeamInfo.cxx
#   Utility.cxx
# )
# 
# ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

