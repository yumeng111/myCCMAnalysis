set(MODULENAME CCMUtils)

set(LOCAL_HEADERS
  MsgLog.h
  PMTInformation.h
  PMTInfoMap.h
  PMTInfoMap.tcc
  CCMBeamInfo.h
  Utility.h
)

set(LOCAL_SRCS
  MsgLog.cxx
  PMTInformation.cxx
  PMTInfoMap.cxx
  CCMBeamInfo.cxx
  Utility.cxx
)

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

