set(MODULENAME CCMSPE)

set(LOCAL_HEADERS
  SPECalibrationVariables.h
  MakeSPECalibrationHists.h
  NearlineSPEDiag.h
  )

set(LOCAL_SRCS
  SPECalibrationVariables.cxx
  MakeSPECalibrationHists.cxx
  NearlineSPEDiag.cxx
  )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

