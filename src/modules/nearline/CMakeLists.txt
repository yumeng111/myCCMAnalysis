set(MODULENAME CCMNearline)

set(LOCAL_HEADERS 
  CCMFindPulses.h
  CCMSPECalc.h
  CCMRateCalc.h
  )

set(LOCAL_SRCS
  CCMFindPulses.cxx
  CCMSPECalc.cxx
  CCMRateCalc.cxx
  )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

