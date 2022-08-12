set(MODULENAME CCMDetSim)

set(LOCAL_HEADERS 
  CCMPMTResponse.h
  CCMStrobeOverlay.h
  CCMTimeShift.h
  CCMQuenchingFactor.h
)

set(LOCAL_SRCS
  CCMPMTResponse.cxx
  CCMStrobeOverlay.cxx
  CCMTimeShift.cxx
  CCMQuenchingFactor.cxx
)

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

