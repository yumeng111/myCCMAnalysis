set(MODULENAME CCMReco)

set(LOCAL_HEADERS
  CCMFindEvents.h
  CCMNa22Cuts.h
  CCMBuildAccumWaveform.h
  CCMSingletTriplet.h
  )

set(LOCAL_SRCS
  CCMFindEvents.cxx
  CCMNa22Cuts.cxx
  CCMBuildAccumWaveform.cxx
  CCMSingletTriplet.cxx
  )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

