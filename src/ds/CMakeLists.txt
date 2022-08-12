set(MODULENAME CCMDS)

set(LOCAL_HEADERS
  data_structures.hh
  Pulses.h
  RawData.h
  SimplifiedEvent.h
  SinglePulse.h
  Events.h
  AccumWaveform.h
  MCTruth.h
  )

set(LOCAL_SRCS
  Pulses.cxx
  RawData.cxx
  SimplifiedEvent.cxx
  SinglePulse.cxx
  Events.cxx
  AccumWaveform.cxx
  MCTruth.cxx
  )

ccm_module("${MODULENAME}" "${LOCAL_HEADERS}" "${LOCAL_SOURCES}")

