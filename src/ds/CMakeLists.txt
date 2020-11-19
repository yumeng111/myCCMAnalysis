set(DS_HEADERS
  data_structures.hh
  Pulses.h
  RawData.h
  SimplifiedEvent.h
  SinglePulse.h
  Events.h
  AccumWaveform.h
	MCTruth.h
  )

set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR} ${CCMAnalysis_INCLUDE_DIRS})

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

link_directories(${LINK_DIRECTORIES})

set(CMAKE_CXX_LINKER_FLAGS "${CMAKE_CXX_LINKER_FLAGS} ${ROOT_LIBRARIES}")

set(DS_SRCS
  Pulses.cxx
  RawData.cxx
  SimplifiedEvent.cxx
  SinglePulse.cxx
  Events.cxx
  AccumWaveform.cxx
	MCTruth.cxx
  )

set(DS_LINKDEF LinkDef.h)
set(DS_DICTIONARY CCMDSDict)
ROOT_GENERATE_DICTIONARY("${DS_DICTIONARY}" "${DS_HEADERS}" LINKDEF "${DS_LINKDEF}")

set(DS_SRCS ${DS_SRCS} ${DS_DICTIONARY}.cxx)

add_library(CCMDS SHARED ${DS_SRCS})
#CREATE_LIBTOOL_FILE(CCMDS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
target_link_libraries(CCMDS ${ROOT_LIBRARIES} CCMUtils)

install(FILES ${DS_HEADERS} DESTINATION $ENV{CCMINSTALLBuild}/include)
